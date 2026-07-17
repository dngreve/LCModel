#!/usr/bin/env fspython
"""
lcm_convert.py

Run with `fspython` (must be on your PATH), not plain `python3` --
`surfa` and this project's other imaging dependencies live in the
`fspython` environment.

Pipeline connecting imaging volumes to LCModel and back, as three
subcommands (not two -- see note below on why `run` is separate). Input/
output volumes may be in **any format surfa supports** (nifti, FreeSurfer
.mgz, etc.) -- format is detected automatically from the file
extension/header by `surfa.load_volume`/`Volume.save`, exactly the same as
any other surfa-based tool; nothing in this file is nifti-specific, and
none of the CLI flags require a particular extension.

`to-lcm`   Read one (met OR ref) real+imaginary volume pair, mask it,
           optionally convert from MIDAS format to LCModel's time-domain
           FID format, optionally prepend a single reference voxel loaded
           from separate volume files, and write an LCModel .RAW-format
           file plus a sidecar .meta.json (needed by `from-csv` if a
           reference voxel was inserted).

`run`      Invoke `lcmodel` with this session's -i/-csv CLI overrides,
           given an existing met.lcm + h2o.lcm pair (each produced by its
           own separate `to-lcm --type met` / `to-lcm --type ref` call --
           `to-lcm` only ever produces one file type at a time, so running
           lcmodel, which needs BOTH files together via `-i met.lcm h2o.lcm`,
           is necessarily a separate step, not something `to-lcm` itself
           can do).

`from-csv` Read the .csv LCModel produced, map each row's 1-based "Col"
           field back to its voxel position in the original mask
           (accounting for any inserted reference voxel from `to-lcm`),
           and reassemble a volume with the mask's geometry, one frame
           per remaining CSV column, non-mask voxels zero. Output format
           is determined by the extension given to --o (surfa dispatches
           automatically), same as any other surfa-based save.

IMPORTANT, UNVERIFIED AGAINST THE REAL LIBRARY: all volume I/O here is
written against surfa's documented API (`surfa.load_volume`, `Volume.data`,
`Volume.save`, `Volume.geom`) as best known, but `surfa` was not installable
in the environment this file was developed in (no network access) and so
these calls have NOT been executed against the real library. Verify the few
`sf.*` call sites (marked "SURFA API:") against your actual surfa version,
under `fspython`, before trusting this in production. Everything else
(MIDAS<->LCM math, mask index bookkeeping, control-file generation,
LCModel invocation, CSV parsing) is pure Python/numpy and has been tested
directly with synthetic data.

Mask voxel ordering: uses Fortran (column-major) flattening throughout, to
match the existing MATLAB pipeline's convention (`find(mask.vol)` in
matlab/lcmodel.m's loadresults() returns column-major linear indices).
This choice is what makes CSV "Col" <-> mask-voxel round-tripping correct;
if you introduce any other mask-flattening code elsewhere, make sure it
also uses Fortran order, or the round trip will silently corrupt voxel
correspondence.
"""

from __future__ import annotations

import argparse
import csv as csv_module
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Optional

import numpy as np

from lcm_control import LCMControl, build_control, parse_control_file, write_control


# -----------------------------------------------------------------------------
# MATLAB-compatible rounding
# -----------------------------------------------------------------------------

def matlab_round(x: float) -> int:
    """MATLAB's round() rounds half-integers away from zero. Python's
    built-in round() uses banker's rounding (round-half-to-even), which
    silently disagrees with MATLAB for exact .5 cases. Replicate MATLAB's
    behavior exactly, since this is used to reproduce midas2lcmodel()'s
    round(nfid/2) center-shift for the water-reference spectrum, and a
    1-sample disagreement there would silently misalign the reference
    spectrum's time-domain data."""
    if x >= 0:
        return math.floor(x + 0.5)
    return math.ceil(x - 0.5)


# -----------------------------------------------------------------------------
# MIDAS -> LCModel FID conversion (mirrors lcmodel.m's midas2lcmodel)
# -----------------------------------------------------------------------------

def midas_to_lcmodel_fid(midas_spectrum: np.ndarray, fidtype: str) -> np.ndarray:
    """
    Convert a MIDAS-format frequency-domain spectrum to an LCModel-format
    time-domain FID, matching matlab/lcmodel.m's midas2lcmodel() exactly.

    midas_spectrum: complex array, shape (nfid, nvox) -- nfid frequency
        points (or time points post-conversion) per voxel, nvox voxels as
        columns.
    fidtype: 'met' or 'ref'. For 'ref' the spectrum is inserted centered
        (shifted by round(nfid/2)) into a double-length buffer before the
        inverse FFT; for 'met' it is inserted at the start (no shift).

    Returns: complex array, shape (nfid, nvox) -- the LCModel time-domain
    FID, truncated back to nfid points after the inverse FFT (matching
    MATLAB's `lcmfid = lcmfid(1:nfid,:)`).
    """
    if fidtype not in ("met", "ref"):
        raise ValueError(f"fidtype must be 'met' or 'ref', got {fidtype!r}")

    nfid, nvox = midas_spectrum.shape
    spec2 = np.zeros((2 * nfid, nvox), dtype=complex)

    if fidtype == "ref":
        shift = matlab_round(nfid / 2)
    else:
        shift = 0
    spec2[shift:shift + nfid, :] = midas_spectrum

    # MATLAB's ifft() operates along columns (dim 1) by default, matching
    # numpy's ifft(..., axis=0) for a (samples, voxels) array.
    lcm_fid = np.fft.ifft(spec2, axis=0)
    lcm_fid = lcm_fid[:nfid, :]
    return lcm_fid


# -----------------------------------------------------------------------------
# Mask <-> flat-index bookkeeping
# -----------------------------------------------------------------------------

def mask_voxel_indices(mask_data: np.ndarray, threshold: float = 0.5) -> np.ndarray:
    """
    Return the 0-based Fortran-order flat indices of voxels where
    mask_data > threshold, in ascending Fortran-order (column-major) index
    order -- i.e. the same order MATLAB's `find(mask.vol)` would produce.

    This ordering is what write_lcm_raw()'s voxel columns are built in, and
    what from-csv's "Col" (1-based) field is interpreted against -- keep
    both ends of the pipeline using this exact function/convention, don't
    reimplement the flattening separately anywhere else.
    """
    flat = mask_data.ravel(order="F")
    return np.flatnonzero(flat > threshold)


# -----------------------------------------------------------------------------
# LCModel .RAW file writing (mirrors lcmodel.m's write_fid, generalized to
# multi-voxel: header written once, then all voxels' data concatenated, per
# the manual's documented multi-voxel .RAW convention, Sec 5.2.3.1)
# -----------------------------------------------------------------------------

def write_lcm_raw(path: str, fid_data: np.ndarray, echot: float, hzpppm: float,
                   ident: str, fmtdat: str, volume: str, tramp: float) -> None:
    """
    Write an LCModel .RAW file: one SEQPAR + NMID namelist header, followed
    by all voxels' complex time-domain data concatenated in column
    (voxel-major) order -- i.e. voxel 1's full nunfil-point data, then
    voxel 2's, etc. -- one real,imag pair per line (fmtdat='(2E15.6)' is the
    manual's own recommended safe choice for multi-voxel files, Sec 5.2.3.1).

    fid_data: complex array, shape (nunfil, nvoxels).
    """
    nunfil, nvox = fid_data.shape
    with open(path, "w") as fp:
        fp.write(" $SEQPAR\n")
        fp.write(f" echot={echot:6.2f}\n")
        fp.write(f" hzpppm={hzpppm:12.4e}\n")
        fp.write(" $END\n")
        fp.write(" $NMID\n")
        fp.write(f" id= {ident}, fmtdat= '({fmtdat})'\n")
        fp.write(f" volume = {volume}\n")
        fp.write(f" tramp= {tramp:g}\n")
        fp.write(" $END\n")
        # Column-major flatten to match MATLAB's `fid(:)`: voxel 1's full
        # column of nunfil points, then voxel 2's, etc.
        flat = fid_data.reshape(-1, order="F")
        for val in flat:
            fp.write(f"{val.real:15.6e}{val.imag:15.6e}\n")


# -----------------------------------------------------------------------------
# SURFA-dependent I/O (unverified against the real library -- see module
# docstring)
# -----------------------------------------------------------------------------

def load_volume_pair(real_path: str, imag_path: str):
    """SURFA API: load a real+imaginary volume pair (any surfa-supported
    format -- nifti, mgz, etc., dispatched automatically by extension), verify
    matching geometry, and return (complex_data_array, geometry_object).

    complex_data_array shape: whatever surfa returns for a single volume's
    .data, typically (X, Y, Z, T) or (X, Y, Z) for a single frame.
    """
    import surfa as sf  # deferred import: only needed for the volume-I/O path

    real_vol = sf.load_volume(real_path)
    imag_vol = sf.load_volume(imag_path)

    if real_vol.data.shape != imag_vol.data.shape:
        raise ValueError(
            f"Real ({real_path}, shape {real_vol.data.shape}) and imaginary "
            f"({imag_path}, shape {imag_vol.data.shape}) volumes must have "
            f"the same shape/geometry.")

    complex_data = real_vol.data.astype(np.float64) + 1j * imag_vol.data.astype(np.float64)
    return complex_data, real_vol.geom


def load_mask(mask_path: str):
    """SURFA API: load the mask volume, return (data_array, geometry_object)."""
    import surfa as sf

    mask_vol = sf.load_volume(mask_path)
    return mask_vol.data, mask_vol.geom


def check_geometry_match(geom_a, geom_b, name_a: str, name_b: str) -> None:
    """SURFA API: raise a clear error if two geometries don't match.

    Uses surfa.transform.image_geometry_equal() (confirmed present in
    surfa 0.6.3, NOT exposed as a top-level sf.* name -- only reachable
    via sf.transform.image_geometry_equal), which compares shape,
    voxsize, center, rotation, shear, and the full vox2world affine
    matrix via np.allclose.

    The original placeholder here compared `geom.shape != geom.shape`
    directly -- confirmed to crash unconditionally (ValueError: "the
    truth value of an array with more than one element is ambiguous"),
    because ImageGeometry.shape returns a numpy array, not a plain
    tuple (unlike ndarray.data.shape, which IS a plain tuple). It
    crashed on a shape MATCH just as reliably as on a mismatch, since
    numpy's array `!=`/`if` ambiguity doesn't depend on the content.

    tol=1e-4, not the function's own default of tol=0.0 (exact
    bit-for-bit match): confirmed against real data (two files from the
    same MRSI pipeline, ku-mrsi/metasurfer/tdad-046-01/vol.lcm/
    mask.nii.gz vs lcm.met.real.nii.gz) that tol=0.0 rejects
    geometrically-identical files over float32 round-off noise alone --
    measured max abs diff there was ~7.6e-6 in vox2world/center and
    ~4.8e-7 in voxsize, consistent with single-precision serialization
    noise between two pipeline steps writing "the same" header, not a
    real misalignment. 1e-4 clears that with over an order of magnitude
    of margin while still catching genuine mismatches (which are
    physically going to be off by far more than float32 noise).
    """
    import surfa as sf

    if not sf.transform.image_geometry_equal(geom_a, geom_b, tol=1e-4):
        def _fmt(g):
            return (f"shape={tuple(int(x) for x in g.shape)} "
                    f"voxsize={tuple(float(x) for x in g.voxsize)} "
                    f"center={tuple(float(x) for x in g.center)}\n"
                    f"    rotation=\n{g.rotation}\n"
                    f"    vox2world.matrix=\n{g.vox2world.matrix}")
        raise ValueError(
            f"Geometry mismatch between {name_a} and {name_b} "
            f"(tolerance 1e-4 exceeded):\n"
            f"  {name_a}: {_fmt(geom_a)}\n"
            f"  {name_b}: {_fmt(geom_b)}")


def save_volume(path: str, data: np.ndarray, geom) -> None:
    """SURFA API: construct and save a volume with the given geometry."""
    import surfa as sf

    vol = sf.Volume(data, geometry=geom)
    vol.save(path)


# -----------------------------------------------------------------------------
# Stage 1: to-lcm
# -----------------------------------------------------------------------------

def cmd_to_lcm(args) -> int:
    real_path, imag_path = args.i
    complex_data, data_geom = load_volume_pair(real_path, imag_path)
    mask_data, mask_geom = load_mask(args.mask)
    check_geometry_match(data_geom, mask_geom, "input data", "mask")

    if complex_data.ndim == 3:
        complex_data = complex_data[..., np.newaxis]
    if complex_data.ndim != 4:
        raise ValueError(
            f"Expected 3D or 4D (X,Y,Z,T) input data, got shape "
            f"{complex_data.shape}")

    nx, ny, nz, nframes = complex_data.shape
    spatial = nx * ny * nz
    # Flatten spatial dims (Fortran order, matching mask_voxel_indices) and
    # keep frames as the FID/time-point axis -- resulting shape
    # (spatial, nframes), then transpose to (nframes, spatial) to match
    # LCModel's (nunfil, nvoxels) convention.
    flat_spatial = complex_data.reshape(spatial, nframes, order="F")
    voxel_idx = mask_voxel_indices(mask_data)
    if voxel_idx.size == 0:
        print("lcm-convert: mask contains no voxels above threshold 0.5",
              file=sys.stderr)
        return 1
    voxel_data = flat_spatial[voxel_idx, :].T  # (nframes, nvoxels_in_mask)

    if not args.no_midas:
        voxel_data = midas_to_lcmodel_fid(voxel_data, args.type)

    n_inserted = 0
    if args.insert_first_voxel:
        ref_real_path, ref_imag_path = args.insert_first_voxel
        ref_complex, ref_geom = load_volume_pair(ref_real_path, ref_imag_path)
        # Expect a single-voxel reference: squeeze all spatial dims to 1.
        ref_flat = ref_complex.reshape(-1, ref_complex.shape[-1]
                                        if ref_complex.ndim == 4 else 1,
                                        order="F")
        if ref_flat.shape[0] != 1:
            # More than one spatial location in the reference file -- error
            # rather than silently averaging or picking one for the user.
            print(f"lcm-convert: --insert-first-voxel file has "
                  f"{ref_flat.shape[0]} spatial voxels, expected exactly 1 "
                  f"(a single reference voxel). Aborting rather than "
                  f"guessing which one you meant.", file=sys.stderr)
            return 1
        ref_voxel = ref_flat.T  # (nframes, 1)
        if not args.no_midas:
            ref_voxel = midas_to_lcmodel_fid(ref_voxel, args.type)
        voxel_data = np.concatenate([ref_voxel, voxel_data], axis=1)
        n_inserted = 1

    nunfil = voxel_data.shape[0]
    write_lcm_raw(
        args.lcm, voxel_data,
        echot=args.echot, hzpppm=args.hzpppm, ident=args.id or "",
        fmtdat=args.fmtdat, volume=args.volume, tramp=args.tramp,
    )

    meta = {
        "n_voxels_in_mask": int(voxel_idx.size),
        "n_inserted_reference_voxels": n_inserted,
        "mask_path": str(Path(args.mask).resolve()),
        "mask_shape": list(mask_data.shape),
        "flatten_order": "F",
        "nunfil": int(nunfil),
    }
    meta_path = str(Path(args.lcm).with_suffix(Path(args.lcm).suffix + ".meta.json"))
    Path(meta_path).write_text(json.dumps(meta, indent=2))
    print(f"lcm-convert: wrote {args.lcm} ({voxel_data.shape[1]} voxels, "
          f"{nunfil} points/voxel) and {meta_path}", file=sys.stderr)
    return 0


# -----------------------------------------------------------------------------
# Stage 2 (of 3): run -- invoke lcmodel now that met.lcm AND h2o.lcm both
# exist (to-lcm only ever produces one file type at a time, per its own
# --type met|ref flag, so running lcmodel is necessarily a separate step
# that needs both outputs already on disk)
# -----------------------------------------------------------------------------

def run_lcmodel(lcmodel_bin: str, met_path: str, h2o_path: Optional[str],
                 control_path: str, csv_path: str) -> int:
    """Invoke lcmodel with -i <met> <h2o> -csv <csv>, control file on stdin.

    Requires the -i/-csv CLI flags built in this project's goal #3 -- an
    lcmodel binary without those flags will reject this invocation with an
    'unrecognized command-line argument' error.

    IMPORTANT: lcmodel's own process exit code is NOT a reliable success
    signal. Confirmed directly (end-to-end test, a real NAMELIST parse
    error): LCModel's fatal-error path for ILEVEL<0 (source/LCModel.f's
    ERRMES, ~line 2950) executes a bare Fortran `STOP` with no explicit
    code, which returns exit status 0 to the OS regardless -- the same
    "reports success while having actually failed" class of issue this
    project's Fortran work already found for ILEVEL>0 fatal errors
    (EXITPS's STOP being commented out there instead). So a 0 returncode
    here does NOT mean the run actually produced results -- it only means
    the process didn't segfault. The actual signal checked is whether
    csv_path was created and is non-empty.
    """
    cmd = [lcmodel_bin, "-i", met_path, h2o_path or "", "-csv", csv_path]
    with open(control_path, "rb") as control_fp:
        result = subprocess.run(cmd, stdin=control_fp, capture_output=True)
    sys.stderr.buffer.write(result.stdout)
    sys.stderr.buffer.write(result.stderr)
    csv_ok = Path(csv_path).is_file() and Path(csv_path).stat().st_size > 0
    if result.returncode != 0:
        return result.returncode
    if not csv_ok:
        print(f"lcm-convert: lcmodel exited 0 but did not produce a "
              f"non-empty {csv_path} -- treating as failure (lcmodel's "
              f"exit code is not reliable for its own fatal-error paths; "
              f"see stderr above for the actual error).", file=sys.stderr)
        return 1
    return 0


def cmd_run(args) -> int:
    return run_lcmodel(args.lcmodel_bin, args.met, args.h2o, args.control,
                        args.csv)


# -----------------------------------------------------------------------------
# Stage 2: from-csv
# -----------------------------------------------------------------------------

def read_lcm_csv(path: str):
    """Read an LCModel .csv output file. Returns (col_indices, field_names,
    values) where col_indices is a list of 1-based mask-space column indices
    (the CSV's own "Col" field) and values[field] is a list of floats, one
    per row, in the same row order as col_indices."""
    with open(path, newline="") as fp:
        reader = csv_module.DictReader(fp)
        fieldnames = reader.fieldnames
        if fieldnames is None:
            raise ValueError(f"{path}: no header row found")
        # LCModel's own CSV header is written as 'Row, Col, <metab1>, ...'
        # (comma-SPACE separated, source/LCModel.f's csv_line construction),
        # so every field name except the first has a leading space -- e.g.
        # the raw header key is literally ' Col', not 'Col'. Re-assigning
        # reader.fieldnames (confirmed to work: DictReader uses whatever
        # fieldnames are current when it builds each row's dict) makes
        # every row dict keyed by the stripped names consistently, instead
        # of stripping only for this presence check and then crashing with
        # a KeyError on every real LCModel CSV when accessing row["Col"].
        fieldnames = [f.strip() for f in fieldnames]
        reader.fieldnames = fieldnames
        if "Col" not in fieldnames:
            raise ValueError(
                f"{path}: expected a 'Col' column, found {fieldnames}")
        data_fields = [f for f in fieldnames if f not in ("Row", "Col")]
        col_indices = []
        values = {f: [] for f in data_fields}
        for row in reader:
            col_indices.append(int(float(row["Col"])))
            for f in data_fields:
                raw = row[f].strip()
                values[f].append(float(raw) if raw not in ("", None) else float("nan"))
    return col_indices, data_fields, values


def check_single_column_grid(control_path: str) -> None:
    """Hard runtime guard: from-csv's Col-is-a-flat-voxel-index mapping
    (see mask_voxel_indices()/cmd_from_csv() below) is only correct
    because idcol reduces to the flat sequential voxel index when
    ndrows=ndslic=1 (confirmed directly against source/LCModel.f's
    ring-scan loop: ivoxel = (idslic-1)*ndrows*ndcols +
    (idrow-1)*ndcols + idcol, collapsing to exactly idcol only in this
    degenerate-grid case). For any other grid shape, Col is just the
    column-within-row coordinate, NOT a flat index -- reusing this same
    mapping would silently misassign every voxel to the wrong mask
    position, with no error, no crash, just wrong output. Error out
    loudly instead of ever letting that happen silently."""
    raw = parse_control_file(control_path)
    ndrows = int(raw.get("ndrows", "1"))
    ndslic = int(raw.get("ndslic", "1"))
    if ndrows != 1 or ndslic != 1:
        print(f"lcm-convert: {control_path} has ndrows={ndrows}, "
              f"ndslic={ndslic} -- from-csv's Col-to-mask-voxel mapping "
              f"is only valid for ndrows=1 AND ndslic=1 (a single column "
              f"of voxels). Reusing it for any other grid shape would "
              f"silently misassign voxel data to the wrong mask "
              f"positions. Refusing to proceed.", file=sys.stderr)
        sys.exit(1)


def cmd_from_csv(args) -> int:
    check_single_column_grid(args.control)
    col_indices, data_fields, values = read_lcm_csv(args.csv)
    mask_data, mask_geom = load_mask(args.mask)
    voxel_idx = mask_voxel_indices(mask_data)  # 0-based flat Fortran indices

    n_offset = 0
    if args.meta:
        meta = json.loads(Path(args.meta).read_text())
        n_offset = meta.get("n_inserted_reference_voxels", 0)
        if meta.get("n_voxels_in_mask") != voxel_idx.size:
            print(f"lcm-convert: warning -- meta file records "
                  f"{meta.get('n_voxels_in_mask')} mask voxels, but current "
                  f"--mask has {voxel_idx.size}. Continuing, but verify this "
                  f"is the same mask used in the to-lcm step.",
                  file=sys.stderr)

    spatial_shape = mask_data.shape
    spatial = int(np.prod(spatial_shape))
    nframes = len(data_fields)
    out_flat = np.zeros((spatial, nframes), dtype=np.float64)

    n_written = 0
    for row_i, col_1based in enumerate(col_indices):
        # col_1based is 1-based into the augmented (possibly
        # reference-voxel-prepended) voxel list this session's LCModel run
        # actually saw. Subtract the inserted-reference offset, then it's
        # 1-based into voxel_idx (the mask's own voxel list).
        mask_pos = col_1based - 1 - n_offset  # 0-based index into voxel_idx
        if mask_pos < 0:
            # This row corresponds to the inserted reference voxel itself,
            # not a real mask voxel -- there's nowhere in the output volume
            # for it to go. Skip it.
            continue
        if mask_pos >= voxel_idx.size:
            print(f"lcm-convert: warning -- CSV row {row_i} has Col="
                  f"{col_1based}, out of range for {voxel_idx.size} mask "
                  f"voxels (with reference-offset {n_offset}); skipping.",
                  file=sys.stderr)
            continue
        flat_idx = voxel_idx[mask_pos]
        for f_i, f in enumerate(data_fields):
            out_flat[flat_idx, f_i] = values[f][row_i]
        n_written += 1

    out_data = out_flat.reshape(*spatial_shape, nframes, order="F")
    save_volume(args.o, out_data, mask_geom)
    print(f"lcm-convert: wrote {args.o} ({n_written}/{voxel_idx.size} mask "
          f"voxels populated, {nframes} frames: {', '.join(data_fields)})",
          file=sys.stderr)
    return 0


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="lcm-convert")
    sub = p.add_subparsers(dest="command", required=True)

    p_to = sub.add_parser("to-lcm", help="Convert one (met or ref) image volume pair "
                                          "real+imaginary pair, masked, to "
                                          "an LCModel .RAW file (see `run` "
                                          "to actually invoke lcmodel once "
                                          "both met and ref .lcm files "
                                          "exist)")
    p_to.add_argument("--i", nargs=2, required=True,
                       metavar=("REAL_VOL", "IMAG_VOL"),
                       help="Real and imaginary input image volumes (any surfa-supported format: nifti, mgz, etc.)")
    p_to.add_argument("--mask", required=True, help="Binary mask image volume "
                       "(voxels with value > 0.5 are converted)")
    p_to.add_argument("--type", choices=["met", "ref"], required=True,
                       help="'met' (metabolite/water-suppressed) or 'ref' "
                            "(water reference) -- matches midas2lcmodel's "
                            "fidtype convention")
    p_to.add_argument("--lcm", required=True, help="Output LCModel .RAW "
                       "file path; a sidecar <path>.meta.json is also "
                       "written, needed by `from-csv` if --insert-first-voxel "
                       "was used")
    p_to.add_argument("--no-midas", action="store_true",
                       help="Input is already LCModel-format time-domain "
                            "FID data; skip MIDAS-to-LCModel conversion")
    p_to.add_argument("--insert-first-voxel", nargs=2, default=None,
                       metavar=("REAL_VOL", "IMAG_VOL"),
                       help="Single-voxel real/imaginary image volume pair to "
                            "prepend as voxel 1, ahead of the mask-derived "
                            "voxels (e.g. to seed LCModel's default "
                            "Bayesian-learning prior-propagation from a "
                            "known-good reference)")
    # Control-parameter passthroughs needed to build the .RAW header and the
    # control file for this run.
    p_to.add_argument("--echot", type=float, default=0.0)
    p_to.add_argument("--hzpppm", type=float, default=0.0)
    p_to.add_argument("--id", type=str, default="")
    p_to.add_argument("--fmtdat", type=str, default="2E15.6")
    p_to.add_argument("--volume", type=str, default="1.0")
    p_to.add_argument("--tramp", type=float, default=1.0)
    p_to.set_defaults(func=cmd_to_lcm)

    p_run = sub.add_parser("run", help="Invoke lcmodel on an existing "
                                        "met.lcm + h2o.lcm pair, using this "
                                        "session's -i/-csv CLI overrides "
                                        "(requires both files already "
                                        "produced by separate to-lcm calls)")
    p_run.add_argument("--met", required=True, help="met .lcm file (from "
                        "`to-lcm --type met`)")
    p_run.add_argument("--h2o", required=True, help="water-reference .lcm "
                        "file (from `to-lcm --type ref`)")
    p_run.add_argument("--control", required=True,
                        help="Control file (build with lcm-control; do not "
                             "put filraw=/filh2o=/filcsv=/ndcols=/icolen= "
                             "in it -- -i/-csv here override those per-run)")
    p_run.add_argument("--csv", required=True, help="Output .csv path")
    p_run.add_argument("--lcmodel-bin", default="lcmodel")
    p_run.set_defaults(func=cmd_run)

    p_from = sub.add_parser("from-csv", help="Reassemble an LCModel .csv "
                                              "fit into an image volume")
    p_from.add_argument("--csv", required=True)
    p_from.add_argument("--control", required=True,
                         help="The control file this CSV was actually "
                              "produced with (required, not optional -- "
                              "used to verify ndrows=1 and ndslic=1, "
                              "since Col only means 'flat mask-voxel "
                              "index' for that specific grid shape; any "
                              "other shape would silently misassign "
                              "voxels if allowed through)")
    p_from.add_argument("--mask", required=True,
                         help="Same mask used in the corresponding to-lcm "
                              "step")
    p_from.add_argument("--meta", default=None,
                         help="Sidecar .meta.json from the to-lcm step "
                              "(needed if --insert-first-voxel was used "
                              "there; safe to omit otherwise)")
    p_from.add_argument("--o", required=True, help="Output image volume path (format determined by extension, e.g. .nii.gz or .mgz)")
    p_from.set_defaults(func=cmd_from_csv)

    return p


def main(argv=None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
