#!/usr/bin/env fspython
"""
lcm_control.py

Run with `fspython` (must be on your PATH) for consistency with
lcm_convert.py. This file has no *unconditional* imaging/surfa dependency
-- it runs under any Python 3 interpreter for ordinary use -- but `surfa`
is imported (lazily, only if `--mask` is given) to derive `ndcols`/`icolen`
from a mask's voxel count; see count_mask_voxels().

LCModel control-file (.control) parameter model, reader, and writer.

Mirrors the parameter set and write_control() logic in matlab/lcmodel.m,
generalized into a reusable Python library plus a `lcm-control` CLI tool.

Precedence when building a control file (highest wins):
    1. Parameter set explicitly on the command line
    2. For ndcols/icolen only: voxel count derived from --mask, if given
    3. Parameter set in an input control file (--i)
    4. Built-in default

Control-file format (Fortran NAMELIST):
    $LCMODL
     key1=value1
     key2='string value'
     ...
    $END

Parameter definitions: see docs/manual.pdf (Provencher), Sec 5.3 in particular.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np


# -----------------------------------------------------------------------------
# Parameter table: (attribute name, NAMELIST key, CLI flag, type, default,
#                    help text)
#
# This single table drives the dataclass fields, the NAMELIST key <-> Python
# attribute mapping, and the argparse CLI -- adding a new control parameter
# only requires adding one row here.
# -----------------------------------------------------------------------------
# type is one of: str, int, float  (control parameters are all scalars here;
# LCModel's array parameters like KEY(20) are not modeled individually --
# use --extra for anything not in this table)

PARAM_TABLE = [
    # (attr,          namelist_key, cli_flag,        type,  default,                 help)
    ("key",           "key",        "--key",         int,   210387309,               "License key"),
    ("title",         "title",      "--title",       str,   "title",                 "Title string for plots/output"),
    ("srcraw",        "srcraw",     "--srcraw",      str,   "",                      "Source-of-raw-file annotation"),
    ("savdir",        "savdir",     "--savdir",      str,   "",                      "Save-directory annotation"),
    ("ppmst",         "ppmst",      "--ppmst",       float, 4.0,                     "Analysis window left (start) edge, ppm"),
    ("ppmend",        "ppmend",     "--ppmend",      float, 1.9,                     "Analysis window right (end) edge, ppm"),
    ("nunfil",        "nunfil",     "--nunfil",      int,   0,                       "Number of complex time-domain points (0 = infer from data)"),
    ("ltable",        "ltable",     "--ltable",      int,   None,                    ".TABLE output unit: 0=suppress (default: 7 if --filtab given, else 0)"),
    ("lps",           "lps",        "--lps",         int,   None,                    ".PS output unit: 0=suppress (default: 8 if --filps given, else 0)"),
    ("lcsv",          "lcsv",       "--lcsv",        int,   11,                      ".CSV output: 0=suppress, 11=make"),
    ("lprint",        "lprint",     "--lprint",      int,   None,                    ".PRINT (detailed) output unit: 0=suppress (default: 6 if --filpri given, else 0)"),
    ("lcoord",        "lcoord",     "--lcoord",      int,   None,                    ".COORD output unit: 0=suppress (default: 9 if --filcoo given, else 0)"),
    ("hzpppm",        "hzpppm",     "--hzpppm",      float, None,                    "REQUIRED (no meaningful default exists): field strength, MHz proton resonance (42.58*B0 Tesla)"),
    ("filtab",        "filtab",     "--filtab",      str,   "",                      "Pathname of .TABLE output file"),
    ("filraw",        "filraw",     "--filraw",      str,   "",                      "Pathname of .RAW (met) input file"),
    ("filh2o",        "filh2o",     "--filh2o",      str,   "",                      "Pathname of water-reference .RAW input file"),
    ("filps",         "filps",      "--filps",       str,   "",                      "Pathname of .PS output file"),
    ("filcsv",        "filcsv",     "--filcsv",      str,   "",                      "Pathname of .CSV output file"),
    ("filbas",        "filbas",     "--filbas",      str,   "",                      "Pathname of .BASIS file"),
    ("filpri",        "filpri",     "--filpri",      str,   "",                      "Pathname of detailed .PRINT output file"),
    ("filcoo",        "filcoo",     "--filcoo",      str,   "",                      "Pathname of .COORD output file"),
    ("neach",         "neach",      "--neach",       int,   99,                      "Number of metabolites for individual plots"),
    ("echot",         "echot",      "--echot",       float, None,                    "REQUIRED (no meaningful default exists): echo time, ms"),
    ("deltat",        "deltat",     "--deltat",      float, None,                    "Sample (dwell) time, s (1/Fs) -- overrides --fs if both given"),
    ("fs",            None,         "--fs",          float, 0.0,                     "Sample frequency, Hz (writes deltat=1/fs; ignored if --deltat given)"),
    ("dows",          "dows",       "--dows",        str,   "T",                     "Do water scaling: T/F"),
    ("doecc",         "doecc",      "--doecc",       str,   "F",                     "Do eddy-current correction: T/F"),
    ("sddegp",        "sddegp",     "--sddegp",      float, None,                    "Expected SD of DEGPPM (default: LCModel's own built-in 20.0 if omitted)"),
    ("nomit",         "nomit",      "--nomit",       int,   None,                    "Number of metabolites to exclude (default: LCModel's own built-in 0 if omitted)"),
    ("dorefs1",       "DOREFS(1)",  "--dorefs1",     str,   "F",                     "DOREFS(1): T=use water peak as landmark"),
    ("atth2o",        "atth2o",     "--atth2o",      float, 1.0,                     "Water-scaling attenuation factor"),
    ("wconc",         "wconc",      "--wconc",       float, 39590.0,                 "NMR-visible water concentration, mM (tissue-dependent)"),
    ("pgnorm",        "pgnorm",     "--pgnorm",      str,   "US",                    "Page style: US or A4/EU"),
    ("ipage2",        "ipage2",     "--ipage2",      int,   1,                       "Second-page behavior: 0=suppress,1=if needed,2=always"),
    ("ndcols",        "ndcols",     "--ndcols",      int,   1,                       "Number of columns in the CSI data set"),
    ("ndrows",        "ndrows",     "--ndrows",      int,   1,                       "Number of rows in the CSI data set"),
    ("ndslic",        "ndslic",     "--ndslic",      int,   1,                       "Number of slices in the CSI data set"),
    ("icolst",        "icolst",     "--icolst",      int,   1,                       "First column to analyze"),
    ("icolen",        "icolen",     "--icolen",      int,   None,                    "Last column to analyze (default: ndcols)"),
    ("irowst",        "irowst",     "--irowst",      int,   1,                       "First row to analyze"),
    ("irowen",        "irowen",     "--irowen",      int,   None,                    "Last row to analyze (default: ndrows)"),
    ("islice",        "islice",     "--islice",      int,   1,                       "Slice number to analyze"),
    ("nvoxsk",        "nvoxsk",     "--nvoxsk",      int,   0,                       "Number of voxels to skip"),
]
# NOTE: volume/tramp/id/fmtdat were removed from this table -- confirmed
# absent from NAMELIST /LCMODL/ (source/nml_lcmodl.inc); they are fields
# of the .RAW file's own $NMID namelist, already correctly handled by
# lcm_convert.py's `to-lcm` subcommand (--volume/--tramp/--id/--fmtdat).
# Including them here caused write_control() to emit them into the
# $LCMODL control file, which crashed LCModel's NAMELIST read
# unconditionally (an unrecognized-variable-name parse error, confirmed
# via direct end-to-end testing) -- not a cosmetic issue, a hard failure
# on every control file this function ever produced.

# Toggles controlling which optional output blocks get written (matching
# lcmodel.m's enableLCM*Output flags). These do not correspond to a single
# NAMELIST key each -- they gate whether filpri/lprint, filcoo/lcoord, and
# neach get written at all.
_TOGGLE_DEFAULTS = {
    "enable_detailed_output": True,   # filpri + lprint
    "enable_spec_output": True,      # filcoo + lcoord
    "enable_all_spec_output": True,  # neach
}


def _attr_names():
    return [row[0] for row in PARAM_TABLE]


def _make_dataclass():
    """Build the LCMControl dataclass fields dynamically from PARAM_TABLE."""
    annotations = {}
    defaults = {}
    for attr, _key, _flag, typ, default, _help in PARAM_TABLE:
        annotations[attr] = Optional[typ] if default is None else typ
        defaults[attr] = default
    for name, default in _TOGGLE_DEFAULTS.items():
        annotations[name] = bool
        defaults[name] = default
    annotations["extra"] = dict
    defaults["extra"] = None  # filled in __post_init__

    ns = {"__annotations__": annotations}
    ns.update(defaults)

    def __post_init__(self):
        if self.extra is None:
            self.extra = {}

    ns["__post_init__"] = __post_init__
    cls = dataclass(type("LCMControl", (), ns))
    return cls


LCMControl = _make_dataclass()
LCMControl.__doc__ = """
LCModel control-parameter set, mirroring matlab/lcmodel.m's property list.

Unrecognized control-file keys (anything not in PARAM_TABLE) are preserved
verbatim in `.extra` (dict of namelist_key -> raw string value) so the tool
stays usable even for control parameters this table doesn't know about yet.
"""


_KEY_TO_ATTR = {row[1]: row[0] for row in PARAM_TABLE if row[1] is not None}
_ATTR_TO_KEY = {row[0]: row[1] for row in PARAM_TABLE if row[1] is not None}


# -----------------------------------------------------------------------------
# Control-file parsing
# -----------------------------------------------------------------------------

import re

_KEY_PATTERN = re.compile(r"[A-Za-z_][A-Za-z0-9_]*(?:\(\d+\))?\s*=")
_QUOTED_SPAN_PATTERN = re.compile(r"'[^']*'")


def _quoted_spans(line: str):
    return [(m.start(), m.end()) for m in _QUOTED_SPAN_PATTERN.finditer(line)]


def parse_control_file(path: str) -> dict:
    """
    Parse an existing .control (NAMELIST LCMODL) file into a flat dict of
    {namelist_key_lower: raw_value_string}.

    Tolerant, line-oriented parser -- not a full Fortran NAMELIST parser.
    Handles:
      - '$LCMODL' / '$END' delimiter lines (ignored)
      - 'KEY=VALUE' or 'KEY = VALUE' (case-insensitive key)
      - single-quoted string values, e.g. FILRAW='foo.raw', including ones
        that happen to contain characters that would otherwise look like a
        new 'key=' assignment (quoted spans are excluded from key scanning)
      - fixed-width/padded values with a space right after '=', e.g.
        'echot= 20.10' (common from Fortran-style %6.2f-ish formatting --
        a naive whitespace-split tokenizer breaks on this; this parser
        locates each 'KEY=' occurrence directly via regex instead, and
        takes everything up to the *next* 'KEY=' occurrence as the value,
        so embedded spaces in an unquoted numeric value are preserved)
      - multiple key=value pairs on one line (space- or comma-separated,
        matching LCModel's own convention, e.g. 'KEY(1)=1 KEY(2)=2')
    Does not evaluate Fortran array-repeat syntax (e.g. 3*1.0).
    """
    raw: dict[str, str] = {}
    text = Path(path).read_text()

    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("$"):
            continue
        if line.lower().startswith("c ") or line.startswith("!"):
            continue

        quoted_spans = _quoted_spans(line)

        def _inside_quote(pos: int) -> bool:
            return any(start <= pos < end for start, end in quoted_spans)

        matches = [m for m in _KEY_PATTERN.finditer(line)
                   if not _inside_quote(m.start())]
        for i, m in enumerate(matches):
            key = m.group(0)[:-1].strip()  # strip trailing '='
            # value spans from the end of this match to the start of the
            # next KEY= match (or end of line)
            val_start = m.end()
            val_end = matches[i + 1].start() if i + 1 < len(matches) else len(line)
            val = line[val_start:val_end].strip()
            val = val.rstrip(",").strip()
            if len(val) >= 2 and val[0] == "'" and val[-1] == "'":
                val = val[1:-1]
            if key:
                raw[key.lower()] = val
    return raw


def _coerce(value: str, typ):
    if typ is int:
        return int(float(value))  # tolerate '11.0' etc.
    if typ is float:
        # Fortran real literals sometimes use D instead of E for exponent
        return float(value.replace("D", "E").replace("d", "e"))
    return value


def apply_raw_dict(ctrl: "LCMControl", raw: dict) -> None:
    """Apply a {namelist_key_lower: raw_str} dict onto an LCMControl instance,
    in place. Recognized keys are coerced to their proper type; unrecognized
    keys are stored verbatim in ctrl.extra."""
    remaining = dict(raw)
    for attr, key, _flag, typ, _default, _help in PARAM_TABLE:
        if key is None:
            continue
        lk = key.lower()
        if lk in remaining:
            val = remaining.pop(lk)
            setattr(ctrl, attr, _coerce(val, typ))
    # anything left over is unrecognized -- preserve verbatim
    ctrl.extra.update(remaining)


# -----------------------------------------------------------------------------
# Precedence resolution
# -----------------------------------------------------------------------------

def build_control(cli_args: dict, input_control_path: Optional[str] = None) -> "LCMControl":
    """
    Build an LCMControl applying, in order (each overwriting the previous):
      1. built-in defaults
      2. --i input control file (if given)
      3. explicit CLI overrides (cli_args: {attr_name: value}, only entries
         that were *actually specified* on the command line -- callers must
         not include argparse defaults that the user didn't type)
    """
    ctrl = LCMControl()
    if input_control_path:
        raw = parse_control_file(input_control_path)
        apply_raw_dict(ctrl, raw)
    for attr, value in cli_args.items():
        if value is not None:
            setattr(ctrl, attr, value)
    return ctrl


# -----------------------------------------------------------------------------
# Control-file writing
# -----------------------------------------------------------------------------

def write_control(ctrl: "LCMControl", out_path: str) -> None:
    """Write ctrl out as a $LCMODL NAMELIST control file, matching the
    field order and conventions of matlab/lcmodel.m's write_control()."""
    # hzpppm/echot have no meaningful default -- unlike nunfil (which
    # to-lcm can derive from the real data it loaded, via --meta) there
    # is no data source either of these could ever be auto-derived from;
    # they're pure acquisition metadata. LCModel's own BLOCK DATA
    # defaults (HZPPPM=84.47, ECHOT=-1.) are not sensible fallbacks --
    # 84.47 MHz is just whatever value happened to be hardcoded, and -1.
    # ms is a plainly invalid sentinel -- so silently omitting either
    # would produce a control file that looks valid but analyzes with
    # the wrong physics. Required, not defaulted: raised here (not just
    # checked in main()) so this guarantee holds for every caller of
    # this function, not only the CLI entry point.
    missing = [name for name, val in
               (("hzpppm", ctrl.hzpppm), ("echot", ctrl.echot))
               if val is None]
    if missing:
        raise ValueError(
            f"write_control: {', '.join(missing)} must be set -- no "
            f"meaningful default exists (LCModel's own BLOCK DATA "
            f"defaults for these, 84.47 and -1. respectively, are not "
            f"safe fallbacks). Pass --{missing[0]} explicitly (or set it "
            f"via an --i input control file).")
    icolen = ctrl.icolen if ctrl.icolen is not None else ctrl.ndcols
    irowen = ctrl.irowen if ctrl.irowen is not None else ctrl.ndrows
    if ctrl.deltat is not None:
        deltat = ctrl.deltat
    elif ctrl.fs:
        deltat = 1.0 / ctrl.fs
    else:
        deltat = None

    # ltable/lps/lcoord/lprint: each is a Fortran UNIT NUMBER (0=disabled,
    # nonzero=enabled AND the actual unit used for that OPEN) paired with a
    # FIL* path variable that defaults to "" (unset). The old code wrote
    # these four unconditionally at their "enabled" PARAM_TABLE defaults
    # (7/8/9/6) regardless of whether the matching path was ever given --
    # a self-contradictory default combination that fatally crashes
    # LCModel's own validation (e.g. source/LCModel.f:1644-1645,
    # `IF (lps.gt.0 .and. filps.eq.' ') call errmes(3,-4,chsubp)`, a real
    # STOP -- confirmed via direct end-to-end testing) unless the user
    # remembers to also supply the path every time.
    #
    # Fixed default: derive "enabled" from whether the path was actually
    # given (None means the user didn't touch this flag at all), so the
    # untouched-defaults case is self-consistent (both off) instead of
    # self-contradictory. An explicit --ltable/--lps/--lcoord/--lprint
    # value (including 0) always wins outright, same as before.
    ltable = ctrl.ltable if ctrl.ltable is not None else (7 if ctrl.filtab else 0)
    lps = ctrl.lps if ctrl.lps is not None else (8 if ctrl.filps else 0)
    lprint = ctrl.lprint if ctrl.lprint is not None else (6 if ctrl.filpri else 0)
    lcoord = ctrl.lcoord if ctrl.lcoord is not None else (9 if ctrl.filcoo else 0)

    lines = [" $LCMODL"]
    lines.append(f" key={ctrl.key}")
    lines.append(f" title='{ctrl.title}'")
    if ctrl.srcraw:
        lines.append(f" srcraw='{ctrl.srcraw}'")
    if ctrl.savdir:
        lines.append(f" savdir='{ctrl.savdir}'")
    lines.append(f" ppmst={ctrl.ppmst:f}")
    lines.append(f" ppmend={ctrl.ppmend:f}")
    if ctrl.nunfil:
        lines.append(f" nunfil={ctrl.nunfil}")
    lines.append(f" ltable={ltable}")
    lines.append(f" lps={lps}")
    lines.append(f" lcsv={ctrl.lcsv}")
    lines.append(f" hzpppm={ctrl.hzpppm:.4e}")
    if ctrl.filtab:
        lines.append(f" filtab='{ctrl.filtab}'")
    if ctrl.filraw:
        lines.append(f" filraw='{ctrl.filraw}'")
    if ctrl.filh2o:
        lines.append(f" filh2o='{ctrl.filh2o}'")
    if ctrl.filps:
        lines.append(f" filps='{ctrl.filps}'")
    if ctrl.filcsv:
        lines.append(f" filcsv='{ctrl.filcsv}'")
    if ctrl.filbas:
        lines.append(f" filbas='{ctrl.filbas}'")
    # enable_detailed_output/enable_spec_output remain available as an
    # explicit force-OFF override (e.g. a path was given but the user
    # still wants that block suppressed) -- but, unlike before, they are
    # no longer the ONLY thing standing between the default state and a
    # guaranteed crash; the ltable/lps/lcoord/lprint derivation above
    # already makes "nothing given" self-consistently safe.
    if ctrl.enable_detailed_output:
        if ctrl.filpri:
            lines.append(f" filpri='{ctrl.filpri}'")
        lines.append(f" lprint={lprint}")
    if ctrl.enable_spec_output:
        if ctrl.filcoo:
            lines.append(f" filcoo='{ctrl.filcoo}'")
        lines.append(f" lcoord={lcoord}")
    if ctrl.enable_all_spec_output:
        lines.append(f" neach={ctrl.neach}")
    lines.append(f" echot={ctrl.echot:6.2f}")
    if deltat is not None:
        lines.append(f" deltat={deltat:11.3e}")
    lines.append(f" dows={ctrl.dows}")
    lines.append(f" doecc={ctrl.doecc}")
    if ctrl.sddegp is not None:
        lines.append(f" sddegp={ctrl.sddegp:g}")
    if ctrl.nomit is not None:
        lines.append(f" nomit={ctrl.nomit}")
    lines.append(f" DOREFS(1)={ctrl.dorefs1}")
    lines.append(f" atth2o={ctrl.atth2o:g}")
    lines.append(f" wconc={ctrl.wconc:g}")
    lines.append(f" pgnorm='{ctrl.pgnorm}'")
    lines.append(f" ipage2={ctrl.ipage2}")
    lines.append(f" ndcols={ctrl.ndcols}")
    lines.append(f" icolst={ctrl.icolst}")
    lines.append(f" icolen={icolen}")
    lines.append(f" ndrows={ctrl.ndrows}")
    lines.append(f" irowst={ctrl.irowst}")
    lines.append(f" irowen={irowen}")
    lines.append(f" ndslic={ctrl.ndslic}")
    lines.append(f" islice={ctrl.islice}")
    lines.append(f" nvoxsk={ctrl.nvoxsk}")
    for key, val in ctrl.extra.items():
        lines.append(f" {key}={val}")
    lines.append(" $END")

    Path(out_path).write_text("\n".join(lines) + "\n")


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="lcm-control",
        description="Build an LCModel .control file from defaults, an "
                     "optional input control file, and CLI overrides. "
                     "Precedence: explicit --ndcols/--icolen (if given) > "
                     "--mask-derived voxel count (if --mask given) > "
                     "--i input file > built-in defaults. Same pattern for "
                     "nunfil: explicit --nunfil (if given) > "
                     "--meta-derived point count (if --meta given) > "
                     "--i input file > built-in defaults (which is 2048, "
                     "NOT '0 = infer' -- LCModel has no such inference; "
                     "omitting nunfil entirely silently uses its own "
                     "BLOCK DATA default of 2048, which will misalign or "
                     "crash on data with any other point count -- always "
                     "pass --nunfil explicitly or, better, --meta pointing "
                     "at the to-lcm .meta.json for the same data). All "
                     "other parameters: CLI > --i input file > built-in "
                     "defaults.",
    )
    p.add_argument("--o", dest="out_path", required=True,
                   help="Output .control file path")
    p.add_argument("--i", dest="in_path", default=None,
                   help="Input .control file to seed values from")
    p.add_argument("--meta", dest="meta_path", default=None,
                   help="Sidecar .meta.json written by lcm-convert's "
                        "to-lcm step (met or ref -- both record the same "
                        "nunfil for a given run) -- sets nunfil to the "
                        "actual point count that step derived from the "
                        "real input data, closing the '0 = infer' "
                        "landmine below. Overridden by an explicit "
                        "--nunfil on this same command line.")
    p.add_argument("--mask", dest="mask_path", default=None,
                   help="Image volume (any surfa-supported format: nifti, "
                        "mgz, etc.) whose voxel count (value > 0.5) sets "
                        "both ndcols and icolen -- overrides whatever "
                        "ndcols/icolen came from --i's control file, but "
                        "is itself overridden by an explicit --ndcols/"
                        "--icolen on this same command line. Requires "
                        "surfa (only imported if --mask is actually used).")
    for attr, _key, flag, typ, default, help_text in PARAM_TABLE:
        p.add_argument(flag, dest=attr, type=typ, default=None,
                        help=f"{help_text} (default: {default})")
    p.add_argument("--no-detailed-output", dest="enable_detailed_output",
                   action="store_false", default=None,
                   help="Suppress filpri/lprint block")
    p.add_argument("--no-spec-output", dest="enable_spec_output",
                   action="store_false", default=None,
                   help="Suppress filcoo/lcoord block")
    p.add_argument("--no-all-spec-output", dest="enable_all_spec_output",
                   action="store_false", default=None,
                   help="Suppress neach")
    p.add_argument("--extra", action="append", default=[],
                   metavar="KEY=VALUE",
                   help="Additional raw control-file key=value pair "
                        "(repeatable); written through verbatim")
    return p


def count_mask_voxels(mask_path: str, threshold: float = 0.5) -> int:
    """Load a mask (any surfa-supported format) and count voxels with
    value > threshold. Threshold must match lcm_convert.py's
    mask_voxel_indices() default (0.5) -- if you ever change one, change
    the other, or --mask-derived ndcols/icolen here will silently disagree
    with the actual number of voxels lcm-convert extracts for the same
    mask, and the LCModel run downstream would see a voxel-count mismatch
    between the control file and the .RAW file.

    Note: this only needs a COUNT, not a per-voxel index mapping, so
    (unlike lcm_convert.py's mask_voxel_indices) it does not need to care
    about Fortran- vs C-order flattening -- a simple count is order-
    independent."""
    import surfa as sf  # deferred: only needed if --mask is actually used

    mask_vol = sf.load_volume(mask_path)
    count = int(np.count_nonzero(mask_vol.data > threshold))
    return count


def main(argv=None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    # Built directly here (not via build_control()'s generic cli_overrides
    # dict) because the ndcols/icolen interaction needs to actively RESET
    # icolen back to "auto-default to ndcols" (None) in some cases, which a
    # plain {attr: value} overrides dict can't express -- a None entry in
    # that dict just means "don't touch", not "clear this back to default".
    ctrl = LCMControl()
    if args.in_path:
        raw = parse_control_file(args.in_path)
        apply_raw_dict(ctrl, raw)

    if args.meta_path:
        meta = json.loads(Path(args.meta_path).read_text())
        if "nunfil" not in meta:
            print(f"lcm-control: --meta {args.meta_path!r} has no 'nunfil' "
                  f"key -- not a to-lcm .meta.json file?", file=sys.stderr)
            return 1
        ctrl.nunfil = int(meta["nunfil"])

    if args.mask_path:
        count = count_mask_voxels(args.mask_path)
        if count == 0:
            print(f"lcm-control: --mask {args.mask_path!r} contains no "
                  f"voxels above threshold 0.5", file=sys.stderr)
            return 1
        ctrl.ndcols = count
        if args.icolen is None:
            # User did not explicitly pin --icolen on this command line.
            # Reset it to "auto-default to ndcols" (None) rather than
            # leaving it at the mask's count -- so if an explicit --ndcols
            # below overrides the mask's count, icolen correctly re-tracks
            # THAT final value instead of staying pinned to the mask's
            # original count (which would silently mismatch ndcols).
            ctrl.icolen = None
        # else: user gave an explicit --icolen; leave ctrl.icolen as-is,
        # the explicit-overrides loop below will set it to args.icolen
        # regardless of what's here now.

    for attr, _key, _flag, _typ, _default, _help in PARAM_TABLE:
        val = getattr(args, attr)
        if val is not None:
            setattr(ctrl, attr, val)  # explicit CLI always wins last
    for name in ("enable_detailed_output", "enable_spec_output",
                 "enable_all_spec_output"):
        val = getattr(args, name)
        if val is not None:
            setattr(ctrl, name, val)

    for item in args.extra:
        if "=" not in item:
            print(f"lcm-control: --extra value {item!r} must be KEY=VALUE",
                  file=sys.stderr)
            return 1
        k, v = item.split("=", 1)
        ctrl.extra[k] = v

    try:
        write_control(ctrl, args.out_path)
    except ValueError as e:
        print(f"lcm-control: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
