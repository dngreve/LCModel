# CLAUDE.md

Project context and rules for Claude Code when working in this repository.

<!-- Maintainer note: keep this file short and stable. Deep, topic-specific
     detail belongs in docs/ or .claude/rules/, not here. Update the
     "Current work" section every time the goal of the fork changes. -->

## What this project is

This repo (`dngreve/LCModel`, forked from `schorschinho/LCModel`) modifies
**LCModel**, Stephen Provencher's linear-combination-modeling program for
quantitative analysis of in-vivo magnetic resonance spectroscopy (MRS) data.
LCModel fits a spectrum as a linear combination of metabolite basis spectra
plus a baseline, and is widely used as a reference/gold-standard analysis
tool in the MRS research community.

The codebase is **legacy Fortran 77**, fixed-form source, largely a single
file (`source/LCModel.f`), originally distributed as a closed-source binary
and later open-sourced. Treat it as scientific/numerical software: the
priority is *numerical correctness and reproducibility*, not code aesthetics.

## Current work: multi-voxel performance project

This fork is a staged performance/architecture project on multi-voxel
(MRSI-style) runs, in priority order. **Do the items in order** — later
items depend on earlier ones being done safely.

1. **Stop re-reading the whole input file per ring.** Confirmed via Claude
   Code trace (see "Per-voxel I/O caching" below for full detail): this is
   not a per-voxel re-read, it's a **per-ring** re-read. The scan pattern
   is an expanding square (Chebyshev distance from center); each grid
   position is fit exactly once, on the one ring matching its distance
   from center, but the inner triple loop (`idslic`/`idrow`/`idcol`) does
   not shrink to just that ring — it walks the *entire* grid on every
   ring pass, calling `DATAIN`→`MYDATA` (and hitting `REWIND LRAW` /
   `REWIND LH2O` at the top of each ring) unconditionally before the
   `skip_voxel` check discards non-boundary positions. Net effect:
   O(rings × grid-size) reads for O(grid-size) actual fits. Fix direction
   is settled (see below) — implementation is not.
2. **Option: derive priors from the first voxel only.** Add a control-file
   option so that, instead of updating priors after every voxel (current
   default behavior), priors are computed once from voxel 1 and reused
   for all subsequent voxels. Must be opt-in — default behavior (rolling
   per-voxel prior updates) has to stay bit-for-bit identical unless the
   new option is explicitly set.
3. **Option: derive priors from a voxel in a separate file.** Extend #2
   so the reference voxel for priors can come from a different input
   file than the one being fit, not just voxel 1 of the current file.

OpenMP parallelization is a possible future goal but is explicitly
**out of scope for now** — don't introduce `!$OMP` directives, threading,
or thread-safety refactors unless asked. Keeping this out of scope
simplifies #1–#3 considerably (no need to reason about concurrent access
while restructuring the I/O and prior logic).

- Constraints: preserve the existing default control-file behavior
  exactly (regression tests must still match `out_ref_build.ps` byte-for-
  byte where no new option is invoked); new options must be additive,
  off by default.
- Definition of done per item: regression protocol below passes for the
  default path, plus a new reference output is captured and reviewed for
  each new opt-in mode (see "New reference outputs" below).

## Repository layout

```
source/       LCModel.f (main F77 source), plus any support files
binaries/     precompiled binaries per OS/arch — build artifacts, not source
test_lcm/     control files + reference outputs for regression testing
              (single-voxel, multi-voxel [100 voxels], multi-voxel-10 [10 voxels])
Makefile      build target plus `multi-voxel` / `multi-voxel-10` test targets
docs/         LCModel manual (PDF) and relevant papers — see "Reference material"
```

## Build (Rocky 9)

`gfortran` is already installed and working in this environment. Just run:

```bash
make            # uses ./Makefile, produces the lcmodel binary
```

If you ever need to invoke the compiler directly, match the flags in the
Makefile rather than inventing new ones — legacy fixed-form code is
sensitive to floating-point flags. In particular, do not add
`-ffast-math`, `-Ofast`, or other flags that relax IEEE floating-point
semantics — they can silently change fit results.

## Testing / regression protocol (critical — do this for every change)

LCModel's correctness is defined by numerical agreement with reference
output, not by unit tests in the usual sense. The Makefile provides
three test targets:

- `make test` (or whatever the default target is) — original single-voxel
  case, diffs against `out_ref_build.ps`.
- `make multi-voxel` — 100-voxel dataset, exercises the multi-voxel loop
  fully. Slower to run; use this before merging/finishing a change.
- `make multi-voxel-10` — same as above but only 10 voxels. Much faster,
  so use this as the quick check while iterating on goal #1/#2/#3 work,
  then confirm with the full `multi-voxel` target before considering the
  change done.

**A change is only considered passing if all three tests pass** — the
single-voxel case guards against regressions in the non-multi-voxel path,
and the two multi-voxel cases guard against the actual behavior this
project is modifying. Don't treat `multi-voxel-10` passing as sufficient
on its own; it's a fast proxy for iteration, not a substitute for the
full `multi-voxel` run.

For any change to `LCModel.f`:

1. Build the binary (`make`).
2. Run all three test targets.
3. Diff each resulting `.ps` output against its reference.
   - Expected differences: build date, version string (`6.3-N` vs `6.3-R`
     style tags).
   - Any difference in fitted concentrations, CRLBs, SNR, or plotted
     spectra is a regression unless the change was *intended* to alter
     the algorithm — flag it explicitly, don't silently accept it.
4. If the change is expected to alter numerical output (e.g. an
   intentional algorithm change), say so up front and show a before/after
   comparison rather than just a passing/failing diff.

If more test datasets exist or are added, extend this protocol to cover
them — don't rely on the three bundled tests as permanently sufficient
coverage once the project grows.

### New reference outputs for this project

Reference outputs for the single-voxel, `multi-voxel`, and
`multi-voxel-10` cases already exist and ship with their respective make
targets — use those as the baseline rather than generating new ones from
scratch.

- For #2 and #3 (new prior-computation modes), generate and commit a new
  reference output for each new option value the first time it's
  implemented correctly and reviewed, ideally against the
  `multi-voxel-10` dataset (fast) with a final confirmation against the
  full `multi-voxel` dataset. These become the regression baseline for
  that mode going forward — don't just eyeball a diff once and move on.

## Fortran-specific conventions

- Fixed-form F77: columns 1–5 are labels, column 6 is continuation,
  code starts at column 7, and columns beyond ~72 are ignored by some
  compilers. Preserve this layout; don't reflow lines.
- Don't convert COMMON blocks, GOTO-driven control flow, or implicit
  typing to "modern" equivalents as a drive-by cleanup. Refactor only
  when it's the actual point of the task, and do it in isolated,
  reviewable commits with regression tests passing before and after.
- Preserve the order of floating-point operations in numerical loops
  (summations, matrix ops). Reordering changes rounding and can shift
  fit results even when "mathematically equivalent."
- When adding new code, F77-style fixed-form is fine for consistency
  with the surrounding file; don't mix in free-form Fortran in the same
  source file without a good reason.

### Per-ring I/O caching (goal #1) — findings confirmed by trace

Full investigation was done via Claude Code in Plan mode before any code
was touched. Confirmed facts, so future sessions don't need to re-derive
these:

**Structure.** The scan is a ring-based expanding square, not a plain
per-voxel loop. `REWIND LRAW`/`REWIND LH2O` (~line 268) sit at the top of
a `DO ioffset = ...` loop wrapping a triple nested loop over
`(idslic, idrow, idcol)` (lines 272–380). Gated by `.not. voxel1`
(`/BLLOG/` COMMON flag, true only before voxel 1 of the whole run) — no
gating on `NUNFIL`/`IAVERG`/etc. `LH2O` additionally requires
`FILH2O .ne. ' '`.

**Cost shape.** Each grid position is fit exactly once — on the single
ring matching its Chebyshev distance (`max(|idrow-irow_center|,
|idcol-icol_center|)`) from center — but the inner triple loop's bounds
are fixed (`1..ndslic`/`1..ndrows`/`1..ndcols`) and don't shrink per
ring. So `DATAIN`→`MYDATA` (lines ~2492, ~2775) is called, and the file
is read, for *every* grid position on *every* ring pass, even positions
that get `skip_voxel`'d immediately after (check is at line ~321, after
the read already happened). Net: O(rings × grid-size) reads for
O(grid-size) actual fits — this, not per-voxel repetition, is the real
cost driver.

**File format constraint.** `LRAW`/`LH2O` are opened with no `ACCESS=`/
`FORM=` clause (line ~1215) → sequential formatted (text), not
direct-access binary. The per-voxel record format (`FMTDAT`) is itself
supplied at runtime via the `NMID` namelist and cached as `FMTDAT_RAW`.
Every voxel's data block is a deterministic, uniform number of text
lines, but Fortran can't seek an arbitrary offset into sequential
formatted access the way it could with `ACCESS='DIRECT'`. **This rules
out a seek/reopen-with-computed-offset fix** — reopening with direct
access against a runtime-determined format string would be a much
riskier change than the alternative below.

**Settled fix direction:** read the file once, sequentially, in natural
voxel order, caching each voxel's `DATAT`/`H2OT` (plus per-voxel
header/format info) in memory; the ring loop indexes into that cache
instead of calling `DATAIN` again. This is the smaller, safer change
given the codebase's conventions — implementation plan (data structure,
call-site changes) is not yet written; do that as its own Plan-mode pass
before touching code.

**Known interaction to protect explicitly:** there's an existing
checkpoint/resume mechanism (`ioffset_current_in`/`nvoxels_done_in`,
read from units 12/13 around lines 185–201) that lets a run resume
mid-scan at an arbitrary ring, not just ring 0. Any in-memory cache must
still work correctly when a run resumes from a checkpoint — don't build
a cache that implicitly assumes the scan always starts at position 1.

**COMMON-block state: confirmed safe, not just assumed.** The per-voxel
data buffers (`DATAT`, `H2OT`, `DATAF`, `H2OF` in `/BLCPLX/`) are freshly
overwritten every call — nothing to protect. The actual cross-voxel
prior state (`DEGPPM`, `SDDEGZ`, etc. in `/BLREAL/`) is threaded through
a completely separate mechanism — `update_priors()` (line ~1582,
subroutine-local `SAVE` arrays) and `restore_settings()` (line ~1465) —
that never touches `LRAW`/`LH2O`. Removing the rewind/re-read does not
risk resetting prior/accumulator state.

**Out-of-scope but related — decide before or during implementation:**
- `check_zero_voxels` does its own separate full pre-scan of the raw
  file before the ring loop even starts — a second read pass over the
  same data.
- `LOADCH` (line ~2117) has a structurally identical rewind-and-rescan
  pattern for the control-file scratch unit, executed once per analyzed
  voxel (cheaper than the ring re-read, but same root anti-pattern).

Decide explicitly whether goal #1's fix should cover only
`LRAW`/`LH2O`, or all three re-read sites in one pass, since they share
a root cause. Don't silently expand scope without noting the decision
here.

## Reference material

The LCModel manual and related MRS methodology papers live under `docs/`
(add them there if not already present: manual PDF from
s-provencher.com/pub/LCModel/manual/manual.pdf, plus any papers you're
working from). Point Claude at specific sections/pages rather than
pasting large excerpts into this file — ask Claude to read the relevant
PDF pages directly when a task depends on the manual's specification
(e.g. control-file parameter semantics, basis-set file format, output
format).

## Licensing

LCModel is distributed under Provencher's own license terms (see
`LICENSE.md`), not a standard OSS license. Don't assume MIT/BSD-style
freedoms — check `LICENSE.md` before adding redistribution, packaging,
or dependency-bundling logic.

## Scope limits — ask before doing these

- Don't rewrite large sections of `LCModel.f` in a single pass; prefer
  small, independently testable diffs.
- Don't remove or "clean up" numerical edge-case handling
  (underflow/overflow guards, `ffpe-summary` suppressions, etc.) without
  understanding why it's there.
- Don't add new external dependencies unless the task explicitly calls
  for it — LCModel's build is intentionally minimal (gfortran + make).
