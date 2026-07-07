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

## Domain terminology / file-naming conventions

**Don't assume software-testing naming conventions apply to MRS domain
filenames** — this has already caused two wrong inferences in this
project and is worth avoiding going forward.

- **`multi-voxel.met`** — **input**: the metabolite spectra data
  (`FILRAW`/`LRAW`). Not an output, despite the extension looking like
  it could be a results table — confirmed directly, don't re-guess this.
- **`multi-voxel.ref`** — **input**: the water reference spectra data
  (`FILH2O`/`LH2O`). "`.ref`" here means "water reference scan," the MRS
  domain meaning — not a "known-good expected output" the way
  `.ref`/`reference` usually means in software test suites.
- **`multi-voxel.csv`** — **output**: what the test run actually
  produces.
- **`test-reference-multi-voxel.csv`** — the real regression baseline
  for this test — **this** is what `multi-voxel.csv` should be diffed
  against, not `.met` or `.ref` (both of which are inputs, not
  baselines).
- **`control.file`** — input control-file configuration (drives what
  LCModel does for a given run — priors options, file paths, etc.).
- **`sim_se_csi_te16.basis`** — a basis-set file used for fitting (the
  metabolite basis spectra LCModel fits the data against).
- **When in doubt about what a file actually contains, ask before
  writing a comparison/diff command** — don't infer purely from
  filename pattern-matching against other tests in this repo or against
  general software-testing conventions. Check `git log` (is it a
  long-committed static file, or does it change every run?), check the
  Makefile rule that produces/consumes it, or just ask directly.

## Current work: multi-voxel performance project

This fork is a staged performance/architecture project on multi-voxel
(MRSI-style) runs, in priority order. **Do the items in order** — later
items depend on earlier ones being done safely.

1. **✅ DONE — Stop re-reading the whole input file per ring.** Confirmed
   via Claude Code trace (see "Per-voxel I/O caching" below for full
   detail): this was not a per-voxel re-read, it was a **per-ring**
   re-read. The scan pattern is an expanding square (Chebyshev distance
   from center); each grid position is fit exactly once, on the one ring
   matching its distance from center, but the inner triple loop
   (`idslic`/`idrow`/`idcol`) did not shrink to just that ring — it
   walked the *entire* grid on every ring pass, calling `DATAIN`→`MYDATA`
   (and hitting `REWIND LRAW`/`REWIND LH2O` at the top of each ring)
   unconditionally before the `skip_voxel` check discarded non-boundary
   positions. Net effect was O(rings × grid-size) reads for
   O(grid-size) actual fits. **Fixed** via an in-memory cache
   (`MODULE RAWCACHE`) built once in `check_zero_voxels`, consumed by
   `MYDATA` with a `raw_cache_ok .and. iaverg .le. 0` guard, falling back
   to the original disk-read path otherwise. All regression tests pass;
   full closure detail in "Goal #1 implementation plan" below.
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
three test targets (these are **file targets**, not phony names — invoke
them by the output file path, using `-B` to force a rebuild if needed):

- `make test_lcm/out.ps` — original single-voxel case, diffs against
  `out_ref_build.ps`.
- `make test_lcm/multi-voxel/multi-voxel.csv` — 100-voxel dataset,
  exercises the multi-voxel/multi-ring loop fully. Slower to run; use
  this before merging/finishing a change. Inputs: `multi-voxel.met`
  (metabolite spectra), `multi-voxel.ref` (water reference spectra),
  `control.file`, `sim_se_csi_te16.basis` (basis set). Output:
  `multi-voxel.csv` — **diff this against
  `test-reference-multi-voxel.csv`**, the actual regression baseline.
  Do not diff `.met`/`.ref` against anything — both are inputs, not
  baselines (see "Domain terminology" above).
- `make test_lcm/multi-voxel-10/multi-voxel.csv` — same as above but
  only 10 voxels. Much faster, so use this as the quick check while
  iterating on goal #1/#2/#3 work, then confirm with the full
  `multi-voxel` target before considering the change done.

**A change is only considered passing if all three tests pass** — the
single-voxel case guards against regressions in the non-multi-voxel path,
and the two multi-voxel cases guard against the actual behavior this
project is modifying. Don't treat `multi-voxel-10` passing as sufficient
on its own; it's a fast proxy for iteration, not a substitute for the
full `multi-voxel` run.

For any change to `LCModel.f`:

1. Build the binary (`make`).
2. Run all three test targets.
3. Diff each resulting output against its actual reference (`out.ps` vs
   `out_ref_build.ps`; each `multi-voxel.csv` vs its own
   `test-reference-multi-voxel.csv`).
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
`multi-voxel-10` cases already exist (`out_ref_build.ps` and each
directory's `test-reference-multi-voxel.csv`) — use those as the
baseline rather than generating new ones from scratch.

- For #2 and #3 (new prior-computation modes), generate and commit a new
  `test-reference-*.csv`-style baseline for each new option value the
  first time it's implemented correctly and reviewed, ideally against
  the `multi-voxel-10` dataset (fast) with a final confirmation against
  the full `multi-voxel` dataset. These become the regression baseline
  for that mode going forward — don't just eyeball a diff once and move
  on.

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
- **Fortran identifiers and keywords are case-insensitive** — `CALL
  MYDATA` and `call mydata` are the same statement. This codebase mixes
  capitalization inconsistently (e.g. `CALL DATAIN ()` vs.
  `call mydata ()` in nearby code). Any grep/search for call sites,
  variable usage, or COMMON-block references **must be case-insensitive**
  (`grep -i`), or it will silently miss real occurrences — this already
  caused one missed call site during the goal #1 investigation (a
  case-sensitive search for `CALL MYDATA` missed `call mydata ()` at
  line 1317, inside `average()`). Always re-verify a "confirmed
  complete" call-site or reference inventory with a case-insensitive
  pass before treating it as exhaustive.

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

### Goal #1 implementation plan — settled, two risks open before coding

Design is settled: build an in-memory cache once inside
`check_zero_voxels` (source/LCModel.f:1198-1248) — it already does one
full sequential pass through `LRAW` in `ivoxel` order and previously
discarded the data; now it captures it instead. `MYDATA`
(source/LCModel.f:2775-2966) branches on a new `raw_cache_ok` flag: if
true, reads from the cache; if false (allocation/read failure), falls
back to the original disk-read path unchanged — no functionality lost.
Uses `ALLOCATABLE` arrays sized to the actual run's `NUNFIL` and grid
dimensions (not the `MUNFIL`/`mvoxel` compiled maxima — a static array
at those maxima would be ~68 TB). `SCALE_RAW`/`SCALE_H2O` computation
stays exactly where it is (only their `TRAMP`/`VOLUME` inputs are
cached) because `fcalib` can be perturbed per-voxel by `LOADCH` and
`SCALE_RAW` reflects voxel-1's active `fcalib` specifically — moving
this earlier would silently change behavior. `LOADCH`'s separate
rewind-and-rescan pattern is confirmed untouched/out of scope for this
change. Checkpoint/resume is unaffected because `check_zero_voxels` runs
before the checkpoint-restore block and always covers the whole grid
regardless of resume point — same as today.

**Two things must be confirmed before/while implementing, not assumed:**

1. **`ivoxel` ordering parity — CONFIRMED, with a real bug caught and
   fixed in the plan before any code was written.** `ivoxel` is reset to
   0 immediately before the triple loop and incremented by exactly 1 at
   the innermost (`idcol`) level, with identical loop-variable order
   (`idslic → idrow → idcol`) and identical bounds at both the
   `check_zero_voxels` build site and the main ring loop's lookup site —
   these two line up 1:1.
   **However:** `average()` (called when `iaverg .ge. 1`, source/LCModel.f
   ~1317) calls `mydata()` **directly**, using its own counter `jvoxel`
   — never syncing COMMON `ivoxel`. `DATAIN`'s existing
   `if (iaverg .le. 0) CALL MYDATA()` gate only stops the *main ring
   loop* from calling `MYDATA` in this mode — it does not make `MYDATA`
   itself unreachable from `average()`'s own loop. Without an explicit
   guard, every call to `mydata()` inside `average()` would silently
   index `DATAT_cache(:, ivoxel)` using whatever stale `ivoxel` value
   `check_zero_voxels` left behind (its final value, since that loop
   runs to completion first) — feeding every averaging channel the same
   wrong cached voxel. **Required fix, now part of the plan:** the
   cache-lookup branch inside `MYDATA` must check `iaverg .le. 0` in
   addition to `raw_cache_ok` before trusting the cache — mirroring
   `DATAIN`'s own gate exactly — so any call reaching `MYDATA` via
   `average()` always takes the disk-read path, unconditionally.
   **Lesson for the rest of this project:** this was found only because
   an "it's out of scope, that path doesn't touch the cache" claim got
   pushed on rather than accepted. Before finalizing any guard/gate in
   this codebase, explicitly enumerate *every* call site of the
   function being gated (grep for all callers) rather than reasoning
   from the two or three sites already in view — `average()` was an
   undocumented third caller nobody had listed going in.
   **CLOSED.** Full call-site inventory confirmed (case-insensitive
   search — see Fortran-conventions note below): exactly three call
   sites total. `DATAIN` (source/LCModel.f:318, main ring loop) has no
   `iaverg` gate itself — it just delegates to `MYDATA`'s internal
   check. `MYDATA` has exactly two callers, which partition cleanly on
   `iaverg`: `DATAIN` (line 2511, fires only when `iaverg .le. 0`) and
   `average()` (line 1317, fires only when `iaverg .ge. 1`, since
   `average()` itself is only called when `iaverg .ge. 1`). The
   `iaverg .le. 0` guard added to `MYDATA`'s cache-lookup branch fully
   and correctly discriminates between the two — no third caller can
   slip past it.

2. **RAW vs. H2O record/voxel-count parity — CLOSED.** Not
   independently verified by file content — guaranteed only by
   construction, and that construction is confirmed sufficient. Both
   reads use the same global `NUNFIL` (no separate `NUNFIL_H2O` exists
   anywhere in the file), and `MYDATA` issues exactly one `LRAW` read
   and one `LH2O` read per call, in lockstep, using the same grid
   dimensions for both — so the two files are already consumed in
   strict 1:1 correspondence today. If the files were physically
   mismatched (wrong length or wrong voxel count), the existing code
   already fails fatally (`end=810` → `errmes` severity 4) rather than
   silently misaligning — there is no existing tolerance for a mismatch
   to preserve or replicate. **Conclusion: the H2O cache reuses the
   identical `(NUNFIL, ivoxel)` indexing as the RAW cache** — giving it
   an independent dimension would invent new mismatch-tolerance behavior
   beyond current scope, contradicting the "preserve existing behavior
   exactly" / "smallest change" constraints.
   **One incidental, documented (not hidden) behavior shift:** today, a
   deficient `LH2O` file's fatal error surfaces whenever the ring loop
   first reaches the affected grid position (first ring pass). With the
   new dedicated `LH2O` pass in `check_zero_voxels`, the same fatal
   error instead surfaces during cache-build, before the ring loop
   starts — same ultimate failure, same voxel identified, only the
   wall-clock timing of the error changes. Flag this explicitly in the
   regression review even though it's benign — it's a real (if trivial)
   behavior difference, not "byte-for-byte identical."

**Both risk items are now closed.** The plan has been through: ring-
structure diagnosis, file-format constraint check, prior-state safety
confirmation, checkpoint/resume safety confirmation, a real stale-cache
bug caught and fixed (`average()`'s direct `mydata()` calls), a complete
call-site inventory, and RAW/H2O parity. This is a reasonable point to
move from Plan mode into writing the actual code — implement, then run
the full regression protocol (all three test targets, plus the two
manual checks: a run with `FILH2O` set, and a checkpoint/resume run).

**First build is the go/no-go checkpoint for whether `gfortran
-std=legacy -O3` accepts the `ALLOCATABLE` extension in this fixed-form
file.** — **Result: `ALLOCATABLE` cannot be a `COMMON`-block member**
(gfortran: "COMMON attribute conflicts with ALLOCATABLE attribute"),
confirmed via an isolated standalone test before touching the real file.
**Design updated accordingly:** the cache lives in a `MODULE RAWCACHE`
instead of a `/BLRAWCACHE/` COMMON block in `lcmodel.inc`. This is a
pure storage-mechanism swap — none of the logic-level guarantees already
vetted (ivoxel ordering parity, the `iaverg .le. 0` guard, RAW/H2O
shared indexing) depended on COMMON specifically, so none of that needed
re-verifying. `check_zero_voxels` (builder) and `MYDATA` (reader) each
need `USE RAWCACHE` as their first statement instead of picking up these
variables from `lcmodel.inc`.
Module-persistence semantics were **empirically confirmed**, not
assumed: the module-level initializer (`raw_cache_ok = .false.`) fires
exactly once at program start; `USE`-association across
`check_zero_voxels` (called once) and `MYDATA` (called once per voxel
per ring, potentially hundreds of times) shares the same persistent
storage location — no reset, no re-triggering per `USE`. Verified with a
standalone test: one builder call, five-plus reader calls including one
from a call site that never invoked the builder, value persisted
correctly throughout.
**Two items still open on this sub-design, not yet confirmed:**
- Module must be placed lexically before any subroutine that `USE`s it
  in the single-file compilation — confirm actual placement relative to
  `PROGRAM LCMODL`/`check_zero_voxels`/`MYDATA`.
- The cached character lengths (`fmtdat_raw_cache*80`, `id_raw_cache*20`)
  are hardcoded literals that must independently track `MCHFMT`/`MCHID`
  in `lcmodel.inc` — two numbers kept in sync only by convention, not by
  the compiler. Prefer referencing `MCHFMT`/`MCHID` directly if feasible,
  rather than duplicating the literals (drift risk if either changes
  later without updating the other).

**Design status: fully vetted, diffs reviewed as text, ready to apply to
disk.** Every risk raised during review is closed: ring-structure
diagnosis, file-format constraint, prior-state safety, checkpoint/resume
safety, the `average()` stale-cache bug (caught and fixed via an
explicit `iaverg .le. 0` guard on every cache-consumption site), complete
call-site inventory, RAW/H2O parity, the COMMON→MODULE pivot (with
empirical persistence confirmation), module placement, `MCHFMT`/`MCHID`
sharing via a new `lcmodel_params.inc`, the SEQPAR-probe/`bascal` gating
equivalence (confirmed identical, not just "probably fine" — the
original authors' own comments in `MYBASI` document that
`echot_raw`/`seq_raw` are expected unset during `bascal` runs), and the
`ivoxel` COMMON-promotion mechanics (no textual change needed — `ivoxel`
was already a bare implicitly-typed local in both `PROGRAM LCMODL` and
`check_zero_voxels`, so adding it to `/BLINT/` via `lcmodel.inc` is
sufficient). Nothing left to verify at the design level — remaining
checkpoints are execution: the build gate, the regression protocol, and
the two manual checks below.

**GOAL #1 STATUS: DONE.** Implemented, correct, and regression-clean on
all three test targets plus all required manual checks. Full closure:

- **All three automated regression targets pass.** `out.ps` matches
  `out_ref_build.ps`; both `multi-voxel.csv` outputs match their
  respective `test-reference-multi-voxel.csv` exactly.
- **`FILH2O` coverage — confirmed, already covered by existing
  fixtures.** Both `test_lcm/multi-voxel/control.file` and
  `test_lcm/multi-voxel-10/control.file` already set
  `filh2o = 'multi-voxel.ref'`, so the regression runs that passed
  already exercised the H2O caching path — no separate fixture was
  needed.
- **Dev-aid fallback check — confirmed.** Forced `raw_cache_ok = .false.`
  after `check_zero_voxels()` and reran all three regression targets:
  all matched the pre-change reference exactly, confirming the
  disk-read fallback path in `MYDATA` was not altered while editing it.
- **A real bug was caught and fixed by this fallback check, not just
  confirmed clean.** `check_zero_voxels`'s own `NMID` read (used to get
  `FMTDAT` for the zero-check pass) never initialized
  `tramp`/`volume`/`id`/`FMTDAT`/`BRUKER`/`SEQACQ` before reading —
  harmless in the original code (those values were discarded after use)
  but now cached and consumed authoritatively by `MYDATA`. Uninitialized
  `BRUKER`/`SEQACQ` (stack garbage, not guaranteed `.FALSE.`) caused an
  unwanted conjugate/`SEQTOT` transform on `DATAT`, producing a real fit
  regression (different noise-SD estimate, different plotted Y-scale)
  despite the raw `DATAT` values themselves being provably identical
  between cache and disk paths. **Fix applied:** the same
  `tramp=1./volume=1./id=' '/FMTDAT=' '/BRUKER=.FALSE./SEQACQ=.FALSE.`
  initialization `MYDATA` always had was added immediately before
  `check_zero_voxels`'s `NMID` read. Found by comparing full-array
  `DATAT` checksums between cache and disk paths (identical) then
  checking what else the cache branch sets — not by guessing from the
  symptom. **Lesson for future work in this file:** when duplicating a
  read/parse pattern from one subroutine into another (as this project
  does repeatedly), don't just copy the read statement — copy the
  initialization that precedes it too, even if the original subroutine
  discarded the values and "didn't need" the initialization at the time.
- **Performance measured honestly — real but not strongly visible on
  the bundled fixtures.** `multi-voxel` (100 voxels, `NUNFIL=1024`):
  73.8s cache-enabled vs. 77.2s cache-disabled (~4% faster).
  `multi-voxel-10`: no measurable difference (6.2–7.3s either way,
  run-to-run noise exceeds the effect). At this scale, per-voxel fit
  computation (`MYBASI`/`TWOREG`/regularized search) dominates
  wall-clock time; the O(rings×grid)→O(grid) read-count reduction is
  architecturally real but too small a fraction of total runtime to
  show clearly on these particular fixtures, which are sized for
  regression correctness, not performance benchmarking.
- **✅ Real-world large-dataset validation — confirmed the performance
  benefit at production scale.** A real 9343-voxel dataset: new
  (cached) version ~1.6 hours; old (uncached) version projected
  8–10 hours (still running at time of comparison) — roughly a 5–6×
  speedup, exactly the scaling behavior predicted (the benefit grows
  with grid size, which is why it was too small to see on the 100-voxel
  fixture). Output confirmed identical between old and new versions for
  the first ~3600 voxels completed by the old run at comparison time.
  **Not yet confirmed for the full 9343 voxels** — the old run hadn't
  finished. Worth a final check once it completes, to convert "matches
  on the first ~3600" into "matches on all 9343," though given the
  fixture-level regression tests already passed bit-for-bit and the
  fallback path was separately verified, this is confirmation rather
  than a live open risk.
  This closes the original motivating complaint ("very slow" multi-voxel
  runs) — the fix now has both correctness and performance validated at
  a scale representative of real use, not just on small regression
  fixtures.
- **Checkpoint/resume run — deliberately SKIPPED, not overlooked.**
  Decision: the checkpoint/resume mechanism (`lcsi_sav_1`/`lcsi_sav_2`,
  `ioffset_current_in`/`nvoxels_done_in`) is not used by this project's
  actual workflows, and hand-constructing synthetic mid-run checkpoint
  state was judged higher-risk than valuable — a wrong synthetic state
  would test "does the program tolerate a file I made up," not "does
  resume actually work," which could give false confidence either way.
  The plan's reasoning for why this should be safe (see "Known
  interaction to protect explicitly," above) still stands as
  **reasoning-only, not empirically confirmed**. If checkpoint/resume
  ever becomes relevant to this project's actual use, the right way to
  test it is to let a real run get interrupted naturally (kill a
  `multi-voxel-10` run mid-scan, after checking where checkpoints
  actually get written) and resume it for real, diffing the final
  output against an uninterrupted run of the same dataset — not by
  hand-building checkpoint files.
- **No user-facing on/off toggle exists for the cache**, by deliberate
  choice — it's an internals-only optimization intended to be
  bit-for-bit transparent, not a new behavior needing a default-off
  option like goals #2/#3. The only way to disable it is editing source
  (`raw_cache_ok = .false.`) and recompiling, which is what the dev-aid
  check above used. Revisit this only if a real debugging or
  memory-constrained-system need for a runtime toggle comes up later.

**Before committing:** clean up the working tree first —
`git status` currently shows several untracked files that don't belong
in the repo: `rawcache.mod` (compiled module artifact), stray backup
files (`.~CLAUDE.md`, `binaries/linux/bak.lcmodel`,
`source/dng.LCModel.f`/`dng2.LCModel.f` if not intentionally kept
elsewhere), and unrelated screenshots. Add a `.gitignore` covering
`*.mod` and editor backup patterns before the next commit so this
doesn't recur.

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
