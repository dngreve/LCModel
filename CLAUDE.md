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

## Guiding principle: control file = shared config, CLI = per-case specifics

The control file is meant to function as a **study-wide configuration**
applied uniformly across many datasets/subjects — the whole point is
that using the same control file guarantees every dataset in a study is
analyzed the same way. **Input files, output files, and other genuinely
per-case specifics should be movable to the command line**, not baked
into the control file, so the same control file can be reused unchanged
across an entire study rather than hand-edited per subject.

This is a standing design principle, not a one-off feature request —
it informs goal #3 (CLI overrides for `-i`/`-csv`) and should inform any
future per-case settings too, even ones not yet identified. When in
doubt about whether a new option belongs in the control file or on the
CLI, ask: does this vary per-dataset (CLI), or should it stay identical
across an entire study for consistency (control file)?

**Mechanism precedent**, established in goal #2 and reused in goal #3:
CLI flags parsed before the control file is read, with a companion
`_set_by_cli` flag per overridable setting, and an active reapplication
step after `MYCONT()` returns for anything CLI-set — not a passive
"set before and hope the control file's NAMELIST read leaves it alone."

**Not yet needed, but worth revisiting if it comes up**: if a fourth or
fifth per-case CLI override gets added on top of goal #3's `-i`/`-csv`,
the current approach (hand-wiring a `_set_by_cli` companion flag and a
reapply step individually for each) will start feeling repetitive.
Don't build a generalized overridable-settings table now — that would
be solving a problem that doesn't exist yet — but that's the point to
reconsider a more systematic mechanism instead of continuing to hand-copy
the same pattern.

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
items depend on earlier ones being done safely. (Note: the original
goal #3 — separate-file prior source — was abandoned; #3 now refers to
a different, unrelated goal, CLI overrides for input/output files. See
below.)

1. **✅ DONE — Stop re-reading the whole input file per ring** (see
   "Per-ring I/O caching" below for full detail, including a resolved
   incident worth reading: a deterministic `NaN`/crash was found and
   root-caused post-implementation — a `COMMON`-block inconsistency
   between `lcmodel.inc` and `LCModel.f`, not a flaw in the caching
   design itself — fixed and independently re-verified). Original
   diagnosis: this was
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
2. **Options: control prior-update behavior via command-line flags.**
   **✅ `-no-prior`: DONE (committed). `-first-prior`: REMOVED, replaced
   by `-prior-file`/`-save-prior` — see "Goal #2b" below, STATUS: DONE,
   not yet committed.** The detail immediately below describes the
   *original* `-first-prior` design, kept for historical record only —
   `-first-prior` was dropped (it inherited a real semantic bug from
   the ring-scan order, plus an administrative burden — a phantom
   voxel that had to be tracked and excluded from CSV output).
   Original design (supersedes the original control-file-keyword idea —
   deliberate mechanism choice, see rationale below):
   - New integer state variable **`UsePrior`** (not boolean — three
     modes now, a fourth reserved for goal #3):
     - `1` = **default** — update priors after every voxel (current,
       unchanged behavior; must stay bit-for-bit identical when no flag
       is given).
     - `2` = **`-no-prior`** — never call `update_priors()`.
     - `3` = **`-first-prior` — DEPRECATED, being removed.** Originally:
       call `update_priors()` once, at voxel 1 only, then never again.
       **Problem discovered:** under the ring-scan order (see goal #1),
       "voxel 1" for this purpose actually meant "the first voxel
       *analyzed*," which is the scan-center voxel, not `ivoxel=1`
       (the first file record) — the plain-language spec and the
       actual mechanism silently disagreed. Separately, freezing a
       prior from an in-map voxel created an administrative burden:
       that voxel's row had to be identified and excluded from the
       output CSV, since it isn't really "a voxel from the map." See
       goal #2b below for the replacement.
     - `4` = now used by **goal #2b** (external fixed prior via
       `-prior-file`), not the original goal #3 (separate-file voxel
       source, itself abandoned — see below).
   - Use named `PARAMETER` constants for these states (e.g.
     `USEPRIOR_DEFAULT=1`, `USEPRIOR_NONE=2`, `USEPRIOR_FIRST=3`
     [being removed], `USEPRIOR_FILE=4` [now goal #2b]), not bare
     integer literals scattered through the code.
   - **If both `-no-prior` and `-prior-file` are given, error out** —
     mutually exclusive, not "last one wins" or silently combined (same
     rule that applied to `-no-prior`/`-first-prior`, carried forward).
   - **Mechanism: real command-line flags, parsed early in
     `PROGRAM LCMODL` before the control file is read** — not a
     control-file keyword. Deliberate choice: pure control-file
     configuration is a real usability pain point, and CLI flags are
     wanted as an extensible pattern for future options too.
     `UsePrior` itself is a plain `INTEGER` (no `ALLOCATABLE` conflict
     like `RAWCACHE` had), so it can live directly in an existing
     `COMMON` block (e.g. `/BLINT/` in `lcmodel.inc`) — no module needed
     for this one.
   - **`GET_COMMAND_ARGUMENT` compile/runtime test — CONFIRMED, with one
     real wrinkle found and fixed.** `gfortran -std=legacy -O3` accepts
     `COMMAND_ARGUMENT_COUNT()`/`GET_COMMAND_ARGUMENT` cleanly in
     fixed-form F77 — same pattern as goal #1's `ALLOCATABLE`/`MODULE`
     check (modern intrinsics accepted under `-std=legacy` without
     special handling). Runtime-tested all six cases: no args (default),
     `-no-prior`, `-first-prior`, both together (error), an unrecognized
     flag (error), and the same flag given twice (correctly *not* an
     error — not a real conflict). **Wrinkle:** a bare `STOP 1` prints
     gfortran's own `STOP 1` runtime annotation before the custom
     `WRITE` error message appears (stdout buffering vs. the runtime's
     termination message) — reads confusingly out of order. Bare `STOP`
     avoids the noise but silently returns exit code 0, which is wrong
     for a CLI validation error (a caller/script needs a nonzero exit to
     detect failure). **Confirmed fix: explicit `FLUSH(6)` before
     `STOP 1`** — correct message ordering and correct nonzero exit
     code, confirmed together. Use this `WRITE`+`FLUSH(6)`+`STOP 1`
     pattern for both the mutual-exclusivity and unrecognized-flag
     errors. Nothing touched in the real source tree during this test —
     scratchpad only.
   - **`update_priors()` call site — CONFIRMED.** `source/LCModel.f:407`,
     `if (.not.single_voxel) call update_priors ()`, gated by
     `single_voxel` (set at line 271:
     `single_voxel = max0(ndslic,ndrows,ndcols) .eq. 1`). The earlier
     CLAUDE.md hint guessed this exact line/condition but hadn't
     verified it — turned out correct, confirmed directly rather than
     trusted.
   - **A real `restore_settings()` bug found and fixed — the most
     consequential finding for this goal.** `restore_settings()`
     (line 324, once per voxel) decides keep-vs-restore `DEGPPM`/
     `DEGZER` via a heuristic: near-zero means "restore saved defaults,"
     non-zero means "`update_priors()` just ran, keep it." But this
     signal is a one-shot echo, not a persistent flag — `STARTV` always
     zeroes `DEGPPM`/`DEGZER` again after each voxel's `PHASTA()`
     consumes them (confirmed sequence: `restore_settings()` [324] →
     `STARTV`→`PHASTA` consumption of `EXDEGZ`/`EXDEGP` [line ~5911→
     ~6744] → zeroing [~6110] → `update_priors()` [407] writes the next
     voxel's values). Under default mode this is fine (`update_priors()`
     refreshes the signal every voxel). **Under `-first-prior`, it
     silently breaks from voxel 3 onward**: voxel 2 correctly gets
     voxel-1's prior (the one-shot echo is still fresh), but voxel 2's
     own `STARTV` zeroes it again, and since `update_priors()` never
     fires a second time, voxel 3's `restore_settings()` sees near-zero
     and silently reverts to control-file defaults — the opposite of
     "frozen after voxel 1," for every voxel from 3 onward. **Fix:** a
     new `first_prior_done` flag, set `.true.` the moment
     `update_priors()` fires under `-first-prior` (only write site, at
     line 407's new branch), wraps `restore_settings()`'s check so it's
     skipped entirely once true — forcing "always keep" from that point
     on. No-op for default/`-no-prior` (`first_prior_done` always
     `.false.` there) — byte-identical behavior preserved for those
     modes. **Verified orthogonal to `STARTV`'s internal `ipass=1`/
     `ipass=2` mechanism** (its own private `_sav_startv` locals,
     unrelated — `restore_settings()` runs exactly once per voxel,
     before either `STARTV` call, never interacting with `STARTV`'s
     internal repeat-pass handling).
   - **`first_prior_done` persistence — verified, not assumed from
     `raw_cache_ok`'s precedent.** Unlike `raw_cache_ok`
     (module+`ALLOCATABLE`, genuine untested-territory needing an
     empirical compile test), this is a plain `LOGICAL` in `COMMON` —
     the same mechanism `voxel1`/`lraw_at_top`/`skip_voxel` already use
     successfully. Still checked directly rather than assumed safe by
     analogy: no name collision (`first_prior_done`/`useprior`,
     case-insensitive grep, zero existing hits); initialization follows
     `voxel1`'s convention (explicit runtime assignment near the top of
     `PROGRAM LCMODL`, not `BLOCK DATA`); single write site (line 407's
     new branch only); and `average()`'s `iaverg .ge. 1` path can't
     reach that write site at all, since `average()` collapses
     `ndslic`/`ndrows`/`ndcols` to 1 before `single_voxel` is computed —
     making `single_voxel` true and skipping the entire gated block,
     the same guarantee goal #1 relied on for the cache guard.
   - **`ERRMES` confirmed unsuitable for CLI-parse-time errors — two
     independent reasons.** (1) `LPRINT` defaults to 0 until `MYCONT()`
     reads it from the control file; `ERRMES`'s diagnostic only prints
     `IF (LPRINT .GT. 0)`, so called this early the message would be
     silently swallowed. (2) More fundamentally: `EXITPS`'s actual
     `STOP` statement is commented out project-wide
     (`source/LCModel.f:~11061`, deliberately, so one bad voxel doesn't
     kill a whole multi-voxel batch) — a "fatal" `ERRMES` call today
     just prints "Error in voxel but continuing" and carries on. For a
     CLI-parsing error that's exactly backwards; a hard stop is needed
     before `MYCONT()`/the rest of the pipeline runs on unparsed or
     contradictory flags. **Replacement:** plain `WRITE`+real `STOP`,
     for both the mutual-exclusivity case and any unrecognized flag
     (confirmed: unrecognized flags should hard-stop, not be silently
     ignored, to catch typos immediately).
   - **No regression baseline exists yet for the non-default modes.**
     The default path (`UsePrior=1`, no flags) must still match all
     three existing regression targets exactly. For `-no-prior` and
     `-first-prior`, there is nothing to diff against yet — the first
     correct, manually-reviewed run becomes the new committed baseline
     (same pattern as "New reference outputs," below). Don't treat this
     as done until that manual review has actually happened — it's not
     just an implementation task, it's an open verification task too.
   - **Forward-looking hook, deferred but cheap to add now:** there's no
     control-file equivalent of `UsePrior` today, so there's no
     CLI-vs-control-file override conflict to resolve yet — don't build
     a general override/precedence mechanism now, it would be designing
     for a hypothetical. **But** since CLI flags are parsed before the
     control file is read, if a control-file `NAMELIST` field for this
     setting is ever added later, an unconditional `NAMELIST` read would
     silently clobber whatever the CLI flag set — an accidental
     parsing-order bug, not a deliberate precedence choice. Cheap
     insurance to add now, while this code is already being written:
     also set a companion flag (e.g. `UsePrior_set_by_cli = .true.`)
     recording that a CLI flag was explicitly given, not just what value
     it produced. Costs one extra logical variable; gives a future
     override rule ("if set via CLI, don't let the control file
     overwrite it") something to check, instead of having to retrofit
     "was this explicitly set via CLI" after the fact.
   - **Design status: fully vetted, diffs reviewed, ready to apply.**
     Every item confirmed, not assumed:
     - Real call sites for both `update_priors()` (407) and
       `restore_settings()` (324) — confirmed exactly one caller each,
       no hidden third path (`average()` was a surprise third caller for
       `MYDATA` in goal #1; swept for here too, came back clean).
     - **`restore_settings()`'s freeze fix corrected during review** —
       the first version only skipped the reset-to-defaults check, which
       left `DEGPPM`/`DEGZER` at whatever `STARTV`'s per-voxel zeroing
       had already set them to (i.e., silently 0 from voxel 3 onward,
       not the frozen value — "worked" for voxel 2 only by coincidence
       of timing). **Corrected:** capture `degppm_frozen`/
       `degzer_frozen` at the moment `update_priors()` fires under
       `-first-prior`, then have `restore_settings()` **actively
       reassign** `degppm`/`degzer` from those frozen values every voxel
       from 2 onward, not just skip a check. Default/`-no-prior` paths
       remain byte-identical (`first_prior_done` never true there).
     - `ERRMES` confirmed unsuitable; `WRITE`+`FLUSH(6)`+`STOP 1`
       confirmed as the working replacement, including a
       message-ordering wrinkle found and fixed (bare `STOP 1` prints
       gfortran's own annotation out of order; `FLUSH(6)` first fixes
       it while keeping the required nonzero exit code).
     - `first_prior_done`/`UsePrior_set_by_cli` persistence verified
       directly (no name collision, correct initialization convention,
       single write site, unreachable from `average()`'s path).
     - `USEPRIOR_*` constants corrected to explicit `INTEGER
       PARAMETER`s — every existing constant in `lcmodel_params.inc`
       starts with `M`/`m` (implicit-`INTEGER` range); `USEPRIOR_*`
       would have been the first to break that convention and rely on
       implicit-`REAL`-but-exactly-representable luck instead.
     - Unit 6/stdout consistency and pre-`STOP` cleanliness (no open
       files/allocated state before the CLI-parsing block, confirmed by
       a full line-by-line read of `PROGRAM LCMODL` through `MYCONT()`)
       both confirmed.
     Nothing left unverified. Remaining work is execution: apply the
     diffs, build, run the existing regression suite (must stay
     unchanged for the default no-flags path), then manually run and
     review `-no-prior`/`-first-prior` to establish their first
     baselines.
3. **ABANDONED — deriving priors from a voxel in a separate file
   (`UsePrior=4`).** Dropped in favor of the goal below. `UsePrior=4` is
   no longer reserved; the slot is free if ever needed again, but
   nothing currently plans to use it.
3. **(new) CLI overrides for input/output files, plus a mandatory-input
   check.** Usage target:
   `lcmodel -i met.lcm h2o.lcm -csv out.csv < lcm.control`
   - **`-i <met-file> <h2o-file>`** — overrides `FILRAW`/`FILH2O` (the
     existing control-file variables). **Both paths are mandatory when
     `-i` is given** — always consumes exactly the next two arguments,
     no optional-second-path parsing needed (decided explicitly, not
     assumed).
   - **`-csv <output-file>`** — overrides the control-file CSV-output
     filename variable (exact name TBD — confirm, don't assume it's
     `FILCSV` by analogy to `FILRAW`/`FILH2O`), and additionally forces
     the internal equivalent of **`LCSV=11`** regardless of what the
     control file says, so CSV output is guaranteed written whenever
     `-csv` is used.
   - **CLI overrides control file, unconditionally, for both.** Neither
     needs to be set in the control file at all if given via CLI.
   - **Mechanism — reuse the established `_set_by_cli` pattern from
     goal #2**, this time putting the forward-looking hook to actual
     use: parse CLI values into staging variables before `MYCONT()`
     (as already established), set companion flags (`filraw_set_by_cli`,
     `filh2o_set_by_cli`, `csv_set_by_cli`), then **actively reassign
     the real control variables *after* `MYCONT()` returns**, for
     whichever were CLI-set — don't just set-before-and-hope the
     control file's `NAMELIST` read leaves them alone.
   - **If no met/raw input is available from either CLI or control
     file, exit with an error.** Find whatever check (if any) already
     exists for a missing `FILRAW` and extend/reuse it rather than
     inventing new error logic.
   - **Open items to confirm before implementing — Plan mode, same
     rigor as goals #1/#2, not yet verified:**
     1. **CONFIRMED — reapply-after-`MYCONT()` is real, necessary
        insurance, not just caution.** `MYCONT()`'s entire body
        (807-1282) never assigns `FILRAW`/`FILH2O`/`FILCSV` — the only
        occurrences of these names are the existing validation checks
        (comparing against blank). So a control file that sets its own
        `filraw=`/`filh2o=`/`filcsv=` would otherwise silently override
        whatever the CLI flags set, unless the post-`MYCONT()` reapply
        step actively overwrites them again. Confirmed directly, not
        inferred from `NAMELIST` semantics alone.
     2. **RESOLVED — real, severe, silent-failure risk confirmed;
        decision reversed to add a guard.** `LCSV` is a Fortran *unit
        number* selector (forced to `11`); collides with `LCOORD`/
        `LTABLE`/`lcoraw`/`LPRINT` if the control file independently
        sets any of those to `11` too. **Traced mechanism, confirmed
        via the F77 standard's `OPEN` semantics (a same-unit re-`OPEN`
        silently closes the previous file — not caught by `ERR=`) and
        this program's exact execution order:** `LCSV`'s `OPEN`
        (~283/285) runs once, early, before the per-voxel loop.
        `LCOORD`'s `OPEN` (~2054, inside `open_output()`, called at
        ~375 — the first statement of voxel 1's loop body) runs
        *after* `LCSV`'s. Since it opens unit 11 second, **`LCOORD`
        silently steals the unit from `LCSV`, not the reverse.** Every
        `WRITE(LCSV,...)` in the file (only at ~10901/10904, inside the
        concentration-table output reached via `FINOUT()` at ~452 —
        after `open_output()`) then lands in `FILCOO` instead. Net
        effect: **`FILCSV` ends up completely empty (0 bytes, not even
        voxel 1's row) for the entire run, while `FILCOO` silently
        accumulates every voxel's CSV-formatted row spliced into its
        own coord-file content.** No `ERRMES`, no console warning, no
        nonzero exit — the run reports success exactly as if everything
        worked; the only traces are passive (an unexpectedly-empty CSV
        file, or foreign-looking lines in the coord file, discoverable
        only if someone happens to look).
        **Decision reversed: add a cheap pre-flight guard, do not ship
        "no guard."** This is exactly the silent-wrong-output failure
        class this project has treated as unacceptable everywhere else
        (the `BRUKER`/`SEQACQ` and `ivoxel`/`COMMON` incidents both
        involved silent, undetected wrongness — this would be worse:
        total, silent data loss with zero signal). **Guard design:**
        right where the `-csv` override-reapply logic forces
        `LCSV=11` (after `MYCONT()` returns, so `LCOORD`/`LTABLE`/
        `lcoraw`/`LPRINT`'s control-file values are already known, and
        before `LCSV`'s own `OPEN` at ~283/285 executes), check whether
        any of those four already equals `11`; if so, fail loudly using
        the same `WRITE`+`FLUSH(6)`+`STOP 1` pattern already proven in
        goal #2, naming which variable collided (e.g. *"Error: -csv
        requires Fortran unit 11, but the control file already assigns
        unit 11 to LCOORD. Change one of these settings."*). Small
        addition (~5-10 lines), reusing an already-verified error
        mechanism, closing off a failure mode that would otherwise be
        essentially undiscoverable. Test plan updated to cover one case
        per potentially-colliding variable, not just review the guard.
     3. **CONFIRMED — `FILCSV`.** `source/lcmodel.inc:43`
        (`FILcsv*(MCHFIL)`), opened at `source/LCModel.f:283/285` —
        matches the `FILRAW`/`FILH2O` naming pattern, but confirmed
        rather than assumed.
     4. **CONFIRMED — existing check at `source/LCModel.f:1270-1271`;
        `ERRMES` unsuitable, for a more precise reason than originally
        hypothesized.** `LPRINT` is a `NAMELIST /LCMODL/` member
        (`source/nml_lcmodl.inc:17`), so its *value* is final by this
        check's location — but that's not what actually matters here:
        `LPRINT`'s **file connection** (the `OPEN`, inside
        `open_output()` at `source/LCModel.f:2033`, called at line 375)
        doesn't exist yet at this earlier point (1270-1271 fires before
        375). `ERRMES`'s normal-case print (`source/LCModel.f:2526`)
        has no unconditional fallback for an unconnected unit. So the
        original framing ("`LPRINT` isn't set yet," from goal #2) and
        this one ("`LPRINT`'s value is set, but its file isn't open
        yet") are different reasons landing on the same conclusion —
        worth keeping the distinction precise rather than reusing the
        old reasoning verbatim. **Same fix: `WRITE`+`FLUSH(6)`+
        `STOP 1`.**

     **Design status: fully vetted, diffs reviewed, ready to apply.**
     All four original items resolved, plus the `LCSV` collision risk
     found, escalated, and reversed to add a guard. Diff review caught
     and fixed three further issues before anything touched disk:
     - **Missing truncation checks** on the `-i`/`-csv` value fetches
       (`GET_COMMAND_ARGUMENT` silently truncates an overlong path
       unless `STATUS=`/`LENGTH=` are checked) — fixed, now errors
       clearly with the actual offending length reported.
     - **Collision-guard list expanded from 4 to 8 variables**, via a
       ground-truth sweep (grep every `open (`/`OPEN (` in the file,
       case-insensitive — a case-*sensitive* first attempt would have
       repeated the exact mistake that missed `average()` as a `MYDATA`
       caller earlier in this project) rather than filtering
       plausible-looking namelist names: `LBASIS`, `LCOORD`, `LCORAW`,
       `LH2O`, `LPRINT`, `LPS`, `LRAW`, `LTABLE`. Two candidates
       correctly excluded with reasoning, not by omission:
       `LCONTR_SCRATCH` (confirmed absent from `NAMELIST /LCMODL/` —
       no control file can ever set it) and `LSCRATCH` (confirmed a
       compile-time `PARAMETER`, not even a `COMMON` variable).
     - **A real, would-have-been-a-compile-error bug caught in this same
       review**: the new `/BLCHAR/` continuation line had a blank
       column 6 instead of a continuation character — gfortran would
       have parsed it as a new statement, not a `COMMON`-list
       continuation. Fixed (continuation character `1` at column 6,
       verified by direct character-position measurement).
     Placement safety (the reapply/guard block always executes after
     `MYCONT()` before anything downstream) reconfirmed from multiple
     angles: no computed `GO TO`/`ASSIGN` anywhere in the file, no
     labels in the relevant line range, `PROGRAM LCMODL` is the file's
     only program unit with no re-entry mechanism.

     **GOAL #3 STATUS: DONE.** Diffs applied (+15 lines `lcmodel.inc`,
     +149 lines `LCModel.f`), clean rebuild from scratch, and verified:
     - All three existing regression baselines pass byte-for-byte with
       no flags — default path untouched.
     - Argument-validation errors (`-i` with one path, `-csv` with no
       path, unrecognized flag) all produce correct nonzero exits.
     - Goal #2's `-no-prior`/`-first-prior` still parse correctly
       through the restructured `DO`→`GOTO` loop.
     - **`-i`'s override genuinely verified, not just assumed**: run
       against a control file with deliberately nonexistent
       `filraw=`/`filh2o=` paths, using `-i` pointing to the real files
       — succeeded, proving the reapply actually overrides rather than
       coincidentally matching.
     - **The collision guard verified to actually fire**: control file
       set `LCOORD=11`, ran with `-csv` — got the exact expected error
       naming `LCOORD`, nonzero exit. Closes the silent-data-loss
       failure mode this whole sub-investigation was about.
     - The missing-input check verified to actually halt (the harmless
       pre-existing `MYCONT`-internal check fires first, as expected
       since `EXITPS`'s `STOP` is disabled, then the new post-reapply
       check correctly stops with a clear message).
     - `-csv` happy path: produced a real, populated CSV with correct
       concentration data.
     - **CLI-driven output confirmed byte-for-byte identical to trusted
       references, not just "looks correct."** `-csv` alone (single-
       voxel, control-file `filraw=` left in place since that dataset
       has no water-reference file to pair with `-i`) matches
       `test-reference-out.csv` exactly. `-i`+`-csv` combined
       (`multi-voxel-10` and `multi-voxel`, both with real `.ref` files,
       control file's `filraw=`/`filh2o=`/`lcsv=`/`filcsv=` lines
       stripped entirely) both match `test-reference-multi-voxel.csv`
       exactly. Diff exit code 0, zero output, in all three cases. This
       confirms the CLI-override mechanism drives the identical
       numerical analysis as the control-file-only path — not merely
       "doesn't crash."
     - **Known, deliberately accepted minor UX rough edge, not a bug**:
       the missing-input error case prints a harmless, uninformative
       pre-existing line (`"Error in voxel but continuing"`, from the
       untouched `MYCONT`-internal check) before the new, clear,
       actionable error message. Not misleading (both agree something's
       wrong; the second is specific), just noise ahead of signal.
       Deliberately left as-is: fixing it would mean touching the
       pre-existing, shared `ERRMES` mechanism that was intentionally
       left alone per this goal's scope, for a cosmetic gain that isn't
       worth the added risk. Revisit only if this actually confuses
       real users in practice.
     Not yet committed as of this writing — commit with the same
     discipline as goals #1/#2 (`git status` before, `git show --stat`
     after).

OpenMP parallelization is a possible future goal but is explicitly
**out of scope for now** — don't introduce `!$OMP` directives, threading,
or thread-safety refactors unless asked. Keeping this out of scope
simplifies #1–#3 considerably (no need to reason about concurrent access
while restructuring the I/O and prior logic).

- Constraints: preserve the existing default behavior exactly
  (regression tests must still match `test-reference-out.csv` and the
  multi-voxel baselines byte-for-byte when no new CLI flag is given —
  see the corrected single-voxel criterion in the Testing section
  below; do not use `out_ref_build.ps` as the pass/fail check);
  new flags must be additive, off by default (`UsePrior=1` unless a
  flag says otherwise).
- Definition of done per item: regression protocol below passes for the
  default path (no flags), plus a new reference output is captured and
  manually reviewed for each new opt-in mode (`-no-prior`, `-prior-file`
  — not `-first-prior`, which is being removed) since none currently
  exists (see "New reference outputs" below).

## Goal #2b: external fixed-prior mechanism (replaces `-first-prior`)

**What "priors" actually are, precisely** (worth stating plainly, since
this wasn't fully nailed down before designing `-first-prior`): a small,
fixed-size set of **phase-correction reference values**, not
concentration priors:
- `DEGPPM` — an expected frequency shift, in ppm.
- `DEGZER` — an expected zero-order phase.
- `SDDEGP`/`SDDEGZ` — standard deviations (uncertainty) around those two
  expected values.
- Possibly `DGPPMN`/`DGPPMX`/`SHIFMN`/`SHIFMX` — range bounds; **not yet
  confirmed whether these are part of the same externalizable set or a
  separate mechanism** — confirm before implementing.

Used inside the fit (`PHASTA`, via `EXDEGZ`/`EXDEGP`) as a Bayesian-style
anchor: the fit is penalized if its own phase/frequency estimate strays
too far from these expected values, weighted by the uncertainty
(`ABS(DEGPPM-EXDEGP) .GT. 4.*SDDEGP`-style check, confirmed during goal
#1's `restore_settings()`/`PHASTA` trace).

**Why the default mode's "voxel order" ever mattered, and why
`-first-prior` was flawed:** the default mode isn't "use voxel 1's
prior" — it's an **adaptive/rolling empirical-Bayes scheme**.
`update_priors()` accumulates running sums/sums-of-squares (subroutine-
local `SAVE` state) across every voxel processed *so far in this run*,
continuously refining the mean/SD estimate. `-first-prior` froze that
rolling process after its first update — which is exactly why it
inherited the ring-scan-order ambiguity (see goal #2's entry above):
"first" only meant anything because the *rolling* mechanism cares about
processing order, and under ring-scan order, the first voxel actually
*analyzed* is the scan-center voxel, not `ivoxel=1`. Separately,
freezing an in-map voxel's prior meant that voxel's row had to be
identified and excluded from the CSV output — an ongoing administrative
burden with no clean fix that didn't touch scan order (a bigger,
riskier change — see the ring-order/default-behavior discussion this
session had, which was explicitly not undertaken because it would
silently change default rolling-prior numerical output for every
existing user, a serious reproducibility concern for scientific
software).

**New design — decided, not yet implemented:**
- **Two-step external workflow, chosen deliberately over any in-map
  voxel selection:**
  1. Run LCModel once on a reference spectrum, single-voxel style,
     letting the default (rolling) mechanism compute
     `DEGPPM`/`DEGZER`/`SDDEGP`/`SDDEGZ` normally, as it already does.
  2. Extract those final values from that run.
  3. Feed them into the real map-analysis run via
     **`-prior-file <file>`** — a new CLI flag reusing the freed
     `UsePrior=4` slot (`USEPRIOR_FILE`) — which sets these values
     once at startup and holds them **fixed for the entire run**: never
     touched by `update_priors()`, never touched by `restore_settings()`'s
     dynamic reset-to-zero-check logic.
- **Two distinct flags, not one flag with dual read/write behavior** —
  a single flag inferring direction from context (e.g. "does the file
  already exist?") was considered and rejected as fragile/ambiguous:
  - **`-prior-file <file>`** — read/apply direction. Sets
    `DEGPPM`/`DEGZER`/etc. fixed for the whole run (`UsePrior=4`), as
    described above.
  - **`-save-prior <file>`** — new, write direction. After the run
    computes its prior normally (default rolling mechanism, unaffected
    otherwise), dump the final values to the named file, in the same
    format `-prior-file` expects to read. Used in the *first* step of
    the two-step workflow (running on the reference spectrum); not
    mutually exclusive with anything else that flag's own run needs,
    but see the two error conditions below.
  - Nothing prevents combining both in one run in principle (e.g. "use
    a fixed prior, but also re-log it for reproducibility") — not
    restricted, though not the primary use case either.
- **`-save-prior`'s mean/SD derivation — FINALIZED, via a dedicated
  independent accumulator, not by touching `update_priors()`.** This
  went through real back-and-forth (see the full reasoning trail
  below if it needs re-deriving) — the final design:
  - **Confirmed fact about `update_priors()` itself (worth recording
    independent of this feature — affects understanding of *default*
    behavior generally):** accumulation (`sum`/`sum2`/`nsamples`) is
    **unconditional**, on every call. But *exposure* — writing the
    computed mean/SD into `DEGPPM`/`DEGZER`/`SDDEGP`/`SDDEGZ`/etc. — is
    gated behind `if (nsamples .lt. max0(2, mnsamp)) return`
    (`mnsamp` defaults to `9`). Below that threshold, these `COMMON`
    variables are untouched no-ops, silently retaining whatever they
    held before (control-file default, or nothing). Once the gate
    opens, **every subsequent call recomputes and overwrites them from
    the full cumulative history** (an expanding-window mean/SD, not
    compute-once) — this is a genuinely continuously-refining
    empirical-Bayes scheme for the entire remainder of a default-mode
    run, not just an initial warm-up. Practical implication: a small
    multi-voxel dataset (e.g. `multi-voxel-10`, 10 voxels) barely
    exceeds the default `mnsamp=9` threshold — most of a small run's
    default-mode prior state is effectively frozen at its initial
    value.
  - This subroutine-local `SAVE` accumulation is invisible outside
    `update_priors()` — there's no way for `-save-prior` to read
    partial/below-threshold accumulated sums from it, even though the
    *mean* formula itself has no mathematical degeneracy below
    threshold (only the `/(n-1)` SD term does). Two ways to close this
    gap were considered: modifying `update_priors()` to also expose
    partial means below its own gate (rejected — touches shared
    default-mode code for every user, disproportionate risk), vs. a
    **separate, independent accumulator** used only when `-save-prior`
    is active (**chosen**) — zero changes to `update_priors()`, zero
    risk to default behavior.
  - **Final logic**, mirroring `update_priors()`'s math in isolated
    state, reusing `mnsamp`'s meaning for consistency:
    - **`DEGZER`/`DEGPPM`**: always computed from whatever's
      accumulated so far by the independent accumulator, even from a
      single voxel — no gate, since the mean has no n=1 degeneracy.
    - **`SDDEGP`/`SDDEGZ`**: computed from the accumulator's own sums
      **only if** `nsamples .ge. max0(2, mnsamp)` (reusing `mnsamp`,
      not an independently-configurable threshold — decided for
      consistency of meaning across the program); **otherwise, fall
      back to the control file's static default `SDDEGP`/`SDDEGZ`**
      values — matching the paper's own convention (Provencher 1993)
      of using a fixed, non-data-derived default when there isn't
      enough data for a real variance estimate.
  - Works uniformly for single-voxel *and* multi-voxel reference runs
    — supersedes the earlier, now-obsolete "restrict to single-voxel
    only" idea. A single voxel just means the accumulator has
    `nsamples=1`: mean from that one value, SD from the control-file
    default (below `mnsamp`'s floor of `max0(2,...)=2`, so always
    true for n=1 regardless of `mnsamp`'s setting).
  - **Not yet designed**: exactly where in the per-voxel flow this new
    accumulator hooks in (presumably alongside wherever `update_priors()`
    itself gets called, but as a fully separate call/state — not a
    modification to that subroutine), and confirming `PHITOT`'s role
    as the correct per-voxel input value in a multi-voxel context too,
    not just the single-voxel case already traced.
- **One error condition still applies:** **`-save-prior` combined with
  `-no-prior` → error.** Under `-no-prior`, nothing resembling a fit
  result worth saving is being computed in the intended sense — error
  out rather than silently saving something meaningless. (The earlier
  "restrict to single-voxel, error on multi-voxel" condition is
  superseded by the hybrid design above and no longer applies.)
- **`-prior-file`/`-save-prior`'s shared format: a single small text
  file** (decided over individual per-value CLI flags) — simple
  `key=value` lines expected (e.g. `degppm=0.05`, `degzer=12.3`,
  `sddegp=...`, `sddegz=...`); exact format/parsing approach, and the
  final variable list (see open item 1 below), not yet designed in
  detail.
- **Mutual exclusivity carried forward**: `-no-prior` and `-prior-file`
  together should still error, same pattern as the old `-no-prior`/
  `-first-prior` rule.

**Open investigation items — Plan mode, not yet resolved:**
1. Confirm the exact, complete list of variables that constitute "the
   prior" — is it just `DEGPPM`/`DEGZER`/`SDDEGP`/`SDDEGZ`, or do
   `DGPPMN`/`DGPPMX`/`SHIFMN`/`SHIFMX` also need externalizing? Trace
   `update_priors()`/`restore_settings()` precisely, don't assume from
   the goal #1-era trace summary alone. This list determines both
   `-prior-file`'s required fields and what `-save-prior` needs to
   write. **Partially informed by the finalized mean/SD design above**
   — `DGPPMN`/`DGPPMX`/`SHIFMN`/`SHIFMX` weren't addressed by that
   design and still need their own sourcing decided (control-file
   passthrough, like the SD fallback? Or something else?).
2. **Is a new "dump priors" feature (`-save-prior`'s core mechanism)
   genuinely new code, or can it reuse/adapt something that already
   prints these values somewhere** (a coord file, the `.ps`, stdout)?
   Check the actual existing output before assuming either way — if
   nothing already surfaces these values, `-save-prior`'s write logic
   is new code, not a small wrapper around existing output. **Note:**
   given the new independent-accumulator design, the values to write
   won't be sitting in `update_priors()`'s output `COMMON` variables at
   all in the below-threshold case — they'll be in the new
   accumulator's own state — so this write logic needs to read from
   wherever that new state lives, not just from `DEGPPM`/`DEGZER`.
3. Design the `-prior-file`/`-save-prior` parsing (format, required vs.
   optional fields, error behavior on malformed/missing fields) and the
   `restore_settings()` branch for `UsePrior=4` (analogous in spirit to
   the old `first_prior_done` wrapper, but simpler — no need to capture
   anything mid-run, just read once at start and hold fixed every
   voxel).
4. Confirm removing `-first-prior`'s code doesn't leave any dangling
   references (call-site sweep, same rigor as before) — `UsePrior`
   comparisons elsewhere in the file that check `.eq. USEPRIOR_FIRST`
   need finding and removing/updating, not just the CLI branch and the
   `restore_settings()` wrapper.
5. **Design the new independent accumulator itself**: where it hooks
   into the per-voxel flow (alongside wherever `update_priors()` is
   normally called, but as fully separate state — not a modification to
   that subroutine), what per-voxel value(s) it accumulates (`PHITOT`
   confirmed correct for the single-voxel case — confirm it's also the
   correct per-voxel input in a multi-voxel context, not just assumed
   to carry over), and where its final state gets read at end-of-run
   for `-save-prior`'s write.

**GOAL #2b STATUS: DONE.** Diffs applied, clean build (17 pre-existing
unrelated warnings only), `-first-prior` fully removed (case-insensitive
sweep confirms zero remaining references except the intentionally-kept
`PARAMETER`), all three existing regression baselines still pass
byte-for-byte with no new flags. Full closure:

- **A real numerical bug caught and fixed during diff review, before
  applying anything**: the first drafted end-of-run write for
  `-save-prior` mirrored `update_priors()`'s raw SD formulas but
  omitted its final clamps (`sddegp = amin1(sddegp_input,
  amax1(sd(2), sddegp_min))`, `sddegz = amax1(sd(1), sddegz_min)`,
  with `sddegp_min=1.0`/`sddegz_min=3.0`). Caught by comparing the new
  code line-by-line against `update_priors()`'s actual source rather
  than accepting "mirrors the formula" at face value — fixed by
  duplicating the two floor/ceiling literals (same accepted
  manual-sync tradeoff as the `msamples=4096` cap and `MCHFMT`/`MCHID`
  in earlier goals).
- **Both mutual-exclusivity checks verified order-independent**, traced
  through actual variable states for both orderings, not asserted from
  symmetry of form: `-no-prior`/`-prior-file` and `-save-prior`/
  `-no-prior` each fire correctly regardless of which flag appears
  first on the command line.
- **Unit-14 collision guard verified to actually fire** (`LCOORD=14`
  test case), same pattern as goal #3's `LCSV=11` guard.
- **`-save-prior` cross-verified against an independent computation
  path, not just "ran successfully"**: single-voxel test's written
  `degzer`/`degppm` matched the `.ps` output's own "Ph:" annotation
  exactly. Multi-voxel-10 test (10 voxels, above the `mnsamp=9`
  threshold) confirmed the threshold-gated real-SD path genuinely
  fires (data-derived `sddegz`, distinct from the below-threshold
  fallback case).
- **A dataset-specific zero explained and independently confirmed, not
  just observed and accepted**: `degppm=0`/`sddegp=0` in the
  multi-voxel-10 output are zero for two genuinely different reasons —
  `sddegp` because the clamp correctly mirrors `update_priors()`'s own
  ceiling against this control file's `SDDEGP=0` setting; `degppm`
  because `PHITOT(2)` is a hard structural zero for every voxel of
  this specific dataset (confirmed by instrumenting unmodified
  default-mode execution — `EXDEGP=0` forces `PHASTA`'s anchor
  uncertainty budget to `4*SDDEGP=0`, so `DEGPPM` never moves off its
  starting point, for any voxel, independent of anything added by this
  goal). Not the same explanation reused for two different variables.
- **`-prior-file`'s mechanism verified end-to-end, not just
  parse-tested**: ran the same multi-voxel-10 dataset with two
  meaningfully different prior files; fitted `DEGPPM`/`DEGZER`
  genuinely shifted toward whichever prior was supplied — a
  substantial, consistent, directionally-correct response, not
  identical output regardless of input. This is the test that would
  have caught a bug where the file parses correctly but the fixed
  values never actually reach `PHASTA`.
- **Debug instrumentation added mid-investigation was fully removed and
  re-verified** — confirmed via grep (zero remaining debug identifiers)
  and a full rebuild + re-run of all three regression baselines after
  removal, not just assumed clean.
- **Not yet tested, deferred explicitly (not silently skipped), same
  pattern as goal #1's checkpoint/resume**: the `msamples=4096` cap and
  its one-shot warning. No >4096-voxel fixture exists in this repo; a
  synthetic test (temporarily lowering the duplicated cap literal,
  running the 100-voxel `multi-voxel` dataset, confirming the warning
  fires exactly once, then reverting the literal) is cheap and
  available if this ever needs closing out — lower priority than
  everything above since it doesn't affect either flag's core
  correctness, just the cap's own warning path.
- **Not yet committed as of this writing.**

## Goal #4: ABANDONED (was: `-nvoxels <N>` CLI flag)

Was scoping a new CLI flag to force `NDCOLS`/`NDROWS`/`NDSLIC`/`ICOLST`/
`ICOLEN`/`IROWST`/`IROWEN`/`ISLICE` for a single-column-of-voxels
workflow (per Provencher's manual, Sec 5.3.2 — `docs/manual.pdf`). No
code was ever applied — this was Plan-mode investigation only, so
there's nothing to revert. **Abandoned before implementation.**

**One thing surfaced during this investigation is real and independent
of this abandoned goal — worth keeping on record:**

**The existing `FATAL` validation on these 8 grid variables is
toothless.** If any of them are `0` or otherwise invalid, the program
prints a message claiming to be a "FATAL ERROR" but doesn't actually
halt — same disabled-abort pattern found earlier in this project
(`EXITPS`'s `STOP` is commented out project-wide, deliberately, so one
bad voxel doesn't kill a multi-voxel batch run). This means a control
file that sets, e.g., `NDCOLS=0` today silently continues running in
some degenerate/undefined state rather than actually stopping, for
*any* control file — not just ones that would have used `-nvoxels`.
**Deliberately left unfixed here** — fixing it means touching a shared
error-handling mechanism (`EXITPS`) that other parts of the program,
including everything built in goals #1-3, currently depend on behaving
the way it does; that's a bigger, riskier change than anything scoped
so far and deserves its own dedicated investigation, not a side effect
of an abandoned goal's diff. **If this ever needs fixing**, start by
mapping everything else that relies on `EXITPS` not actually stopping,
before touching it.

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

- `make test_lcm/out.ps` (or the corresponding `.csv` target, once
  wired in — see correction below) — original single-voxel case.
  **⚠️ CORRECTED CRITERION — read before trusting any `out.ps` diff:**
  `out_ref_build.ps` (upstream's checked-in reference) was discovered to
  differ from this fork's own legitimate, pre-existing output for
  reasons unrelated to any regression — it predates this fork's own
  divergence from upstream, so diffing against it was never actually
  the right test. **The correct baseline is this fork's own output at
  commit `621487d`** (confirmed start-of-fork baseline, immediately
  before `MODULE RAWCACHE`), captured as `test-reference-out.csv` via
  the same control-file CSV-output option the multi-voxel tests use —
  not a `.ps` diff, which also carries an unavoidable build-date
  timestamp that makes a naive `diff` always show a difference even
  with zero real change. The `.ps` output can still be generated as a
  secondary human-readable sanity check, but it is **not the pass/fail
  criterion** — the `.csv` diff against `test-reference-out.csv` is.
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
3. Diff each resulting `.csv` output against its actual reference
   (single-voxel `.csv` vs `test-reference-out.csv`; each
   `multi-voxel.csv` vs its own `test-reference-multi-voxel.csv`).
   **Do not rely on a raw `.ps` diff as the pass/fail signal** — see the
   corrected criterion above.
   - Any difference in fitted concentrations, CRLBs, SNR, or plotted
     spectra is a regression unless the change was *intended* to alter
     the algorithm — flag it explicitly, don't silently accept it.
     the algorithm — flag it explicitly, don't silently accept it.
4. If the change is expected to alter numerical output (e.g. an
   intentional algorithm change), say so up front and show a before/after
   comparison rather than just a passing/failing diff.

If more test datasets exist or are added, extend this protocol to cover
them — don't rely on the three bundled tests as permanently sufficient
coverage once the project grows.

### New reference outputs for this project

Reference outputs for the single-voxel (`test-reference-out.csv`,
generated from commit `621487d`, this fork's own confirmed baseline —
**not** upstream's `out_ref_build.ps`, which was found to diverge from
this fork for reasons unrelated to any regression), `multi-voxel`, and
`multi-voxel-10` cases (each directory's `test-reference-multi-voxel.csv`)
all exist — use those as the baseline rather than generating new ones
from scratch.

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

**✅ INCIDENT RESOLVED — root cause found and fixed, re-verified clean.
Read this before trusting "GOAL #1 STATUS: DONE" below; it explains why
that status was briefly, genuinely wrong, and why it's now trustworthy
again.**

**What happened:** commit `9f7b9e8` (the goal #1 commit) turned out to
be an **internally inconsistent snapshot** — `source/LCModel.f` at that
commit already contained `MODULE RAWCACHE`/`check_zero_voxels`/`MYDATA`
code that assumes `ivoxel` is a shared `/BLINT/` `COMMON` member, but
the committed `source/lcmodel.inc` at that same commit did **not**
actually contain `ivoxel` in `/BLINT/`'s member list. Silent Fortran
fallback: with no matching `COMMON` declaration, `gfortran` gave
`PROGRAM LCMODL` and `MYDATA` **two independent implicit-local**
`ivoxel` variables instead of one shared one — no compile error, no
warning. The main loop's real `ivoxel` (correctly `1`, for the
single-voxel test) and `MYDATA`'s completely disconnected local
`ivoxel` (uninitialized, containing stack garbage — confirmed via gdb:
`1431020856`) were different memory locations entirely. `MYDATA` then
indexed `datat_cache(1:NUNFIL, ivoxel)` — a 1-column array — with that
garbage value, reading essentially random memory, which became `DATAT`,
propagated through the entire fit, and surfaced as a deterministic
`NaN` (same garbage value every run, since it was `COMMON`/BSS-adjacent
memory, not stack noise) that crashed `CHREAL`'s axis-label formatting.

**How it was found:** bisected between `621487d` (clean) and `9f7b9e8`
(crashes) to confirm the bug was introduced by goal #1's changes; built
a debug binary (`-g -O0 -fbacktrace -ffpe-trap=invalid,zero,overflow`)
to get an exact crash location instead of tracing by hand; used gdb to
compare `&ivoxel`'s actual memory address at the ring loop's call site
vs. inside `MYDATA` — **different stack addresses** proved they weren't
the same variable at all, which led directly to finding the missing
`/BLINT/` member.

**Why "regression-clean, all three tests pass, 9343-voxel validated"
had been reported earlier despite this bug existing:** unresolved with
certainty, but the most likely explanation, pieced together from the
session: mid-session, an unrelated `git stash` (run during the
`raw_cache_ok`-forcing bisection attempt for a *different* diagnostic
purpose) was never popped, so several turns of "confirmed" debugging
were actually run against **bare committed HEAD**, not the working
tree's actual (correct) code — the stash, once recovered, was found to
already contain the correct `ivoxel` `COMMON` member. The earlier
"all three pass" report likely reflected a state that was never
actually captured cleanly as its own commit before other edits and a
stash operation overlapped it. **Lesson captured for this project:**
after any `git stash`, verify `git status` shows a clean/expected state
before trusting further test results — an un-popped stash makes "what
you're testing" and "what you think you're testing" silently diverge,
exactly as happened here.

**Fix applied and independently re-verified, not just restored:**
`ivoxel` correctly added to `/BLINT/` in `lcmodel.inc`. Verified via
gdb: `&ivoxel` at both the ring loop and inside `MYDATA` now resolves to
the **identical address**, reported by gdb as `<blint_+40836>` — an
offset into the actual `/BLINT/` linker symbol, direct confirmation of
real shared `COMMON` storage, not two disconnected locals. Value `= 1`
at both breakpoints for the single-voxel test, matching expectation.
Full clean rebuild performed (stale `.mod`/`.o`/binary deleted first).
**All three regression targets re-run from scratch and confirmed
passing** (single-voxel: `NaN`/crash gone, matched `out_ref_build.ps` at
the time — **note:** `out_ref_build.ps` was itself later found to be
the wrong baseline for this fork, see the corrected criterion in
"Testing / regression protocol" above; the crash/`NaN` fix itself is
unaffected by this — a crashing run and a clean run are distinguishable
regardless of which reference file you diff against — but the "matches
its reference exactly" framing at this point in the log used a
reference that was subsequently replaced);
both multi-voxel CSVs match their `test-reference-multi-voxel.csv`).

**GOAL #1 STATUS: DONE — for real this time, independently re-verified
after the above incident, not just restored from the original claim.**
Implemented,
correct, and regression-clean on
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
- **✅ Real-world large-dataset validation — confirmed complete, both
  performance and correctness, at production scale.** A real
  9343-voxel dataset: new (cached) version ~1.6 hours; old (uncached)
  version projected 8–10 hours — roughly a 5–6× speedup, exactly the
  scaling behavior predicted (the benefit grows with grid size, which
  is why it was too small to see on the 100-voxel fixture). **Output
  confirmed identical between old and new versions on the complete
  9343-voxel dataset** — initially only the first ~3600 voxels had
  been checked (the old run hadn't finished at comparison time); this
  was later completed and reconfirmed on the full dataset. Both the
  speedup and the output match were **independently re-verified a
  second time**, after all of this goal's later fixes (the `ivoxel`/
  `COMMON` incident, goal #2's changes) — an initial re-check appeared
  to show the speedup had disappeared, but this turned out to be a
  false alarm caused by accidentally running a stale/old binary, not a
  real regression; re-run with the correct current binary confirmed
  both the ~5–6× speedup and exact output match hold on current HEAD,
  not just on some earlier commit.
  This closes the original motivating complaint ("very slow" multi-voxel
  runs) — the fix has both correctness and performance validated at
  a scale representative of real use, confirmed twice, not just on small
  regression fixtures or a single early check.
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
