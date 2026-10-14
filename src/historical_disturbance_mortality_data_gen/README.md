# historical_disturbance_mortality_data_gen

Builds the historical disturbance mortality rates datacube used as a calibration target.
Pipeline entry point: [00_run_all.jl](00_run_all.jl) (steps `02`-`06`; see note on `01` below).

[base.jl](base.jl) is included at the start of every numbered step and, besides the usual
config/path setup, guards three prerequisite pulls - each only regenerated if its target file
doesn't already exist (delete the file to force a fresh pull):

- [base/generate_manta_tow_cache.jl](base/generate_manta_tow_cache.jl) writes
  `data/manta_tow_raw_cache.csv` (guarded on that file's own existence): a local cache of
  `ReefMonitoring.get_manta_tow`'s raw per-reef, per-year output, read by both
  [03_a_create_cover_manta_csv.jl](03_a_create_cover_manta_csv.jl) and
  `generate_clean_growth_intervals.jl` below instead of each hitting the API independently for
  the same underlying data on every run.
- [base/generate_calibration_split.jl](base/generate_calibration_split.jl) writes
  `data/calibration_split.csv` (guarded on that file's own existence) from a bare ADRIA
  domain, sidestepping the circular dependency on `historical_disturbance_mortality_rates.nc`
  described in the file's header comment.
- [base/generate_clean_growth_intervals.jl](base/generate_clean_growth_intervals.jl) writes
  `data/clean_growth_intervals.csv` (guarded on `data/background_growth_rates.csv`'s
  existence, not its own - see the file's header comment) via a ReefMonitoring disturbance-
  history pull plus the manta tow cache above. Also applies an automatic suspicious-cover-drop
  exclusion filter at generation time - see Decisions below.

**[01_create_disturbance_years_csv.jl](01_create_disturbance_years_csv.jl) is deliberately not
in the routine run.** It only writes `template_disturbance_years_{manta,transect}.csv` -
throwaway drafts for manual curation into `disturbance_years_manta.csv`/
`disturbance_years_transect.csv` (see its header comment) - nothing later in the pipeline reads
its output, and it never overwrites the already-curated real files, only the templates. Run it
standalone only when new disturbances have appeared in ReefMonitoring since those real files
were last curated: regenerate the templates, diff against what's already curated, and manually
pick `year_before`/`year_after` (by eye, from the Reef Monitoring Dashboard) for whatever's
new before merging it in. Not needed just to rerun the pipeline on already-curated data.

## Decisions

- Disturbance survival is estimated as `cover_after / projected_undisturbed_cover`, where
  `projected_undisturbed_cover` is `cover_before` grown forward across the survey gap using a
  background (undisturbed) growth rate - i.e. mortality is measured against what the reef
  would plausibly look like had the disturbance not happened, not against raw pre-disturbance
  cover.
- **Background growth is modelled ONCE, on TOTAL manta tow cover (population-level), not per
  functional group.** This is the final form after a multi-stage redesign:
  - Originally additive (absolute cover-proportion increment, stratified by
    `(functional_group, gap_years)`) - wrong operator for a multiplicative process; a pooled
    absolute increment dominated the correction for low-cover reefs and produced physically
    impossible `survival_rate > 1` for >30% of `Small_Massives`/`Large_Massives` records.
    [diagnostics/A_growth_rate_by_gap_length.jl](diagnostics/A_growth_rate_by_gap_length.jl)
    and
    [diagnostics/B_starting_cover_correlation.jl](diagnostics/B_starting_cover_correlation.jl)
    (extended mid-redesign to report both absolute and relative growth) motivated the switch.
  - Switched to relative growth, applied multiplicatively; found strong, consistent negative
    correlation between `cover_start` and relative growth (real density dependence). A linear
    per-year rate applied over multi-year gaps then blew up (projected growth up to +462% for
    long gaps) - redefined as an ANNUALIZED/compounding rate,
    `(cover_end / cover_start)^(1/gap) - 1`, applied as `cover_before * (1 + rate)^gap`.
  - Still didn't move `survival_rate > 1` for the massives groups (stayed ~30-37%). Root cause:
    growth was estimated and applied independently PER FUNCTIONAL GROUP, on cover already split
    out of the manta tow total via photo-transect composition proportions - each group's
    before/after numbers carried their own independent single-point transect sampling noise,
    and numerator/denominator used two different (before vs after) transect samples not
    guaranteed to resurvey the same patch of reef. No growth-model change could fix that.
  - Fixed by modelling growth ONCE on the trustworthy, reef-wide manta tow total
    ([02_estimate_background_growth_rates.jl](02_estimate_background_growth_rates.jl)), a
    single continuous regression `relative_growth ~ intercept + slope_cover_start *
    cover_start + slope_gap * gap_years` with a cluster-bootstrap-by-reef CI. Transect
    composition is used ONLY as a weighting to split the corrected total into
    per-functional-group baselines
    ([05_create_survival_rates_csv.jl](05_create_survival_rates_csv.jl)), never as an
    independent basis for its own growth model.
  - [diagnostics/C_growth_by_cover_quintile.jl](diagnostics/C_growth_by_cover_quintile.jl)'s
    finding that growth turns negative in the top cover quintile for `Small_Massives`, and
    [diagnostics/D_gap_cover_confound.jl](diagnostics/D_gap_cover_confound.jl)'s check that
    the gap-length effect isn't a spurious gap/cover confound, both predate the population-level
    switch but motivated treating `cover_start` and `gap_years` as joint linear predictors
    rather than binning.
- **Fixed a fractional-date boundary bug in the clean-interval window check.** Disturbance
  dates from the API are fractional calendar years (e.g. a March cyclone ≈ `2009.18`), but
  were compared against integer survey years with a raw `y1 <= d <= y2` - silently missing
  any disturbance after Jan 1st of the window's end-year. Fixed to compare calendar years via
  `floor` in
  [base/generate_clean_growth_intervals.jl](base/generate_clean_growth_intervals.jl). Removed
  243 of the original 1291 candidate intervals (≈18.8%) once corrected.
- **Suspicious cover-drop exclusion is automatic, not a hand-curated list.** Even with the
  date fix, the ReefMonitoring disturbance flags can have false negatives. A >30% manta tow
  cover drop (checked on raw manta cover, not per functional group) in a supposedly "clean"
  window is excluded at generation time in
  [base/generate_clean_growth_intervals.jl](base/generate_clean_growth_intervals.jl). A
  static exclusion list was tried first and rejected - it would silently drift out of sync as
  source data changes. Threshold pre-registered, not tuned post-hoc.
- **Clean intervals are restricted to `year_start >= 2008`
  (`MIN_YEAR_START`, [02_estimate_background_growth_rates.jl](02_estimate_background_growth_rates.jl)),
  matching the calibration period itself** - older survey-era data may not be comparable.
- **An "exclude every decreasing interval" filter was tried and reverted.** Unlike the
  pre-registered 30% suspicious-drop filter (which targets drops implausible as measurement
  noise), excluding ANY decrease - down to a fraction of a percent - has no biological
  justification: real survey-to-survey noise straddles zero symmetrically, so removing only its
  negative half while keeping the positive half is a one-sided truncation that biases the
  estimated background growth rate upward (visible in retained data as `cover_start` falling
  and `cover_end` rising with gap length). Decision: keep only the 30% threshold, which is
  targeted at genuinely-implausible drops rather than all decreases.
- **The growth-rate fit is WEIGHTED per interval, `1 / (intervals from that reef)`.** The
  double loop in `generate_clean_growth_intervals.jl` emits every pairwise `gap <= 6`
  combination of a reef's surveyed years, not just adjacent-survey pairs, so a reef with a
  long, gap-free survey run generates far more (highly overlapping, non-independent) rows than
  one with the same growth signal but sparser surveys - confirmed empirically: one reef with 11
  real survey years produced 32 intervals, ~6% of the whole unweighted fit's influence, vs a
  median of ~5 intervals/reef. Unweighted OLS let survey density, not growth rate, help drive
  the fit. Fixed in
  [02_estimate_background_growth_rates.jl](02_estimate_background_growth_rates.jl) by weighting
  each row so every reef contributes equal total weight regardless of its interval count
  (weights fixed from the unresampled data, so the cluster bootstrap's resampling multiplicity
  is preserved rather than normalized away).
- **The growth-rate fit uses CALIBRATION reefs only, not calibration + test.**
  `_target_loc_ids` in [base.jl](base.jl) covers every reef (both `USAGE` values, since
  mortality itself must be computed for test reefs too), but using it to build the growth
  intervals meant test reefs were feeding the correction their own mortality is later
  judged against - 146 of 530 growth observations (≈28%) came from test reefs before this
  was caught. Fixed by adding a separate `_calibration_loc_ids` (filtered on
  `USAGE == "calibration"`), used only by
  [base/generate_clean_growth_intervals.jl](base/generate_clean_growth_intervals.jl); every
  other step still uses the full reef set.
- **`survival_rate` is clamped to `[0.05, 1.0]`**, matching
  [06_build_dist_mortality_datacube.jl](06_build_dist_mortality_datacube.jl)'s existing
  `clamp(1 .- survival_rate, 0.0, 0.95)` mortality ceiling exactly. Values above 1 are a real,
  expected residual of applying a population-level growth correction to individual noisy reef
  trajectories - not eliminated, just bounded. Currently ~12.3% of the 830 (reef, disturbance,
  functional group) records land above 1 pre-clamp (this number moved up slightly, in stages,
  as each bias fixed above - the decrease-filter reversion, the interval weighting, and the
  calibration-only restriction - was corrected; read as the estimate becoming progressively
  less biased, not worse).
- **`NaN` survival values (0/0 - a functional group recorded as absent both before and after a
  disturbance) are replaced with `survival_rate = 1.0`**: no coral of that group was present to
  have died, so there's no evidence of mortality to record. Confirmed this is the only source of
  non-finite values - total manta cover (the numerator/denominator's real anchor) is never zero.
- **Composition taxonomy mismatch between the live pull and the frozen dataset, fixed.**
  `generate_clean_growth_intervals.jl` used to build transect composition proportions from a
  live `ReefMonitoring.multiple_location_comparison` pull, which returned a degenerate 50/50
  Tabular/Corymbose Acropora split disagreeing with the frozen `coral_composition.nc` used by
  [03_b_create_cover_transect_csv.jl](03_b_create_cover_transect_csv.jl) (e.g. Agincourt Reef
  No.1 2010: 0.066 vs 0.126 in the frozen data, identical in the live pull). Fixed by reading
  the same frozen `coral_composition.nc` in both places. The population-level growth regression
  itself was never affected (it only reads total manta cover), but the per-group diagnostic
  columns in `clean_growth_intervals.csv` (consumed by diagnostics B/C) were. Bulk-reads the
  per-reef composition slice once rather than one disk read per (year, taxa) pair - the
  scalar-at-a-time access pattern was found to crash the HDF5 driver outright
  (`EXCEPTION_ACCESS_VIOLATION` in `H5F_addr_decode`) after enough reefs.
- **CV contamination bug fixed in [04_fill_missing_transect_disturbances.jl](04_fill_missing_transect_disturbances.jl).**
  The global (all-reef) coefficient-of-variation comparison used to decide "global mean vs
  bioregion mean" for imputing missing transect cover was computed over rows that already
  included the all-zero placeholder rows for the missing records themselves, inflating the
  global CV and systematically biasing the choice toward the bioregion mean. The bioregion-side
  CV already excluded them. Fixed by masking both sides the same way.
- **Manta cover masking bug fixed in
  [03_a_create_cover_manta_csv.jl](03_a_create_cover_manta_csv.jl).** Disturbance records were
  matched to manta cover by `reef_name + cyclone_name` alone; repeated same-named (mostly
  "unknown") disturbances at a reef shared one mask and silently overwrote each other's cover
  data. Fixed to match on `reef_name + year`, with an `@assert` guard against future
  reintroduction.
- **CI usability gating (relative/absolute bootstrap CI width thresholds) is explicitly not
  implemented.** `background_growth_rates.csv` carries `lo`/`hi` columns regardless of width;
  nothing downstream currently checks them before use. Conscious call, not an oversight.
- **Circular calibration_split.csv dependency fixed via a bare ADRIA domain**, not a change to
  the shared `ADRIACoralBloxCalibration.load_domain` wrapper (a wrapper-level fix was drafted and
  explicitly reverted in favor of scoping the fix to
  [base/generate_calibration_split.jl](base/generate_calibration_split.jl) only). Verified
  byte-identical output to a split built from a fully historical-data-attached domain.
- **Bioregion stratification (was Q2) is explicitly not investigated.** The bioregion-variance
  diagnostic that would have answered it was deleted along with the per-reef-imbalance sanity
  check - same conscious-drop treatment as CI usability gating, not an oversight.
- **[diagnostics/E_final_datacube_sanity_dashboard.jl](diagnostics/E_final_datacube_sanity_dashboard.jl)**
  compares the final datacube against real disturbance data: datacube mortality by functional
  group, the intermediate `survival_rate` distribution against the `[0.05,1.0]` clamp bounds,
  and a wiring check confirming the CSV-implied mortality matches the datacube exactly.
