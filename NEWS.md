# MVMRcMLSuSiE 0.1.0.9001 (development)

* Declare MRcML and qualify external calls so the offline workflow runs in a
  fresh session without manually attaching dependencies.
* Keep TwoSampleMR optional for online retrieval, with an informative missing
  dependency error. The separate MVMRcML post-selection analysis remains in
  the original README and is outside the initial offline check.
* Bound the outer SuSiE reweighting loop with `max.outer.iter` (default 200).
  Preserve the existing update equations, inner settings, and convergence
  criterion. Stop explicitly on failure to converge; record diagnostics on
  successful fits.
* Preserve matrix dimensions after invalid-IV exclusion and reject malformed
  numerical inputs before fitting.
* Correct the README's bundled outcome field names from `.mat` to `.vec`.
* Add a small offline runner using bundled data and focused regression tests.
  This is a runtime repair, not a new scientific validation or a full rerun
  of all 249 exposure screens.
* Add `inst/scripts/offline-full.R` to recompute the complete bundled offline
  chain: all 249 screens, all selected step-2 fits, and step 3 using those new
  results. Save per-exposure checkpoints and historical comparison tables;
  stop if the new selection differs from the joint input available in the bundle.
* Record successful small and full offline runs on Apple Silicon macOS with
  R 4.3.2 and susieR 0.12.35. The full run retained the historical 43-exposure
  selection and 10-invalid-instrument set, completing in 2308.244 seconds.
  Its individual warning messages were not retained and remain to be reviewed.
* Add optional `--workers=N` parallel step-1 fits with ordered collection and
  the existing per-fit random seed. Keep the default serial and steps 2-3
  unchanged. Capture warnings by exposure/stage and save a separate summary
  without truncating them to the first 50 conditions.
