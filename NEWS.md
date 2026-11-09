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
* Record the successful small-check run on Apple Silicon macOS with R 4.3.2
  and susieR 0.12.35. This does not establish completion of the full Mac run.
