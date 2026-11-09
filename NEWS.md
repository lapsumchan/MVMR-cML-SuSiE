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
