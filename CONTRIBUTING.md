# Contributing

Contributions should preserve the project's fail-closed evidence model.

1. Open a focused issue or pull request.
2. Add a source locator and unit conversion for every physical parameter.
3. Keep different adsorbents and gas systems in separate adapters.
4. Add tests for balances, positivity, topology and failure behavior.
5. Declare calibration data, holdouts and thresholds before inspecting
   holdout results.
6. Do not commit publisher PDFs, restricted datasets, credentials or
   proprietary simulator files.
7. Run `python -m pytest -q` before submission.

A numerical pass must not silently open design, DWSIM, TEA/LCA or training-data
release. Those decisions belong to explicit, reviewable release audits.
