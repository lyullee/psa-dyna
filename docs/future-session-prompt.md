# Future-session bootstrap prompt

The text below can be pasted into a new research or coding session.

> Continue the PSA-Dyna research from <https://github.com/lyullee/psa-dyna>.
> First clone or inspect the repository and read `docs/research-handoff.md`,
> `docs/research-state.json`, `docs/decision-log.md`,
> `docs/known-issues.md`, `docs/validation-playbook.md` and
> `docs/roadmap.md`. Verify the exact Git commit and run the complete tests
> before changing anything. Treat `validation/validation_summary.json` as the
> public metric source. Preserve the default-deny evidence model: a numerical
> pass must not open design, DWSIM, TEA/LCA, optimization or training release.
> Do not mix parameters across papers, gases or adsorbents. Do not treat Wang
> single-bed breakthrough as cyclic validation. Preserve the negative Boer
> apparatus-identifiability result and the closed S-015 design gate. The
> highest-priority next task is to obtain or identify a source-complete
> independent cyclic H2 PSA dataset, preregister calibration/holdout metrics,
> and run a no-refit comparison. Keep copyrighted PDFs, figures, proprietary
> simulator files, credentials and private local paths outside the public
> repository. Document every changed gate in the machine-readable state,
> validation summary, tests and changelog.

Before accepting a task as complete, answer these questions:

1. Which evidence level was reached?
2. Which exact hash-bound files support it?
3. Which release flags remain false?
4. Did any parameter or threshold change after viewing a holdout?
5. Can a fresh environment reproduce the result from public files alone?
6. Is any required source legally unavailable for redistribution?

