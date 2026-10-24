# Data, copyright and provenance policy

## Included

- original PSA-Dyna source code;
- original documentation and tests;
- small case definitions that cite their parameter sources;
- derived validation metrics and release decisions;
- hashes that identify the local evidence used during development.

## Excluded

- publisher PDFs and article figures;
- supplementary files without clear redistribution permission;
- proprietary Aspen Plus or DWSIM binaries;
- unpublished author data;
- large internal calculation folders and checkpoints;
- local absolute paths, user identifiers and credentials.

## Reproducing literature comparisons

Obtain the cited source through the publisher, author or official repository.
Respect its license and access terms. Place user-obtained files outside the Git
repository and record their hashes locally. Do not commit them unless their
license explicitly permits redistribution.

Digitized figure data are treated as derived research evidence. This public
snapshot reports comparison metrics and methods but does not redistribute the
publisher figure or the full digitized coordinate set.

## Contributions

Every new material model or validation adapter should state:

- source DOI and exact table/figure locator;
- adsorbent grade and binder when reported;
- gas system, temperature and component partial-pressure range;
- unit conversions;
- fitted and fixed parameters;
- calibration/holdout split;
- claim boundary and downstream release flags.
