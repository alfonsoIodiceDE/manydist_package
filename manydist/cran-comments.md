## Test environments

- local macOS Tahoe 26.5.1, R 4.6.0
- win-builder, R-devel
- win-builder, R-release

## R CMD check results

0 errors | 0 warnings | 1 note

## Notes

The CRAN incoming feasibility note concerns the change of maintainer from Angelos Markos to Alfonso Iodice D'Enza. This change is intentional and agreed among the package authors.

The spell-check note concerns author names and a proper name/technical term: "Cavicchia", "D'Enza", "Iodice", "Markos", and "Gower".

## Maintainer change

The maintainer has changed from Angelos Markos to Alfonso Iodice D'Enza.

Angelos Markos has confirmed the maintainer change to CRAN from the previous maintainer email address.

## Update

This is an update to version 0.5.0.

Main changes:

- Renamed `distance_cat` to `method_cat`.
- Renamed `scaling_cont` to `method_num`.
- Removed the user-facing `distance_cont` argument; the numerical metric is now determined by the selected preset.
- Renamed the `"euclidean_onehot"` preset to `"euclidean"`.
- Updated documentation, examples, and method-specification helpers accordingly.

The long-form vignette is not included in this CRAN submission.
