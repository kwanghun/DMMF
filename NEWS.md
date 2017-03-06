## News for Package `DMMF`
#### Changes in DMMF version 0.2.7.9000
- Change `NAMESPACE` to avoid note messages when build in unstable developement version of R in win-builder.
- Remove line `exportPattern("^[[:alpha:]]+")`
- Add line `export("DMMF", "DMMF_Simple", "SinkFill")`
#### Changes in DMMF version 0.2.5
- Change `NAMESPACE` from `useDynLib(DMMF)` to `useDynLib(DMMF, .registration = TRUE)` to remove note messages when build in unstable developement version of R in win-builder.
#### Changes in DMMF version 0.2.4
- Add `\dontrun{}` option for `SinkFill` and `DMMF_simple` example to avoid error messages in win-builder. 
#### Changes in DMMF version 0.2.4
- Add `\dontrun{}` option for `DMMF` example to avoid error messages in win-builder. 
#### Changes in DMMF version 0.2.3
- Modifying the source code of MMF Algorithm from `SW_t = dmax1( 0.0, SW - SW_fc )` to `SW_t = dmax1( 0.0d0, SW - SW_fc )` to make 0.0 as double.
#### Changes in DMMF version 0.2.2
- Modifying NAMESPACE as import `rgdal` package (import("rgdal"))and import `is` from `methods` package (importfrom("methods", "is")).
#### Changes in DMMF version 0.2.1
- Bug fixed: All non-ASCII characters in the manuals are corrected to ASCII character to remove warnings when check the package with `R CMD CHECK --as-CRAN`.
- Bug fixed: Title and package description parts of the `/man/DMMF-package.Rd` is changed to remove warnings when check the package with `R CMD CHECK --as-CRAN`.
- README.md and NEWS.md files are added to the package.
#### Changes in DMMF version 0.2.0
- The new function named `DMMF_Simple` is added to the package which is the simpler version of `DMMF` function.
- The `/src/2_2_BoundaryChecker.f95` is changed to convert logical matrix to real value matrix of 1 and 0 smoothly without warning message of "implicit conversion from logical to integer".
#### Changes in DMMF version 0.1.0
- The manual of functions in the package are added to the package.
- Datasets from two potato fields (`Potato.Concave` and `Potato.Convex`) are added to the package.


