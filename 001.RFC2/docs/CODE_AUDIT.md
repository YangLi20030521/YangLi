# Privacy and data-integrity code audit

## Privacy and filesystem operations

The public scripts no longer configure a personal working directory or a local input location. No direct email address, local account-directory path, embedded authentication token or password literal was found in the executable-source scan. This is a source-code review, not a guarantee about unseen data files or upstream packages. GEO accession/sample IDs are public scientific identifiers; they are retained.

R disk input calls were replaced with named in-memory objects. Automatic file writes, save.image, plot exports, directory operations and source calls were disabled. Plotting/analysis code remains, but a caller may need to manage graphical devices. Third-party packages can access public network resources internally; the CellPhoneDB wrapper explicitly delegates resource and output handling to its library using caller-supplied arguments.

## Additional transformations requiring disclosure

- Missingness handling includes complete-case exclusion, gene-mean imputation and, inside univariate Cox screening, NA-to-zero replacement.
- Discovery annotation manually reassigns selected clusters and Cycling_T labels to Prolif; malignant-cell analysis excludes GSM4089153. These require scientific justification and provenance.
- TCGA PFS and other-cohort RFS were renamed OS during modelling, and original time values were multiplied by 365; reconcile endpoint and unit semantics.
- The clinical response workflow uses cohort-standardized linear scores, whereas stored prognostic scores use exp(X beta). Their scales and interpretation differ.
- Mean cross-cohort model selection, within-cohort Cox direction refitting, cohort-optimized cutoffs, nominal P filtering and selected pathway displays require the limitations documented in KNOWN_LIMITATIONS.md.

These operations modify analysis objects and can alter inference; they should not all be called falsification. No complete forensic guarantee or validation of published results is made. Original scientific transformations other than the three blocked resampling loops were not silently repaired.
