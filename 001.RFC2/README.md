# RFC2-associated prostate cancer analysis code

This code-only repository contains nine analysis modules, recovered shared functions and in-memory utilities. Input files, results, working-directory settings and automatic file imports/exports are not included.

## Modules

1. `scripts/01_discovery`: single-cell discovery and RGDP signature.
2. `scripts/02_subtyping`: multi-omics preparation and clustering.
3. `scripts/03_modeling`: prognostic model development.
4. `scripts/04_clinical`: clinical and immune associations.
5. `scripts/05_gene_screen`: model gene prioritization.
6. `scripts/06_rfc2_bulk`: RFC2 bulk analysis and supplied CancerSEA signatures.
7. `scripts/07_rfc2_spatial`: single-cell and spatial analysis.
8. `scripts/08_drug_response`: predicted drug sensitivity.
9. `scripts/09_cell_states`: cell states, TF activity and communication.

Open and evaluate R/inputs.R, R/historical_helpers.R and the chosen module's helpers.R in the R editor. Supply the required preprocessed objects, then evaluate the relevant main.R sections. No launcher changes the working directory or imports data. Required input keys and expected object descriptions are listed in docs/input_objects.csv.

```r
# After evaluating the function definitions in the editor:
prad_set_inputs(list(input_001 = prepared_object))
# prad_input retrieves the registered object. A bundle is a named list of objects.
# Replace input_001 with the key used by the relevant code section.
```

This is an analysis template, not a standalone executable pipeline. Repeated historical imports of different datasets now have distinct input keys. Correctly prepared objects must be supplied for every evaluated section. File exports were disabled so the R scripts no longer directly overwrite disk files. Python prediction uses an already prepared classifier; CellPhoneDB still requires externally provided resources/output arguments and can perform library-managed I/O. No local directories or resource filenames are configured in those wrappers.

## Code audit

Three historical loops repeatedly sampled patients until P values crossed a threshold. They are blocked with explicit errors in the public clinical script. This is a substantive scientific concern, not a privacy-only edit. Relevant comparisons and final figures need full-cohort reanalysis. See docs/CODE_AUDIT.md; do not describe the code as an unchanged reproduction of all publication results.

Missingness handling, annotation, cohort filtering, score transformations and other original scientific logic remain visible; known risks are described in docs/KNOWN_LIMITATIONS.md. Dependency inventories are under config/. No historical environment lock or complete independent validation is claimed. Rights status is described in RIGHTS.md.

## In-memory helpers

`R/historical_helpers.R` also defines `score_model(expression, coefficients)` for linear and exponentiated model scores and `audit_models(scores)` for ranking stored cohort C-index values. Neither function reads or writes files. There is no separate tools directory.
