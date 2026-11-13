# Source and cleanup provenance

The study code originated in nine project modules. source_provenance.csv preserves project-relative source names and original SHA-256 hashes; these hashes identify historical originals, not the edited public copies. Current public-file hashes are in SHA256SUMS.txt.

Shared functions were recovered from a saved modelling workspace; deparsing does not preserve original closure objects. During the privacy/I/O audit, R code was parsed and reformatted. Disk reads were replaced with named in-memory input requests; bundle loading was replaced with named object injection; working-directory changes, automatic sourcing, file exports and direct filesystem side effects were disabled. Directory/path-only setup utilities were retired. Python entrypoints were rewritten as parameterized wrappers; CancerSEA now accepts supplied gene sets. Three significance-driven resampling loops were blocked, not concealed.

Detailed before/after records and the unmodified pre-cleanup snapshot are retained outside the public repository. No raw patient data or private local paths are distributed. A source-code audit does not establish that source datasets, historical figures or reported numerical results were valid. See CODE_AUDIT.md and KNOWN_LIMITATIONS.md.
