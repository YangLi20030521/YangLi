# Validation scope

All current R files were syntax-parsed and Python files AST-parsed after cleanup. Pure scoring and input-registry functions were checked with synthetic in-memory objects. Source scans cover direct local file readers/writers, directory setup, private path patterns and common embedded secret literals. No whole-stage analysis, original patient data, final figure, or library-internal behavior was validated. Three significance-selected subsampling loops are intentionally blocked. See docs/CODE_AUDIT.md.
