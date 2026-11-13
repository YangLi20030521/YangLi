# Work with gene sets already supplied in memory. No downloads or file reads.
prepare_signatures <- function(signature_sets) {
    states <- c("Angiogenesis", "Apoptosis", "Cell Cycle", "Differentiation", "DNA damage", "DNA repair",
                "EMT", "Hypoxia", "Inflammation", "Invasion", "Metastasis", "Proliferation", "Quiescence", "Stemness")
    if (!is.list(signature_sets) || is.null(names(signature_sets))) stop("Supply named gene sets")
    lapply(signature_sets[intersect(states, names(signature_sets))], function(x) unique(x[!is.na(x) & x != ""]))
}
