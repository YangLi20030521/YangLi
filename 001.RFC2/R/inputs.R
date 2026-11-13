# In-memory input registry. No disk, directory, or environment-variable access.
.prad_inputs <- new.env(parent = emptyenv())
prad_set_inputs <- function(values) {
    if (!is.list(values) || is.null(names(values)) || any(names(values) == ""))
        stop("Supply a named list of input objects")
    for (key in names(values)) assign(key, values[[key]], envir = .prad_inputs)
    invisible(NULL)
}
prad_input <- function(key) {
    if (!exists(key, envir = .prad_inputs, inherits = FALSE))
        stop(paste("Required in-memory input is missing:", key))
    get(key, envir = .prad_inputs, inherits = FALSE)
}
prad_restore <- function(key) {
    bundle <- prad_input(key)
    if (!is.list(bundle) || is.null(names(bundle)) || any(names(bundle) == ""))
        stop("An object bundle must be a named list")
    list2env(bundle, envir = parent.frame())
    invisible(NULL)
}
