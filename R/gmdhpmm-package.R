#' gmdhpmm: GMDH with PMM coefficient estimation
#'
#' MIA-GMDH whose inner coefficient step dispatches automatically between LSE,
#' PMM2 and PMM3 based on bootstrap-stabilized residual cumulants. See
#' \code{\link{gmdh_pmm}} for the entry point and \code{../algorithm-spec.md}
#' for the design.
#'
#' @keywords internal
"_PACKAGE"

# NULL-coalescing helper used throughout.
`%||%` <- function(a, b) if (is.null(a)) b else a
