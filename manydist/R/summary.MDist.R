#' @export
summary.MDist <- function(object, ...) {
  # call the R6 method; return invisibly to avoid echo
  object$summary(...)
  invisible(object)
}
