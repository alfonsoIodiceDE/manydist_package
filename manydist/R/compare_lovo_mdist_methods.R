#' @export
print.MDistLOVOCompare <- function(x, ...) x$print(...)

#' @export
summary.MDistLOVOCompare <- function(object, ...) object$summary(...)

#' @importFrom ggplot2 autoplot
#' @exportS3Method ggplot2::autoplot MDistLOVOCompare
autoplot.MDistLOVOCompare <- function(object, ...) object$autoplot(...)
