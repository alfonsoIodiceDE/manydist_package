#' @export
print.MDistLOVO <- function(x, ...) x$print(...)

#' @export
summary.MDistLOVO <- function(object, ...) object$summary(...)

# If you have $autoplot():
#' @importFrom ggplot2 autoplot
#' @exportS3Method ggplot2::autoplot MDistLOVO
autoplot.MDistLOVO <- function(object, ...) object$autoplot(...)
