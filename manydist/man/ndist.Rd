\name{ndist}
\alias{ndist}
\title{Calculation of Pairwise Distances for Continuous Data}
\description{
Computes a distance matrix with options for scaling, weighting, and dimensionality reduction using principal component analysis (see, van de Velden et al., 2024).
}
\usage{
ndist(x, commensurable = FALSE, method = "manhattan",
      scaling = "none", ncomp = ncol(x), threshold = NULL,
      weights = rep(1, ncol(x)))
}
\arguments{
  \item{x}{A data frame or matrix of input variables.}
  \item{commensurable}{Logical. If \code{TRUE}, standardizes each variable's distance matrix by dividing by its mean distance. Default is \code{FALSE}.}
  \item{method}{Character string specifying the distance metric to use (e.g., \code{"manhattan"} or \code{"euclidean"}). Default is \code{"manhattan"}.}
  \item{scaling}{Character string specifying the scaling method to use on variables. Options include \code{"none"}, \code{"std"} for standard scaling, \code{"range"} for min-max scaling, \code{"pc_scores"} for PCA-based dimensionality reduction, and \code{"robust"} for robust scaling. Default is \code{"none"}.}
  \item{ncomp}{Number of components to retain when \code{scaling = "pc_scores"}. Default is the number of columns in \code{x}.}
  \item{threshold}{Proportion of variance to retain when \code{scaling = "pc_scores"}. If specified, overrides \code{ncomp}. Default is \code{NULL}.}
  \item{weights}{Numeric vector of weights for each variable. Default is a vector of ones with length equal to the number of variables in \code{x}.}
}
\details{
  The \code{ndist} function calculates a distance matrix for the variables in \code{x} based on a specified distance \code{method} and optional scaling transformations. Different weighting and distance calculation options are available depending on the \code{commensurable} parameter.

  \strong{Warning:} A warning message is issued if \code{weights} is not specified by the user and differs from a vector of all ones. Additionally, if the length of \code{weights} differs from the number of variables (\code{nvars}), an error is raised.
}
\value{
A distance matrix where each entry represents the distance between variables in the data according to the specified method and scaling.
}

\examples{
library(palmerpenguins)
penguins_cont <- palmerpenguins::penguins[, c("bill_length_mm",
"bill_depth_mm", "flipper_length_mm", "body_mass_g")]

penguins_cont <- penguins_cont[complete.cases(penguins_cont), ]

# Example 1: Manhattan Distance
dist_matrix <- ndist(penguins_cont)

# Example 2: Using commensurability
dist_matrix <- ndist(penguins_cont, commensurable = TRUE)

# Example 3: Custom weights with matching length to variables
dist_matrix <- ndist(penguins_cont, weights = c(0.5, 1, 1.5, 2))

# Example 4: Specifying "pc_scores" scaling with number of components
dist_matrix <- ndist(penguins_cont, scaling = "pc_scores", ncomp = 2)

# Example 5: Using threshold with PCA-based scaling
dist_matrix <- ndist(penguins_cont, scaling = "pc_scores", threshold = 0.9)

# Example 6: Manhattan distance with robust scaling
dist_matrix <- ndist(penguins_cont, scaling = "robust")

# Example 7: Euclidean distance with custom weights
dist_matrix <- ndist(penguins_cont, method = "euclidean", weights = c(1, 0.5, 2, 1))

}
\references{
van de Velden, M., Iodice D'Enza, A., Markos, A., Cavicchia, C. (2024). (Un)biased distances for mixed-type data. \emph{arXiv preprint}. Retrieved from \url{https://arxiv.org/abs/2411.00429}.
}

\seealso{
\code{\link{mdist}}, \code{\link{cdist}}
}
