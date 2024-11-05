\name{cdist}
\alias{cdist}
\title{Calculation of Pairwise Distances for Categorical Data}
\description{
Computes a distance matrix for categorical variables with optional commensurability and custom scaling. Supports multiple methods for distance calculation and can handle different methods per variable (see, van de Velden et al., 2024).
}
\usage{
cdist(x, y = NULL, method = "tot_var_dist", commensurable = FALSE,
       weights = 1)
}
\arguments{
  \item{x}{A data frame or matrix of categorical variables.}
  \item{y}{Optional additional data for joint distance calculations with \code{x}. Default is \code{NULL}.}
  \item{method}{Character string specifying the distance metric for categorical variables. Options include \code{"tot_var_dist"} (default), \code{"HL"}, \code{"HLeucl"}, \code{cat_dis}, \code{mca}, \code{st_dev},  \code{"matching"}, \code{"eskin"}, \code{"iof"}, \code{"of"}, \code{"goodall_3"}, \code{"goodall_4"}, \code{"gifi_chi2"}, \code{"lin"}, \code{"var_entropy"}, \code{"var_mutability"}, \code{"supervised"}, \code{"supervised_full"}, \code{"le_and_ho"} and all the options in the package philentropy.}
  \item{commensurable}{Logical. If \code{TRUE}, applies commensurability adjustments to standardize each variable's contribution to the distance matrix. Default is \code{FALSE}.}
  \item{weights}{Numeric vector or matrix of weights to apply to each variable or level in the calculation. Default is \code{1} (no weighting).}
}
\details{
The \code{cdist} function computes a distance matrix tailored for categorical variables, with support for different methods for each variable. This function also provides an option for commensurability to standardize the distance contributions from each variable.

\strong{Warning:} If the length of \code{weights} differs from the number of variables, the function will either expand \code{weights} or produce an error.
}
\value{
  A list with the following components:
  \item{distance_mat}{The computed distance matrix.}
  \item{delta}{Matrix of pairwise level distances for each variable.}
  \item{delta_ms}{Any additional output relevant to the specific scaling method, if available.}
  \item{delta_names}{Names of the levels used in the \code{delta} matrix.}
}
\examples{
# Load a sample categorical dataset from the palmerpenguins package
library(palmerpenguins)
penguins_cat <- palmerpenguins::penguins[, c("species", "island", "sex")]
penguins_cat <- penguins_cat[complete.cases(penguins_cat), ]

# Example 1: Total Variation Distance
result <- cdist(penguins_cat)

# Example 2: Using commensurability with Total Variation Distance
result <- cdist(penguins_cat, commensurable = TRUE)

# Example 3: Custom weights for each variable
result <- cdist(penguins_cat, weights = c(0.5, 1.5, 1))

# Example 4: Specifying different methods for each variable
result <- cdist(penguins_cat, method = c("tot_var_dist", "matching", "goodall_3"))
}
\references{
van de Velden, M., Iodice D'Enza, A., Markos, A., Cavicchia, C. (2024). (Un)biased distances for mixed-type data. \emph{arXiv preprint}. Retrieved from \url{https://arxiv.org/abs/2411.00429}.
}

\seealso{
\code{\link{mdist}}, \code{\link{ndist}}
}
