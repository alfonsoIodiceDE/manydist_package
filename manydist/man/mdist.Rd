\name{mdist}
\alias{mdist}
\title{Calculation of Pairwise Distances for Mixed-Type Data}
\description{
Computes pairwise distances between observations described by numeric and/or categorical attributes. The function provides options for computing independent, dependent, and practice-based distances as defined in van de Velden et al. (2024), with support for various continuous and categorical distance metrics, scaling, and commensurability adjustments.
}
\usage{
mdist(x, response = NULL, distance_cont = "manhattan", distance_cat = "tot_var_dist",
      commensurable = FALSE, scaling = "none",
      ncomp = ncol(x), threshold = NULL, preset = "custom")
}
\arguments{
  \item{x}{A dataframe or tibble containing continuous (coded as numeric), categorical (coded as factors), or mixed-type variables.}
  \item{response}{An optional factor for supervised distance calculation in categorical variables, applied only if \code{distance_cat = "supervised"}. Default is \code{NULL}.}
  \item{distance_cont}{Character string specifying the distance metric for continuous variables. Options include \code{"manhattan"} (default) and \code{"euclidean"}.}
  \item{distance_cat}{Character string specifying the distance metric for categorical variables. Options include \code{"tot_var_dist"} (default), \code{"HL"}, \code{"HLeucl"}, \code{cat_dis}, \code{mca}, \code{st_dev},  \code{"matching"}, \code{"eskin"}, \code{"iof"}, \code{"of"}, \code{"goodall_3"}, \code{"goodall_4"}, \code{"gifi_chi2"}, \code{"lin"}, \code{"var_entropy"}, \code{"var_mutability"}, \code{"supervised"}, \code{"supervised_full"}, \code{"le_and_ho"} and all the options in the package philentropy.}
  \item{commensurable}{Logical. If \code{TRUE}, the function adjusts each variable's contribution to ensure equal average influence in the overall distance. Default is \code{FALSE}.}
  \item{scaling}{Character string specifying the scaling method for continuous variables. Options include \code{"none"} (default), \code{"std"}, \code{"range"}, \code{"pc_scores"}, and \code{"robust"}.}
  \item{ncomp}{Integer specifying the number of components to retain when \code{scaling = "pc_scores"}. Default is \code{ncol(x)}.}
  \item{threshold}{Numeric value specifying the percentage of variance explained by retained components when \code{scaling = "pc_scores"}. Overrides \code{ncomp} if specified. Default is \code{NULL}.}
  \item{preset}{Character string specifying pre-defined combinations of arguments. Options include:
    \itemize{
      \item{\code{"custom"} (default)}
      \item{\code{"gower"}}
      \item{\code{"unbiased_dependent"}}
      \item{\code{"euclidean_onehot"}}
      \item{\code{"catdissim"}}
    }
  }
}
\value{
  A list with the following components:
  \item{distance_mat}{The computed pairwise distance matrix.}
}

\examples{
library(palmerpenguins)
pengmix <- palmerpenguins::penguins[complete.cases(palmerpenguins::penguins),]

# Example 1: Mixed-type dataset with 'gower' preset and commensurability
distance_matrix1 <- mdist(pengmix, preset = "gower", commensurable = TRUE)

# Example 2: Euclidean metric on one-hot encoded categorical and standardized continuous variables
distance_matrix <- mdist(pengmix, preset = "euclidean_onehot")

# Example 3: Custom preset with Manhattan for continuous and matching for categorical
distance_matrix <- mdist(pengmix, preset = "custom", distance_cont = "manhattan",
                           distance_cat = "matching", commensurable = TRUE)

# Example 4: Custom preset with Manhattan for
# standardized continuous and category dissimilarity
# for categorical
distance_matrix <- mdist(
    pengmix,
    distance_cont = "manhattan",
    distance_cat = "cat_dis",
    commensurable = FALSE,
    scaling = "std"
)

# Example 5: Euclidean for standardized continuous and HLeucl for categorical
 distance_matrix <- mdist(pengmix, distance_cont = "euclidean",
distance_cat = "HLeucl", commensurable = FALSE,
                           scaling = "std")

# Example 6: Using principal component scores with 85% variance threshold
distance_matrixPCA <- mdist(pengmix, distance_cont = "euclidean",
                              distance_cat = "tot_var_dist", scaling = "pc_scores",
                              threshold = 0.85)

# Example 7: Distance for categorical variables only
distance_matrix <- mdist(pengmix[, c(1, 2, 7)], distance_cat = "tot_var_dist")

# Example 8: Distance for continuous variables only
distance_matrix <- mdist(pengmix[, 3:6], distance_cont = "manhattan")
}

\references{
van de Velden, M., Iodice D'Enza, A., Markos, A., Cavicchia, C. (2024). (Un)biased distances for mixed-type data. \emph{arXiv preprint}. Retrieved from \url{https://arxiv.org/abs/2411.00429}.
}

\seealso{
\code{\link{cdist}}, \code{\link{ndist}}
}
