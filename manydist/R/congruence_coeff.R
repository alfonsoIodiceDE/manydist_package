#' Congruence coefficient between two configurations
#'
#' Computes the congruence coefficient between two data configurations
#' using the Frobenius inner product of their pairwise distance matrices.
#'
#' @param L1 A numeric matrix or data frame (rows = observations)
#' @param L2 A numeric matrix or data frame with the same number of rows as \code{L1}
#'
#' @return A scalar in \eqn{[-1,1]} measuring similarity between the two configurations.
#'
#' @details
#' The congruence coefficient is defined as
#' \deqn{
#'   \frac{\langle D_1, D_2 \rangle_F}
#'        {\sqrt{\langle D_1, D_1 \rangle_F \langle D_2, D_2 \rangle_F}}
#' }
#' where \eqn{D_1} and \eqn{D_2} are the pairwise distance matrices derived from
#' \code{L1} and \code{L2}.
#'
#' @export
congruence_coeff <- function(L1, L2) {
  L1 <- stats::dist(L1) |> as.matrix()
  L2 <- stats::dist(L2) |> as.matrix()

  num <- sum(diag(t(L1) %*% L2))
  den <- sqrt(
    sum(diag(t(L1) %*% L1)) *
      sum(diag(t(L2) %*% L2))
  )

  num / den
}
