#' generate toy mixed datasets for the supplementary material: figures 1 to 3
#'
#' @param n number of observations
#' @param porig number of original informative continuous variables
#' @param pn number of continuous variables (total, before adding noise variables)
#' @param pnnoise number of extra numeric noise variables
#' @param pcnoise number of extra categorical noise variables
#' @param sigma sd for noise added to informative variables
#' @param qoptions number of bins for categorization (vector if per_variable, scalar if shared)
#' @param seed optional seed
#' @param mode either "per_variable" or "shared"
#' @export
generate_dataset <- function(n, porig, pn, pnnoise, pcnoise, sigma, qoptions, seed = NULL,
                             mode = "per_variable") {
  # mode <- match.arg(mode)
  if (!is.null(seed)) set.seed(seed)

  # Generate latent 2D structure
  T <- cbind(runif(n, -2, 2), runif(n, -2, 2))

  SVDT <- svd(T)
  Y <- SVDT$u



  # Create informative projection
  R <- matrix(runif(2 * porig, -2, 2), nrow = 2)


  X <- as.data.frame(Y %*% R)


  colnames(X) <- paste0("X", seq_len(porig))

  # Xorig without noise

  Xorig <- X


  # Add noise to true variables
  X[1:porig] <- purrr::map_df(X[1:porig], ~ .x + rnorm(n, 0, sigma))


  # Xorig with noise

  R2_Xorig <- X




  # Add numerical and categorical noise
  N <- tibble::tibble(.rows = n)

  if (pnnoise > 0) {
    num_noise <- purrr::map_dfc(
      seq_len(pnnoise),
      ~ tibble::tibble(!!paste0("N", .x) := rnorm(n))
    )
    N <- dplyr::bind_cols(N, num_noise)
  }

  if (pcnoise > 0) {
    qrand <- sample(3:9, 1)
    cat_noise <- purrr::map_dfc(
      seq_len(pcnoise),
      ~ tibble::tibble(!!paste0("C", .x) := as.factor(round(runif(n, 1.5, qrand + 0.5))))
    )
    N <- dplyr::bind_cols(N, cat_noise)
  }

  if (ncol(N) > 0) {
    X     <- dplyr::bind_cols(N, X)
    Xorig <- dplyr::bind_cols(N, Xorig)
  }
  # Categorize: determine which variables to cut
  p <- ncol(X)
  pnum <- pn + pnnoise
  pcat <- p - pnum


  # --- Here is the conditional behavior ---
  if (mode == "per_variable") {
    stopifnot(length(qoptions) == pcat)
    for (i in (pnum + 1):p) {
      q <- qoptions[i - pnum]
      X[[i]] <- cut(X[[i]], breaks = q, labels = FALSE) |> as.factor()
    }
  } else if (mode == "shared") {
    stopifnot(length(qoptions) == 1)
    q <- qoptions
    for (i in (pnum + 1):p) {
      X[[i]] <- cut(X[[i]], breaks = q, labels = FALSE) |> as.factor()
    }
  }

  return(list(X = X, Xorig = Xorig,R2_Xorig=R2_Xorig, truth=Y))
}
