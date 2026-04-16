cat_custom_delta <- function(ZZod, Z, Z_y, Z_list, zm, Q, nvar, method, Qs) {

  .x <- NULL

  n <- nrow(Z)

  if (method %in% c("tvd", "tot_var_dist")) {
    method <- "tvd"
  }

  # ------------------------------------------------------------------
  # helpers
  # ------------------------------------------------------------------

  Qlist <- split(Q, seq_along(Q))
  blocks <- Matrix::bdiag(purrr::map(.x = Qlist, ~ matrix(1, nrow = .x, ncol = .x)))

  create_delta <- function(profile_mat, n_profiles, p, chi2 = FALSE) {
    if (!chi2) {
      if (n_profiles != 1) {
        Dsel <- as.matrix(stats::dist(profile_mat, method = "minkowski", p = p)) /
          ((3 - p) * (n_profiles - 1))
      } else {
        Dsel <- as.matrix(stats::dist(profile_mat, method = "minkowski", p = p)) /
          ((3 - p) * n_profiles)
      }
    } else {
      if (n_profiles != 1) {
        Dsel <- as.matrix(stats::dist(profile_mat, method = "minkowski", p = p)^2) /
          ((3 - p) * (n_profiles - 1))
      } else {
        Dsel <- as.matrix(stats::dist(profile_mat, method = "minkowski", p = p)) /
          ((3 - p) * n_profiles)
      }
    }
    Dsel
  }

  get_supervised_profiles <- function() {
    if (is.null(Z_y)) {
      stop("Supervised profiles requested but `Z_y` is NULL.", call. = FALSE)
    }
    sweep(t(Z) %*% Z_y, 1, zm, "/")
  }

  get_unsupervised_profiles <- function() {
    ZZod
  }

  # ------------------------------------------------------------------
  # response-irrelevant methods
  # ------------------------------------------------------------------

  if (method == "matching") {
    full_delta <- as.matrix(blocks)
    diag(full_delta) <- 0
    return(full_delta)
  }

  if (method == "eskin") {
    full_delta <- Matrix::bdiag(
      purrr::map(.x = Z_list, .f = function(x) {
        Qi <- ncol(x)
        sk <- Qi^2 / (Qi^2 + 2)
        matrix((1 / sk - 1), nrow = Qi, ncol = Qi)
      })
    ) |>
      as.matrix()

    diag(full_delta) <- 0
    return(full_delta)
  }

  if (method == "goodall_3") {
    full_delta <- (blocks * (
      rep(1, Qs, Qs) -
        diag(1 - (zm * (zm - 1)) / (n * (n - 1)), nrow = Qs, ncol = Qs)
    )) / nvar
    return(as.matrix(full_delta))
  }

  if (method == "goodall_4") {
    full_delta <- (blocks * (
      rep(1, Qs, Qs) -
        diag(((zm * (zm - 1)) / (n * (n - 1))), nrow = Qs, ncol = Qs)
    )) / nvar
    return(as.matrix(full_delta))
  }

  if (method == "iof") {
    full_delta <- 1 / (1 + (log(zm) %*% t(log(zm))))
    diag(full_delta) <- 1
    full_delta <- (blocks * (1 / full_delta - 1)) / nvar
    return(as.matrix(full_delta))
  }

  if (method == "of") {
    full_delta <- 1 / (1 + (log(n / zm) %*% t(log(n / zm))))
    diag(full_delta) <- 1
    full_delta <- (blocks * (1 / full_delta - 1)) / nvar
    return(as.matrix(full_delta))
  }

  if (method == "lin") {
    prop <- as.matrix(zm / n)

    pv <- prop
    pr <- pv %*% rep(1, Qs)
    pc <- rep(1, Qs) %*% t(pv)
    pp <- pr + pc
    pplog <- log(pr) + log(pc)
    diag(pp) <- prop

    linsim2 <- (pplog - 2 * log(pp)) / 2 * log(pp)
    linsim2[is.nan(linsim2)] <- 0

    full_delta <- as.matrix(blocks * linsim2)

    start_index <- 1
    for (j in seq_along(Q)) {
      if (Q[j] == 2) {
        end_index <- start_index + 1
        simple_matching_block <- matrix(1, nrow = 2, ncol = 2) - diag(2)
        full_delta[start_index:end_index, start_index:end_index] <- simple_matching_block
      }
      start_index <- start_index + Q[j]
    }

    return(full_delta)
  }

  if (method == "var_entropy") {
    prop <- as.matrix(zm / n)
    full_delta <- matrix(0, Qs, Qs)
    pos <- 0

    for (i in seq_len(nvar)) {
      sk <- -1 / log(Q[i])
      plogp <- prop[(pos + 1):(pos + Q[i])] %*% log(prop[(pos + 1):(pos + Q[i])])
      full_delta[(pos + 1):(pos + Q[i]), (pos + 1):(pos + Q[i])] <- 1 - sk * diag(rep(plogp, Q[i]))
      pos <- pos + Q[i]
    }

    full_delta <- full_delta / nvar
    return(full_delta)
  }

  if (method == "var_mutability") {
    prop <- as.matrix(zm / n)
    full_delta <- matrix(0, Qs, Qs)
    pos <- 0

    for (i in seq_len(nvar)) {
      sk2 <- Q[i] / (Q[i] - 1)
      pp <- 1 - prop[(pos + 1):(pos + Q[i])] %*% prop[(pos + 1):(pos + Q[i])]
      full_delta[(pos + 1):(pos + Q[i]), (pos + 1):(pos + Q[i])] <- 1 - sk2 * diag(rep(pp, Q[i]))
      pos <- pos + Q[i]
    }

    full_delta <- full_delta / nvar
    return(full_delta)
  }

  # ------------------------------------------------------------------
  # profile-based custom methods
  # ------------------------------------------------------------------

  if (method == "tvd") {
    profile_mat <- if (is.null(Z_y)) get_unsupervised_profiles() else get_supervised_profiles()
    n_profiles  <- if (is.null(Z_y)) nvar else 3 / 2

    full_delta <- blocks * create_delta(profile_mat, n_profiles = n_profiles, p = 1, chi2 = FALSE)
    return(as.matrix(full_delta))
  }

  if (method == "gifi_chi2") {
    if (is.null(Z_y)) {
      profile_mat <- t(t(get_unsupervised_profiles()) / sqrt(zm / n))
      n_profiles  <- nvar
    } else {
      sup_prof <- get_supervised_profiles()
      sup_marg <- rowSums(t(Z) %*% Z_y)
      profile_mat <- t(t(sup_prof) / sqrt(sup_marg / n))
      n_profiles  <- 3 / 2
    }

    full_delta <- blocks * create_delta(profile_mat, n_profiles = n_profiles, p = 2, chi2 = TRUE)
    return(as.matrix(full_delta))
  }

  stop("Unsupported custom categorical method: ", method, call. = FALSE)
}
