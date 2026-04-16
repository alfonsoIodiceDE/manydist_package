cat_delta <- function(x, y = NULL, method = NULL, method_cat = "tvd", mkw_p = 1) {

  a = NULL
  b = NULL
  blocks = NULL
  .x = NULL
  id = NULL

  is_scalar_vector <- function(x) {
    is.atomic(x) && length(x) == 1
  }

  method = method_cat
  catdiss = method
  x = purrr::map_df(x, forcats::fct_drop)

  if (is.null(dim(x))) {
    Q = nlevels(x)
  } else {
    Q = purrr::map_dbl(x, nlevels)
  }

  n = nrow(x)
  nvar = length(Q)

  z_prep = z_preproc(x = x, y = y, Q = Q)

  Z_names = colnames(z_prep$Z)
  Z = z_prep$Z %>% data.matrix()
  ZZod = z_prep$ZZod
  zm = z_prep$zm
  Z_list = z_prep$Z_list

  level_pos = data.table::data.table(
    start = c(1, cumsum(Q)[-length(Q)] + 1),
    stop  = cumsum(Q)
  )

  if (!is.null(y)) {
    Z_y = z_prep$Z_y
  } else {
    Z_y = NULL
  }

  Qs = ncol(Z)

  if (method %in% c("tvd", "tot_var_dist")) {
    method <- "tvd"
  }

  response_irrelevant_methods <- c(
    "matching", "eskin", "goodall_3", "goodall_4",
    "iof", "of", "lin", "var_entropy", "var_mutability"
  )

  custom_profile_methods <- c("tvd", "gifi_chi2")

  profile_based_methods <- c(
    custom_profile_methods,
    "le_and_ho",
    philentropy::getDistMethods()
  )

  if (method %in% c(response_irrelevant_methods, custom_profile_methods)) {

    full_delta = cat_custom_delta(
      ZZod = ZZod, Z = Z, Z_y = Z_y, Z_list = Z_list,
      zm = zm, Q = Q, nvar = nvar, method = method, Qs = Qs
    )


  } else {

    crs = tidyr::crossing(a = 1:nvar, b = 1:nvar) %>%
      dplyr::filter(a != b)

    blocks_id_a = level_pos[crs$a, ]
    blocks_id_b = level_pos[crs$b, ] %>%
      dplyr::rename(end_start = start, end_stop = stop)

    block_ids = cbind(blocks_id_a, blocks_id_b)

    pull_block <- function(start = 1, stop = 1, end_start = 1, end_stop = 1) {
      ZZod[start:stop, end_start:end_stop]
    }

    pull_conditional_block <- function(start, stop, end_start, end_stop) {
      joint_block <- ZZod[start:stop, end_start:end_stop]
      conditional_block <- sweep(joint_block, 2, colSums(joint_block), FUN = "/")
      conditional_block[is.nan(conditional_block)] <- 0
      conditional_block[is.infinite(conditional_block)] <- 0
      conditional_block
    }

    # metrics that really want probability vectors
    prob_methods <- c(
      "hellinger", "kullback-leibler", "jensen-shannon", "bhattacharyya"
    )

    sanitize_block <- function(block, method_name) {
      block <- base::as.matrix(block)

      # basic cleanup
      block[is.na(block)] <- 0
      block[!is.finite(block)] <- 0

      # numerical tolerance: tiny negatives -> 0
      block[block < 0 & block > -1e-12] <- 0

      prob_methods <- c(
        "hellinger", "kullback-leibler", "jensen-shannon", "bhattacharyya"
      )

      if (method_name %in% prob_methods) {
        rs <- rowSums(block)

        # if a row is all zero, replace with uniform probabilities
        zero_rows <- rs <= .Machine$double.eps
        if (any(zero_rows)) {
          block[zero_rows, ] <- 1 / ncol(block)
          rs <- rowSums(block)
        }

        # normalize rows to sum to 1
        block <- sweep(block, 1, rs, "/")

        # final cleanup
        block[is.na(block)] <- 0
        block[!is.finite(block)] <- 0
      }

      block
    }

    run_philentropy <- function(block, method_name, p = 1) {
      tmp <- NULL

      suppressMessages(
        suppressWarnings(
          invisible(capture.output(
            capture.output(
              tmp <- philentropy::distance(
                x = block,
                method = method_name,
                mute.message = TRUE,
                p = p
              ),
              type = "message"
            )
          ))
        )
      )

      out <- base::as.matrix(tmp)

      if (length(out) == 1) {
        out <- matrix(as.numeric(out), 2, 2)
      }

      off_diag <- out
      diag(off_diag) <- 0

      if (anyNA(off_diag) || any(!is.finite(off_diag))) {
        stop(
          "philentropy produced non-finite off-diagonal values for method '",
          method_name, "'.",
          call. = FALSE
        )
      }

      diag(out) <- 0
      out
    }

    get_profile_blocks <- function(method_name, supervised = FALSE) {
      if (!supervised) {
        if (method_name == "le_and_ho") {
          tibble::tibble(
            row_ind = crs$a,
            col_ind = crs$b,
            blocks = purrr::pmap(
              block_ids,
              ~ pull_conditional_block(
                start = ..1, stop = ..2,
                end_start = ..3, end_stop = ..4
              )
            )
          )
        } else {
          tibble::tibble(
            row_ind = crs$a,
            col_ind = crs$b,
            blocks = purrr::pmap(
              block_ids,
              ~ pull_block(
                start = ..1, stop = ..2,
                end_start = ..3, end_stop = ..4
              )
            )
          )
        }
      } else {
        if (is.null(Z_y)) {
          stop("Supervised profile blocks requested but `Z_y` is NULL.", call. = FALSE)
        }

        prof_sup <- sweep(t(Z) %*% Z_y, 1, zm, "/")
        prof_sup[is.na(prof_sup)] <- 0
        prof_sup[is.infinite(prof_sup)] <- 0

        tibble::tibble(
          row_ind = seq_len(nvar),
          col_ind = NA_integer_,
          blocks = purrr::map(
            seq_len(nvar),
            function(i) {
              prof_sup[level_pos$start[i]:level_pos$stop[i], , drop = FALSE]
            }
          )
        )
      }
    }

    supervised <- !is.null(y)
    distance_blocks <- get_profile_blocks(method_name = catdiss, supervised = supervised)

    if (catdiss == "le_and_ho") {

      distance_blocks <- distance_blocks %>%
        dplyr::mutate(
          block_dist = purrr::map(
            .x = blocks,
            .f = function(x) {
              x <- sanitize_block(x, "kullback-leibler")

              dist_long <- tidyr::crossing(a = 1:nrow(x), b = 1:nrow(x)) %>%
                dplyr::filter(a != b) %>%
                dplyr::mutate(
                  kl = purrr::map2_dbl(
                    .x = a, .y = b,
                    ~ philentropy::kullback_leibler_distance(
                      P = x[.x, ],
                      Q = x[.y, ],
                      unit = "log2",
                      testNA = FALSE,
                      epsilon = 1e-8
                    )
                  )
                )

              phil_dist <- matrix(0, nrow(x), nrow(x))
              for (i in 1:nrow(dist_long)) {
                phil_dist[dist_long$a[i], dist_long$b[i]] <- dist_long$kl[i]
              }
              phil_dist <- phil_dist + t(phil_dist)
              diag(phil_dist) <- 0

              if (anyNA(phil_dist) || any(!is.finite(phil_dist))) {
                stop("`le_and_ho` produced non-finite values.", call. = FALSE)
              }

              phil_dist
            }
          )
        )

    } else {

      distance_blocks <- distance_blocks %>%
        dplyr::mutate(
          block_dist = purrr::map(
            .x = blocks,
            .f = function(x) {
              x <- sanitize_block(x, catdiss)
              run_philentropy(x, catdiss, p = mkw_p)
            }
          )
        )
    }

    if (!supervised) {

      delta_blocks <- tibble::tibble(id = as.list(1:nvar)) %>%
        dplyr::mutate(
          diag_delta = purrr::map(
            .x = id,
            .f = ~{
              mats <- distance_blocks %>%
                dplyr::filter(row_ind == .x) %>%
                dplyr::pull(block_dist)

              dims <- vapply(mats, function(m) paste(dim(m), collapse = "x"), character(1))
              if (length(unique(dims)) != 1) {
                stop("Non-conformable philentropy block distances for variable ", .x, ".", call. = FALSE)
              }

              Reduce(`+`, mats)
            }
          )
        )

      full_delta <- Matrix::bdiag(delta_blocks$diag_delta) %>% as.matrix()
      full_delta <- full_delta / (nvar - 1)

    } else {

      full_delta <- Matrix::bdiag(distance_blocks$block_dist) %>% as.matrix()

    }
}
  out = list()
  out$delta_names = Z_names
  out[[method_cat]] = full_delta %>% as.matrix()
  out$Z = Z
  out
}
