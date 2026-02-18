###  numsep =.1 and catsep=0.6 generation


```{r}
nrep <- 100
k_true <- 4
q <- 9
numsep <- 0.1
catsep <- 0.6
clustSizeEq <- 50

param_generator_grid <- tibble::tibble(
  numsignal = c(0,0,4,4,4,8,8,8),
  numnoise  = c(8,8,4,4,4,0,0,0),
  catsignal = c(4,8,0,4,8,0,4,8),
  catnoise  = c(4,0,8,4,0,8,4,0)
) |>
  tidyr::crossing(replicate = 1:nrep)

mixed_datasets_grid <- param_generator_grid |>
  dplyr::mutate(
    data = purrr::pmap(
      dplyr::pick(numsignal, catsignal, numnoise, catnoise, replicate),
      \(numsignal, catsignal, numnoise, catnoise, replicate) {
        gen_mixed(
          k_true      = k_true,
          clustSizeEq = clustSizeEq,
          numsignal   = numsignal,
          numnoise    = numnoise,
          catsignal   = catsignal,
          catnoise    = catnoise,
          q           = q,
          q_err       = q,
          numsep      = numsep,
          catsep      = catsep,
          seed        = replicate
        )
      }
    )
  )

methods <- c("naive","hl","hl_add","gower","gudmm","dkss","mod_gower","u_ind","u_dep","u_mix")

sim_methods_grid <- mixed_datasets_grid |>
  tidyr::crossing(method = methods) |>
  dplyr::left_join(manydist:::mdist_method_lookup, by = "method")

stopifnot(!any(is.na(sim_methods_grid$mdist_type)))

ari_pam_results_01_06 <- sim_methods_grid |>
  dplyr::mutate(
    ari = purrr::pmap_dbl(
      dplyr::pick(data, mdist_type, mdist_preset, param_set),
      \(data, mdist_type, mdist_preset, param_set) {
        df <- data$df
        D <- run_mdist_within(df, mdist_type, mdist_preset, param_set, outcome = "y")
        pam_fit <- cluster::pam(D, k = k_true, diss = TRUE)
        mclust::adjustedRandIndex(df$y, pam_fit$clustering)
      }
    )
  )

# ari_pam_results_01_06 = ari_pam_results

save(ari_pam_results_01_06, file = "ari_pam_experiment_numsep_01_catsep_06.RData")
```



###  numsep =.01 and catsep=0.3 generation


```{r}
nrep <- 100
k_true <- 4
q <- 9
numsep <- 0.01
catsep <- 0.3
clustSizeEq <- 50

param_generator_grid <- tibble::tibble(
  numsignal = c(0,0,4,4,4,8,8,8),
  numnoise  = c(8,8,4,4,4,0,0,0),
  catsignal = c(4,8,0,4,8,0,4,8),
  catnoise  = c(4,0,8,4,0,8,4,0)
) |>
  tidyr::crossing(replicate = 1:nrep)

mixed_datasets_grid <- param_generator_grid |>
  dplyr::mutate(
    data = purrr::pmap(
      dplyr::pick(numsignal, catsignal, numnoise, catnoise, replicate),
      \(numsignal, catsignal, numnoise, catnoise, replicate) {
        gen_mixed(
          k_true      = k_true,
          clustSizeEq = clustSizeEq,
          numsignal   = numsignal,
          numnoise    = numnoise,
          catsignal   = catsignal,
          catnoise    = catnoise,
          q           = q,
          q_err       = q,
          numsep      = numsep,
          catsep      = catsep,
          seed        = replicate
        )
      }
    )
  )

methods <- c("naive","hl","hl_add","gower","gudmm","dkss","mod_gower","u_ind","u_dep","u_mix")

sim_methods_grid <- mixed_datasets_grid |>
  tidyr::crossing(method = methods) |>
  dplyr::left_join(manydist:::mdist_method_lookup, by = "method")

stopifnot(!any(is.na(sim_methods_grid$mdist_type)))

ari_pam_results_001_03 <- sim_methods_grid |>
  dplyr::mutate(
    ari = purrr::pmap_dbl(
      dplyr::pick(data, mdist_type, mdist_preset, param_set),
      \(data, mdist_type, mdist_preset, param_set) {
        df <- data$df
        D <- run_mdist_within(df, mdist_type, mdist_preset, param_set, outcome = "y")
        pam_fit <- cluster::pam(D, k = k_true, diss = TRUE)
        mclust::adjustedRandIndex(df$y, pam_fit$clustering)
      }
    )
  )


save(ari_pam_results_001_03, file = "ari_pam_experiment_numsep_001_catsep_03.RData")
```


