
```{r}
x_cat <- purrr::map_df(x_cat, forcats::fct_drop)

Q <- purrr::map_dbl(x_cat, nlevels)

z_prep <- z_preproc(x = x_cat, y = NULL, Q = Q)

ZZod <- z_prep$ZZod

level_pos <- data.table::data.table(
  start = c(1, cumsum(Q)[-length(Q)] + 1),
  stop  = cumsum(Q)
)

crs <- tidyr::crossing(a = 1:length(Q), b = 1:length(Q)) |>
  dplyr::filter(a != b)

blocks_id_a <- level_pos[crs$a, ]
blocks_id_b <- level_pos[crs$b, ] |>
  dplyr::rename(end_start = start, end_stop = stop)

block_ids <- cbind(blocks_id_a, blocks_id_b)
```

```{r}
test_block <- ZZod[
  block_ids$start[1]:block_ids$stop[1],
  block_ids$end_start[1]:block_ids$end_stop[1]
]

test_block
dim(test_block)
```

```{r}
tmp <- philentropy::distance(
  x = test_block,
  method = "euclidean",
  mute.message = TRUE,
  p = 1
)

tmp
class(tmp)
length(tmp)
dim(tmp)

base::as.matrix(tmp)
dim(base::as.matrix(tmp))
```

```{r}

x_cat <- purrr::map_df(x_cat, forcats::fct_drop)

Q <- purrr::map_dbl(x_cat, nlevels)
nvar <- length(Q)

z_prep <- z_preproc(x = x_cat, y = NULL, Q = Q)

ZZod <- z_prep$ZZod

level_pos <- data.table::data.table(
  start = c(1, cumsum(Q)[-length(Q)] + 1),
  stop  = cumsum(Q)
)

crs <- tidyr::crossing(a = 1:nvar, b = 1:nvar) |>
  dplyr::filter(a != b)

blocks_id_a <- level_pos[crs$a, ]
blocks_id_b <- level_pos[crs$b, ] |>
  dplyr::rename(end_start = start, end_stop = stop)

block_ids <- cbind(blocks_id_a, blocks_id_b)

pull_block <- function(start = 1, stop = 1, end_start = 1, end_stop = 1) {
  ZZod[start:stop, end_start:end_stop]
}

distance_blocks <- tibble::tibble(
  row_ind = crs$a,
  col_ind = crs$b,
  blocks = purrr::pmap(
    block_ids,
    ~ pull_block(start = ..1, stop = ..2, end_start = ..3, end_stop = ..4)
  )
)

catdiss <- "euclidean"
mkw_p <- 1

distance_blocks_checked <- distance_blocks |>
  dplyr::mutate(
    block_dist = purrr::map(
      blocks,
      function(x) {
        x[is.na(x)] <- 0
        tmp <- philentropy::distance(
          x = x,
          method = catdiss,
          mute.message = TRUE,
          p = mkw_p
        )
        base::as.matrix(tmp)
      }
    ),
    nr = purrr::map_int(block_dist, nrow),
    nc = purrr::map_int(block_dist, ncol),
    cls = purrr::map_chr(block_dist, ~ paste(class(.x), collapse = ","))
  )

distance_blocks_checked |>
  dplyr::select(row_ind, col_ind, cls, nr, nc) |> print(n=30)


```

```{r}
catdiss <- "kullback-leibler"

db <- distance_blocks |>
  dplyr::mutate(
    block_dist = purrr::map(
      blocks,
      function(x) {
        x[is.na(x)] <- 0
        base::as.matrix(
          philentropy::distance(
            x = x,
            method = catdiss,
            mute.message = TRUE,
            p = 1
          )
        )
      }
    ),
    any_na  = purrr::map_lgl(block_dist, ~ any(is.na(.x))),
    any_inf = purrr::map_lgl(block_dist, ~ any(is.infinite(.x)))
  )

db |>
  dplyr::filter(any_na | any_inf)

sum(db$any_na)
sum(db$any_inf)
```
```{r}
purrr::map_lgl(1:nvar, function(i) {
  mats <- db |>
    dplyr::filter(row_ind == i) |>
    dplyr::pull(block_dist)

  tryCatch({
    Reduce(`+`, mats)
    TRUE
  }, error = function(e) {
    message("row_ind = ", i, " failed: ", conditionMessage(e))
    FALSE
  })
})

```
```{r}
methods_to_try <- c(
  "euclidean", "manhattan", "squared_euclidean",
  "canberra", "chebyshev", "minkowski",
  "kullback-leibler", "jensen-shannon",
  "hellinger", "bhattacharyya", "soergel"
)

method_checks <- purrr::map_dfr(methods_to_try, function(catdiss) {
  out <- tryCatch({

    db <- distance_blocks |>
      dplyr::mutate(
        block_dist = purrr::map(
          blocks,
          function(x) {
            x[is.na(x)] <- 0
            base::as.matrix(
              philentropy::distance(
                x = x,
                method = catdiss,
                mute.message = TRUE,
                p = 1
              )
            )
          }
        ),
        any_na  = purrr::map_lgl(block_dist, ~ any(is.na(.x))),
        any_inf = purrr::map_lgl(block_dist, ~ any(is.infinite(.x)))
      )

    ok_rows <- purrr::map_lgl(1:nvar, function(i) {
      mats <- db |>
        dplyr::filter(row_ind == i) |>
        dplyr::pull(block_dist)

      tryCatch({
        Reduce(`+`, mats)
        TRUE
      }, error = function(e) FALSE)
    })

    tibble::tibble(
      method = catdiss,
      ok = all(ok_rows),
      any_na = any(db$any_na),
      any_inf = any(db$any_inf),
      failed_rows = paste(which(!ok_rows), collapse = ", ")
    )

  }, error = function(e) {
    tibble::tibble(
      method = catdiss,
      ok = FALSE,
      any_na = NA,
      any_inf = NA,
      failed_rows = conditionMessage(e)
    )
  })

  out
})

method_checks |> print(n=11)
```


```{r}
x_cat <- df |>
  dplyr::select(-Name) |>
  dplyr::select(where(is.factor)) |>
  purrr::map_df(forcats::fct_drop)

Q <- purrr::map_dbl(x_cat, nlevels)
nvar <- length(Q)

z_prep <- z_preproc(x = x_cat, y = NULL, Q = Q)

ZZod <- z_prep$ZZod

level_pos <- data.table::data.table(
  start = c(1, cumsum(Q)[-length(Q)] + 1),
  stop  = cumsum(Q)
)

crs <- tidyr::crossing(a = 1:nvar, b = 1:nvar) |>
  dplyr::filter(a != b)

blocks_id_a <- level_pos[crs$a, ]
blocks_id_b <- level_pos[crs$b, ] |>
  dplyr::rename(end_start = start, end_stop = stop)

block_ids <- cbind(blocks_id_a, blocks_id_b)

pull_block <- function(start = 1, stop = 1, end_start = 1, end_stop = 1) {
  ZZod[start:stop, end_start:end_stop]
}

distance_blocks <- tibble::tibble(
  row_ind = crs$a,
  col_ind = crs$b,
  blocks = purrr::pmap(
    block_ids,
    ~ pull_block(start = ..1, stop = ..2, end_start = ..3, end_stop = ..4)
  )
)
bad_block <- distance_blocks$blocks[[15]]
bad_block
rowSums(bad_block)
conflicted::conflicts_prefer(base::as.matrix)
bad_block_fixed <- as.matrix(bad_block)
bad_block_fixed[is.na(bad_block_fixed)] <- 0
bad_block_fixed[!is.finite(bad_block_fixed)] <- 0
bad_block_fixed[bad_block_fixed < 0] <- 0

rs <- rowSums(bad_block_fixed)
zero_rows <- rs <= .Machine$double.eps
if (any(zero_rows)) {
  bad_block_fixed[zero_rows, ] <- 1 / ncol(bad_block_fixed)
  rs <- rowSums(bad_block_fixed)
}
bad_block_fixed <- sweep(bad_block_fixed, 1, rs, "/")

bad_block_fixed
rowSums(bad_block_fixed)

tmp <- philentropy::distance(
  x = bad_block_fixed,
  method = "hellinger",
  mute.message = TRUE,
  p = 1
)

tmp
class(tmp)
dim(tmp)
```
