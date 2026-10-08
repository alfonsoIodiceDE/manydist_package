commensurability_example <- function() {
  data.frame(
    z = c(0, 2, 8, 10), w = c(1, 5, 2, 7),
    g = factor(c("a", "a", "a", "b")),
    h = factor(c("x", "x", "y", "y"))
  )
}

test_that("training means match square-matrix means without allocating a square", {
  for (x in list(c(0, 2, 8, 10), c(1, 1, 1), c(-9, 2, 2, 100), 5)) {
    expect_equal(.mean_absolute_training_distance(x), mean(abs(outer(x, x, "-"))))
  }
  Z <- cbind(c(1, 1, 1, 0), c(0, 0, 0, 1))
  delta <- matrix(c(0, 2, 2, 0), 2)
  expect_equal(.mean_categorical_training_distance(Z, delta),
               mean(Z %*% delta %*% t(Z)))
})

test_that("variable contributions use training-only means before aggregation", {
  train <- commensurability_example()[, c("z", "w", "g")]
  test <- train[c(1, 4), ]
  test$z <- c(1, 9)
  test$w <- c(3, 6)
  expected <- matrix(0, nrow(test), nrow(train))
  for (column in c("z", "w")) {
    expected <- expected + abs(outer(test[[column]], train[[column]], "-")) /
      mean(abs(outer(train[[column]], train[[column]], "-")))
  }
  expected <- expected + outer(test$g, train$g, "!=") /
    mean(outer(train$g, train$g, "!="))
  actual <- mdist(train, new_data = test, preset = "u_indep")$distance
  expect_equal(unname(as.matrix(actual)), expected)
})

test_that("commensurable presets are independent of test-batch composition", {
  train <- commensurability_example()
  test <- train[c(1, 4), ]
  test$z <- c(1, 9)
  extra <- test[1, , drop = FALSE]
  extra$z <- 100
  extra$w <- -100
  cases <- list(numerical = c("z", "w"), categorical = c("g", "h"),
                mixed = names(train))
  for (columns in cases) {
    tr <- train[, columns, drop = FALSE]
    te <- test[, columns, drop = FALSE]
    expanded <- rbind(te, extra[, columns, drop = FALSE])
    for (preset in c("u_indep", "u_dep", "u_dep_bw", "u_mix")) {
      reference <- as.matrix(mdist(tr, new_data = te, preset = preset)$distance)
      extended <- as.matrix(mdist(tr, new_data = expanded, preset = preset)$distance)
      single <- as.matrix(mdist(tr, new_data = te[1, , drop = FALSE], preset = preset)$distance)
      reversed <- as.matrix(mdist(tr, new_data = te[2:1, ], preset = preset)$distance)
      expect_equal(unname(extended[1:2, , drop = FALSE]), unname(reference))
      expect_equal(unname(single), unname(reference[1, , drop = FALSE]))
      expect_equal(unname(reversed[2:1, , drop = FALSE]), unname(reference))
      expect_equal(unname(as.matrix(mdist(tr, new_data = tr, preset = preset)$distance)),
                   unname(as.matrix(mdist(tr, preset = preset)$distance)))

      fitted <- recipes::recipe(~., data = tr) |>
        step_mdist(recipes::all_predictors(), preset = preset,
                   output = "distance_to_training") |>
        recipes::prep(training = tr)
      bake_dist <- function(data) {
        baked <- recipes::bake(fitted, new_data = data)
        as.matrix(baked[, grep("^dist_", names(baked)), drop = FALSE])
      }
      expect_equal(unname(bake_dist(te)), unname(reference))
      expect_equal(unname(bake_dist(expanded)[1:2, , drop = FALSE]), unname(reference))
      expect_equal(unname(bake_dist(te[1, , drop = FALSE])),
                   unname(reference[1, , drop = FALSE]))
    }
  }
})

test_that("custom categorical methods use training means", {
  train <- commensurability_example()[, c("g", "h")]
  test <- train[c(1, 4), ]
  expanded <- rbind(test, test[1, ])
  for (method in c("matching", "tvd", "none", "st_dev", "HL", "cat_dis",
                   "HLeucl", "mca")) {
    distances <- function(new_data) {
      as.matrix(mdist(train, new_data = new_data, preset = "custom",
                      method_cat = method, commensurable = TRUE)$distance)
    }
    expect_equal(unname(distances(expanded)[1:2, , drop = FALSE]),
                 unname(distances(test)), info = method)
    expect_equal(unname(distances(train)),
                 unname(as.matrix(mdist(train, method_cat = method)$distance)),
                 info = method)
  }
})

test_that("custom robust preprocessing is also training-only", {
  train <- commensurability_example()[, c("z", "w")]
  test <- train[c(1, 4), ]
  extra <- test[1, ]
  extra$z <- 100
  extra$w <- -100
  for (comm in c(TRUE, FALSE)) {
    distances <- function(new_data) as.matrix(mdist(
      train, new_data = new_data, method_num = "robust",
      commensurable = comm
    )$distance)
    expect_equal(unname(distances(rbind(test, extra))[1:2, ]),
                 unname(distances(test)))
    expect_equal(unname(distances(test[1, , drop = FALSE])),
                 unname(distances(test)[1, , drop = FALSE]))
  }
})

test_that("a constant unscaled numerical contribution remains finite and zero", {
  train <- data.frame(z = rep(2, 4))
  expect_equal(unname(as.matrix(mdist(train, method_num = "none")$distance)),
               matrix(0, 4, 4))
  expect_equal(unname(as.matrix(mdist(train, new_data = train[1, , drop = FALSE],
                                     method_num = "none")$distance)), matrix(0, 1, 4))
})
