cross_distance_example <- function() {
  i <- seq_len(36)
  data.frame(
    z = i + sin(i), w = cos(i / 2) * 7 + i / 3,
    g = factor(rep(c("a", "b", "c"), 12)),
    h = factor(rep(c("x", "x", "y", "z"), 9)),
    j = factor(rep(c("u", "v"), 18))
  )
}

test_that("Euclidean cross-distances match independent rowwise calculations", {
  for (n in c(1, 3, 6, 36)) {
    train <- cbind(seq_len(n), sin(seq_len(n)))
    test <- rbind(c(-3, 7), c(100, -2), train[1, ])
    expected <- matrix(0, nrow(test), n)
    for (i in seq_len(nrow(test))) {
      expected[i, ] <- sqrt(rowSums(sweep(train, 2, test[i, ], "-")^2))
    }
    expect_equal(.cross_euclidean_distance(test, train), expected)
    expect_equal(.cross_euclidean_distance(train, train),
                 unname(as.matrix(stats::dist(train))))
    expect_equal(.cross_euclidean_distance(test[1, , drop = FALSE], train),
                 expected[1, , drop = FALSE])
  }
  train <- matrix(c(1e10, 1e10 + 1), ncol = 1)
  expect_equal(.cross_euclidean_distance(train, train), matrix(c(0, 1, 1, 0), 2))
  expect_equal(dim(.cross_euclidean_distance(matrix(numeric(), 0, 1), train)),
               c(0L, 2L))
})

test_that("non-commensurable numerical operators are fitted on training data", {
  train <- cross_distance_example()[, c("z", "w")]
  test <- train[c(3, 8, 17), ]
  test$z <- c(-50, 11, 100)
  test$w <- c(2, -30, 75)
  for (method in c("none", "std", "range", "robust", "pc_scores")) {
    tr <- as.matrix(train)
    te <- as.matrix(test)
    if (method == "std") {
      scales <- apply(tr, 2, stats::sd)
      tr <- sweep(tr, 2, scales, "/")
      te <- sweep(te, 2, scales, "/")
    } else if (method == "robust") {
      centers <- apply(tr, 2, stats::median)
      scales <- apply(tr, 2, stats::IQR)
      tr <- sweep(sweep(tr, 2, centers, "-"), 2, scales, "/")
      te <- sweep(sweep(te, 2, centers, "-"), 2, scales, "/")
    } else if (method == "range") {
      centers <- apply(tr, 2, min)
      scales <- apply(tr, 2, max) - centers
      tr <- sweep(sweep(tr, 2, centers, "-"), 2, scales, "/")
      te <- sweep(sweep(te, 2, centers, "-"), 2, scales, "/")
      # Match the existing step_range() policy for out-of-range observations.
      te[] <- pmin(1, pmax(0, te))
    } else if (method == "pc_scores") {
      centers <- colMeans(tr)
      scales <- apply(tr, 2, stats::sd)
      tr <- sweep(sweep(tr, 2, centers, "-"), 2, scales, "/")
      te <- sweep(sweep(te, 2, centers, "-"), 2, scales, "/")
      rotation <- stats::prcomp(tr, center = FALSE, scale. = FALSE)$rotation
      tr <- tr %*% rotation
      te <- te %*% rotation
    }
    expected <- Reduce(`+`, lapply(seq_len(ncol(tr)), function(column) {
      abs(outer(te[, column], tr[, column], "-"))
    }))
    observed <- mdist(train, new_data = test, method_num = method,
                      commensurable = FALSE)
    expect_equal(unname(as.matrix(observed$distance)), unname(expected), info = method)
    if (method == "pc_scores") {
      reduced <- mdist(train, new_data = test, method_num = method,
                       ncomp = 1, commensurable = FALSE)
      expect_equal(unname(as.matrix(reduced$distance)),
                   unname(abs(outer(te[, 1], tr[, 1], "-"))))
    }
  }
})

test_that("HLeucl preserves its trained representation on larger datasets", {
  train <- cross_distance_example()[, c("g", "h", "j")]
  test <- train[c(3, 8, 17), ]
  rec <- recipes::recipe(~., data = train) |>
    recipes::step_dummy(recipes::all_nominal_predictors(), one_hot = TRUE) |>
    recipes::prep(training = train)
  eta <- unlist(lapply(train, function(x) {
    rep(fpc::distancefactor(cat = nlevels(x), catsizes = table(x)), nlevels(x))
  }))
  Z <- as.matrix(recipes::bake(rec, new_data = NULL)) %*% diag(eta)
  expected_training <- as.matrix(stats::dist(Z))
  expected_test <- expected_training[c(3, 8, 17), , drop = FALSE]
  for (comm in c(FALSE, TRUE)) {
    distance <- function(data) as.matrix(mdist(
      train, new_data = data, method_cat = "HLeucl", commensurable = comm
    )$distance)
    original <- as.matrix(mdist(train, method_cat = "HLeucl",
                                commensurable = comm)$distance)
    observed <- distance(test)
    expect_gt(max(observed), 0)
    expect_equal(unname(distance(train)), unname(original))
    expect_equal(unname(observed), unname(original[c(3, 8, 17), ]))
    expect_equal(unname(distance(rbind(test, train[1, ]))[1:3, ]), unname(observed))
    expect_equal(unname(distance(test[1, , drop = FALSE])),
                 unname(observed[1, , drop = FALSE]))
    if (!comm) expect_equal(unname(observed), unname(expected_test))

    fitted <- recipes::recipe(~., data = train) |>
      step_mdist(recipes::all_predictors(), method_cat = "HLeucl",
                 commensurable = comm) |>
      recipes::prep(training = train)
    baked <- recipes::bake(fitted, new_data = test)
    expect_equal(unname(as.matrix(baked[, grep("^dist_", names(baked))])),
                 unname(observed))
  }
})

test_that("Gower and hl reproduce training distances through new_data", {
  train <- cross_distance_example()
  cases <- list(numerical = c("z", "w"), categorical = c("g", "h", "j"),
                mixed = names(train))
  for (columns in cases) {
    tr <- train[, columns, drop = FALSE]
    test <- tr[c(3, 8, 17), ]
    for (preset in c("gower", "hl")) {
      original <- as.matrix(mdist(tr, preset = preset)$distance)
      distance <- function(data) as.matrix(mdist(tr, new_data = data, preset = preset)$distance)
      observed <- distance(test)
      expect_equal(unname(distance(tr)), unname(original))
      expect_equal(unname(observed), unname(original[c(3, 8, 17), ]))
      expect_equal(unname(distance(rbind(test, tr[1, ]))[1:3, ]), unname(observed))
      expect_equal(unname(distance(test[1, , drop = FALSE])),
                   unname(observed[1, , drop = FALSE]))
      fitted <- recipes::recipe(~., data = tr) |>
        step_mdist(recipes::all_predictors(), preset = preset) |>
        recipes::prep(training = tr)
      baked <- recipes::bake(fitted, new_data = test)
      expect_equal(unname(as.matrix(baked[, grep("^dist_", names(baked))])),
                   unname(observed))
      if (preset == "gower") {
        expect_equal(unname(observed),
                     unname(as.matrix(cluster::daisy(tr, metric = "gower"))[c(3, 8, 17), ]))
        summed <- mdist(tr, new_data = test, preset = "gower", gower_average = FALSE)
        expect_equal(unname(as.matrix(summed$distance)), unname(observed * ncol(tr)))
      }
    }
  }
})
