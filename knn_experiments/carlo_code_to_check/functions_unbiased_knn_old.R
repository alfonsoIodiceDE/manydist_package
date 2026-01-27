# to compute knn predictions (with distance matrix as input)
knn_dist <- function(dist_matrix, train_labels, k) {
  n_test <- ncol(dist_matrix)
  predictions <- vector("character", n_test)
  for (i in 1:n_test) {
    dist_col <- dist_matrix[, i]
    nn_indices <- order(dist_col)[1:k]
    nn_labels <- train_labels[nn_indices]
    predictions[i] <- names(sort(table(nn_labels), decreasing = TRUE))[1]
  }
  return(predictions)
}

# to compute the accuracy
compute_accuracy <- function(predicted, actual, show_confusion = FALSE) {
  predicted <- factor(predicted)
  actual <- factor(actual, levels = levels(predicted))
  acc <- mean(predicted == actual)
  if (show_confusion) {
    print(table(Predicted = predicted, Actual = actual))
  }
  return(acc)
}
compute_all_dist_train_test <- function(tr_df, ts_df) {
  list(
    g    = t(mdist(tr_df, new_data  = ts_df, preset = "gower")$distance |> as.matrix()),
    oh   = t(mdist(tr_df, new_data = ts_df, preset = "euclidean_onehot")$distance |> as.matrix()),
    hl   = t(mdist(tr_df, new_data  = ts_df, preset = "custom",
                   distance_cont = "euclidean", distance_cat = "HLeucl",
                   scaling_cont  = "std")$distance |> as.matrix()),
    hla  = t(mdist(tr_df, new_data = ts_df, preset = "custom",
                   distance_cont = "manhattan", distance_cat = "HL",
                   scaling_cont  = "std")$distance |> as.matrix()),
    uind = t(mdist(tr_df, new_data = ts_df,
                   distance_cont = "manhattan", distance_cat = "cat_dis",
                   commensurable = TRUE)$distance |> as.matrix()),
    udep = t(mdist(tr_df, new_data = ts_df, preset = "unbiased_dependent")$distance |> as.matrix())
  )
}

compute_all_dist_full <- function(X_full) {
  list(
    g    = as.matrix(mdist(X_full, preset = "gower")$distance |> as.matrix()),
    oh   = as.matrix(mdist(X_full, preset = "euclidean_onehot")$distance |> as.matrix()),
    hl   = as.matrix(mdist(X_full, preset = "custom",
                           distance_cont = "euclidean",
                           distance_cat  = "HLeucl",
                           scaling_cont  = "std")$distance |> as.matrix()),
    hla  = as.matrix(mdist(X_full, preset = "custom",
                           distance_cont = "manhattan",
                           distance_cat  = "HL",
                           scaling_cont  = "std")$distance |> as.matrix()),
    uind = as.matrix(mdist(X_full,
                           distance_cont  = "manhattan",
                           distance_cat   = "cat_dis",
                           commensurable  = TRUE)$distance |> as.matrix()),
    udep = as.matrix(mdist(X_full, preset = "unbiased_dependent")$distance |> as.matrix())
  )
}
