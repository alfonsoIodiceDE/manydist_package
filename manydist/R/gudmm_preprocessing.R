gudmm_preprocessing <- function(df, no_f_cont){
  # Convert tibble to data.frame if needed
  if (is_tibble(df)) {
    df <- as.data.frame(df)
  }
  # Ensure no_f_cont is provided
  
  
  
  X_data <- df
  X_processed <- X_data
  
  # Handle missing no_f_cont (set default to 0)
  if (is.na(no_f_cont) || is.null(no_f_cont)) {
    no_f_cont <- 0
  }
  
  # Convert continuous features to numeric
  if (no_f_cont > 0 && no_f_cont <= ncol(X_processed)) {
    for (i in 1:no_f_cont) {
      col_data <- X_processed[, i]
      
      if (is.factor(col_data)) {
        # Convert factor levels to numeric if they represent numbers
        char_vals <- as.character(col_data)
        numeric_vals <- suppressWarnings(as.numeric(char_vals))
        
        # If conversion creates all NAs, use factor indices instead
        if (all(is.na(numeric_vals))) {
          X_processed[, i] <- as.numeric(col_data)
        } else {
          X_processed[, i] <- numeric_vals
        }
      } else if (is.character(col_data)) {
        numeric_vals <- suppressWarnings(as.numeric(col_data))
        
        # If conversion fails, convert to factor first then to numeric
        if (all(is.na(numeric_vals))) {
          X_processed[, i] <- as.numeric(as.factor(col_data))
        } else {
          X_processed[, i] <- numeric_vals
        }
      } else {
        X_processed[, i] <- as.numeric(col_data)
      }
    }
  }
  
  # Convert categorical features to integers
  if (no_f_cont < ncol(X_processed)) {
    cat_cols <- (no_f_cont + 1):ncol(X_processed)
    
    for (i in cat_cols) {
      col_data <- X_processed[, i]
      
      if (is.factor(col_data)) {
        X_processed[, i] <- as.integer(col_data)
      } else if (is.character(col_data)) {
        X_processed[, i] <- as.integer(as.factor(col_data))
      } else if (is.logical(col_data)) {
        X_processed[, i] <- as.integer(col_data) + 1  # Convert TRUE/FALSE to 2/1
      } else {
        # Force to integer, ensuring positive values
        int_vals <- as.integer(col_data)
        min_val <- min(int_vals, na.rm = TRUE)
        if (!is.na(min_val) && min_val <= 0) {
          int_vals <- int_vals - min_val + 1
        }
        X_processed[, i] <- int_vals
      }
    }
  }
  
  # Convert to matrix - ensure all columns are numeric
  X_matrix <- as.matrix(X_processed)
  
  # Ensure matrix is numeric
  if (storage.mode(X_matrix) != "double") {
    storage.mode(X_matrix) <- "double"
  }
  
  # Scale continuous features (min-max normalization)
  if (no_f_cont > 0 && no_f_cont <= ncol(X_matrix)) {
    for (i in 1:no_f_cont) {
      col_vals <- X_matrix[, i]
      
      # Replace infinite values with NA
      col_vals[is.infinite(col_vals)] <- NA
      col_vals_clean <- col_vals[!is.na(col_vals)]
      
      if (length(col_vals_clean) > 1) {
        min_val <- min(col_vals_clean)
        max_val <- max(col_vals_clean)
        
        if (max_val > min_val) {
          scaled_vals <- (col_vals_clean - min_val) / (max_val - min_val)
          X_matrix[!is.na(col_vals), i] <- scaled_vals
        } else {
          # All values are the same
          X_matrix[!is.na(col_vals), i] <- 0
        }
      } else if (length(col_vals_clean) == 1) {
        # Only one non-NA value
        X_matrix[!is.na(col_vals), i] <- 0
      }
    }
  }
  
  # Ensure categorical features are positive integers
  if (no_f_cont < ncol(X_matrix)) {
    cat_cols <- (no_f_cont + 1):ncol(X_matrix)
    
    for (i in cat_cols) {
      col_vals <- X_matrix[, i]
      col_vals_clean <- col_vals[!is.na(col_vals)]
      
      if (length(col_vals_clean) > 0) {
        min_val <- min(col_vals_clean)
        if (min_val <= 0) {
          adjustment <- 1 - min_val
          X_matrix[!is.na(col_vals), i] <- col_vals_clean + adjustment
        }
        
        # Ensure they are integers
        X_matrix[!is.na(col_vals), i] <- round(X_matrix[!is.na(col_vals), i])
      }
    }
  }
  
  # Remove rows that are all NA
  complete_rows <- apply(X_matrix, 1, function(x) !all(is.na(x)))
  
  if (sum(complete_rows) < 3) {
    stop("Too few complete rows after preprocessing")
  }
  
  if (sum(complete_rows) < nrow(X_matrix)) {
    X_matrix <- X_matrix[complete_rows, , drop = FALSE]
  }
  
  # Final validation
  if (any(is.na(X_matrix))) {
    # Replace remaining NAs with appropriate values
    for (j in 1:ncol(X_matrix)) {
      na_indices <- is.na(X_matrix[, j])
      if (any(na_indices)) {
        if (j <= no_f_cont) {
          # For continuous: use median
          replacement <- median(X_matrix[!na_indices, j], na.rm = TRUE)
          if (is.na(replacement)) replacement <- 0
        } else {
          # For categorical: use mode (most frequent value)
          tab <- table(X_matrix[!na_indices, j])
          replacement <- as.numeric(names(tab)[which.max(tab)])
          if (is.na(replacement)) replacement <- 1
        }
        X_matrix[na_indices, j] <- replacement
      }
    }
  }
  return(X_matrix)
}

