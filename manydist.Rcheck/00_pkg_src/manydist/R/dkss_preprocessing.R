dkss_preprocessing<- function(df, no_f_cont){
  if (is_tibble(df)) {
    df <- as.data.frame(df)
  }
  
  # Ensure proper data types for KDML
  # Continuous variables should be numeric
  if (no_f_cont > 0 && no_f_cont <= ncol(df)) {
    for (i in 1:no_f_cont) {
      if (!is.numeric(df[, i])) {
        df[, i] <- as.numeric(as.character(df[, i]))
      }
    }
  }
  
  # Categorical variables should be factors
  if (no_f_cont < ncol(df)) {
    cat_cols <- (no_f_cont + 1):ncol(df)
    for (i in cat_cols) {
      if (!is.factor(df[, i])) {
        df[, i] <- as.factor(df[, i])
      }
    }
  }
  return(df)
}
