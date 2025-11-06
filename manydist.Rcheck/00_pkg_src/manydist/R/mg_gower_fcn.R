mg_gower_fcn <- function(x, y, rng = NULL, KR.corr = TRUE, robcb = "iqr") {
  nx <- length(x)
  ny <- length(y)
  cx <- class(x)
  cy <- class(y)
  delta <- matrix(1, nx, ny)
  if (!identical(cx, cy)) 
    stop("the x and y object are of different type")
  if (is.logical(x)) {
    dd <- abs(outer(X = x, Y = y, FUN = "-"))
    delta[outer(x == FALSE, y == FALSE, FUN = "&")] <- 0
    delta[outer(is.na(x), is.na(y), FUN = "|")] <- 0
  }
  else if (is.character(x) || (is.factor(x) && !is.ordered(x))) {
    if (is.factor(x) && !identical(levels(x), levels(y))) 
      stop("x and y have different levels")
    
    # Create a unified set of levels
    unified_levels <- union(levels(factor(x)), levels(factor(y)))
    
    # Create dummy variables using the unified levels
    x_dummies <- model.matrix(~ factor(x, levels = unified_levels) - 1)
    y_dummies <- model.matrix(~ factor(y, levels = unified_levels) - 1)
    
    
    # Initialize distance matrix for dummy variables
    dd_dummies <- matrix(0, nrow(x_dummies), nrow(y_dummies))
    
    # Calculate distance for each dummy variable
    for (i in 1:ncol(x_dummies)) {
      #p <- mean(x_dummies[, i], na.rm = TRUE)
      #std_dev <- sqrt(p * (1 - p))
      scaled_x <- x_dummies[, i] 
      scaled_y <- y_dummies[, i] 
      dd_dummies <- dd_dummies + abs(outer(scaled_x, scaled_y, FUN = "-"))
    }
    
    # Scale the summed distances
    # Each variable will be x, and calculate with one of the rest y
    dd <- dd_dummies / ncol(x_dummies)
    # Apply a softening function like sigmoid
    # sigmoid <- function(x) 1 / (1 + exp(-x))
    # dd <- sigmoid(dd_dummies)
    delta[outer(is.na(x), is.na(y), FUN = "|")] <- 0
  }
  else if (is.ordered(x)) {
    if (KR.corr) {
      x <- as.numeric(x)
      y <- as.numeric(y)
      if (is.null(rng) || is.na(rng)) 
        rng <- max(x, y, na.rm = TRUE) - 1
      if(rng==0) {
        dd <- matrix(0, nx, ny)
        delta[outer(is.na(x), is.na(y), FUN = "|")] <- 0
      }    
      else{
        zx <- (x - 1)/rng
        zy <- (y - 1)/rng
        dd <- abs(outer(X = zx, Y = zy, FUN = "-"))/(max(zx, 
                                                         zy, na.rm=TRUE) - min(zx, zy, na.rm=TRUE))
        delta[outer(is.na(zx), is.na(zy), FUN = "|")] <- 0
      }
    }
    else {
      x <- as.numeric(x)
      y <- as.numeric(y)
      if (is.null(rng) || is.na(rng)) 
        rng <- max(x, y, na.rm=TRUE) - 1
      if(rng==0) dd <- matrix(0, nx, ny)
      else dd <- abs(outer(X = x, Y = y, FUN = "-"))/rng
      delta[outer(is.na(x), is.na(y), FUN = "|")] <- 0
    }
  }
  else {
    if (is.null(rng) || is.na(rng)) rng <- max(x, y, na.rm=TRUE) - min(x, y, na.rm=TRUE)
    if(!is.null(robcb)){
      if(tolower(robcb)=="iqr") {
        rng <- IQR(x=c(x,y), na.rm = TRUE)
      }
      if(tolower(robcb)=="idr") {
        rng <- c(quantile(x = c(x,y), probs=0.9, na.rm=TRUE) - 
                   quantile(x = c(x,y), probs=0.1, na.rm=TRUE))
      }
    }  
    
    if(rng==0) dd <- matrix(0, nx, ny)
    else dd <- abs(outer(X = x, Y = y, FUN = "-"))/rng
    dd[dd>1] <- 1
    
    delta[outer(is.na(x), is.na(y), FUN = "|")] <- 0
  }
  list(dist = dd, delta = delta)
}
