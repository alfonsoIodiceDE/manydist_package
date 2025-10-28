mg_gower_dist_modify <- function(data.x, data.y = data.x, rngs = NULL, KR.corr = TRUE,
                              var.weights = NULL, robcb=NULL){
  
  ####################    
  source("R/mg_gower_fcn.r")
  
  ######## END  gower_fcn() ###################       
  if (is.null(dim(data.x)) && is.null(dim(data.y))) {
    out.gow <- mg_gower_fcn(x = data.x, y = data.y, rng = rngs, 
                         KR.corr = KR.corr, robcb = robcb)
    out <- (out.gow$dist * out.gow$delta)/out.gow$delta
  }
  else if (is.null(dim(data.x)) && !is.null(dim(data.y))) {
    p <- ncol(data.y)
    if (length(data.x) != p) 
      stop("data.x should be of the same length of the no. of cols of data.y")
    num <- array(0, c(1, nrow(data.y)))
    den <- array(0, c(1, nrow(data.y)))
    if(is.null(var.weights)) var.weights <- rep(1, p)
    for (k in 1:p) {
      
      if (is.null(rngs)) 
        rng.k <- NULL
      else rng.k <- rngs[k]
      w.k <- var.weights[k]
      out.gow <- mg_gower_fcn(x = data.x[, k], y = data.y[,k],
                           rng = rng.k, KR.corr = KR.corr, robcb = robcb)
      n <- out.gow$dist * out.gow$delta * w.k
      
      n[is.na(n)] <- 0
      num <- num + n
      d <- out.gow$delta * w.k
      d[is.na(d)] <- 0
      den <- den + d
    }
    out <- num/den
  }
  else {
    p <- ncol(data.y)
    if (ncol(data.x) != p) 
      stop("data.x and data.y must have the same no. of cols")
    num <- array(0, c(nrow(data.x), nrow(data.y)))
    den <- array(0, c(nrow(data.x), nrow(data.y)))
    if(is.null(var.weights)) var.weights <- rep(1, p)
    for (k in 1:p) {
      
      if (is.null(rngs)) 
        rng.k <- NULL
      else rng.k <- rngs[k]
      w.k <- var.weights[k]
      out.gow <- mg_gower_fcn(x = data.x[, k], y = data.y[, k], rng = rng.k,
                           KR.corr = KR.corr, robcb = robcb)
      n <- out.gow$dist * out.gow$delta * w.k
      
      n[is.na(n)] <- 0
      num <- num + n
      d <- out.gow$delta * w.k
      d[is.na(d)] <- 0
      den <- den + d
    }
    out <- num/den
  }
  out
}    
