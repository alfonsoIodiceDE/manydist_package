congruence_coeff <- function(L1,L2){
  L1<-dist(L1) |> as.matrix()
  L2<-dist(L2)|> as.matrix()
    CC<-sum(diag(t(L1) %*% L2)) / (sqrt(sum(diag(t(L1) %*% L1))*sum(diag(t(L2) %*% L2))))
    return(CC)
  }
