library(mlbench)
library(clusterSim)
library(tibble)
n=1000
names=c('Udep Int.','Udep','naive','Gower','mod G.', 'dkss' )
pnum=4
pcat=4
p=pnum+pcat
n1=n2=500
niter=50
t_N2_NC=matrix(0,niter,6)
ARIM_N2_NC=matrix(0,niter,6)

for(ite in 1:niter){
  
  data2 <- mlbench.spirals(n, 2, 0.05)  # alternative: spirals
  d <- data2$x[order(data2$classes), ]
  moon <- shapes.two.moon(n/2)
  # pnumN=pnum-4
  # sig=matrix(0.5,pnumN,pnumN) ####correlation for the numerical
  # diag(sig)=1
  # nor1=rmvnorm(n1,rep(0,pnumN),sigma=sig) ###usi
  # sig=matrix(-0.5,pnumN,pnumN) ####correlation for the numerical
  # diag(sig)=1
  # nor2=rmvnorm(n2,rep(4,pnumN),sigma=sig) ###using mean 0
  data=cbind(d,moon$data)#,rbind(nor1,nor2))
  truth=moon$clusters
  pairs(data, col =truth, pch = 16,cex=0.8)
  ##Choice of number of binary variables 
  pbin=pcat #all categoricals binary for now
  #betaBin=matrix(c(1,rep(2,p)),pbin,(p+1),byrow = TRUE)## the beta0 is 1 all the others 2
  Xn=data[which(truth==1),]
  #XX=cbind(1,Xn) ##for intercept
  Xcat=matrix(0,n1,pcat)
  #XX=cbind(XX,Xcat)
  for(i in 1:n1){
    Xcat[i,]=rbinom(pbin,1,p=rep(0.85,pbin))
    # for(j in 1:pbin){
    #   Xcat[i,j]=rbinom(1,1,p=1/(1+exp(betaBin[j,1:(j+pbin-1)]%*%XX[i,1:(j+pbin-1)])))
    # }
  }
  #  boxplot(Xn[,1]~Xcat[,1])## just to visualize one association
  # table(Xcat[,1],Xcat[,2])## just to visualize one association
  #  table(Xcat[,1],Xcat[,3])## just to visualize one association
  X1=cbind(Xn,Xcat)
  
  ##Choice of number of binary variables 
  #betaBin=matrix(c(1,rep(-2,p)),pbin,(p+1),byrow = TRUE)## the beta0 is 1 all the others 2
  Xn=data[which(truth==2),]
  #XX=cbind(1,Xn) ##for intercept
  Xcat=matrix(0,n2,pcat)
  #XX=cbind(XX,Xcat)
  
  for(i in 1:n2){
    Xcat[i,]=rbinom(pbin,1,p=rep(0.15,pbin))
    # for(j in 1:pbin){
    #   
    #   Xcat[i,j]=rbinom(1,1,p=1/(1+exp(betaBin[j,1:(j+pnum-1)]%*%XX[i,1:(j+pnum-1)])))
    # }
  }
  
  #  boxplot(Xn[,1]~Xcat[,1])## just to visualize one association
  #  table(Xcat[,1],Xcat[,2])## just to visualize one association
  # table(Xcat[,1],Xcat[,3])## just to visualize one association
  X2=cbind(Xn,Xcat)
  
  
  X=rbind(X1,X2)
  # l=c(rep(1,n1),rep(2,n2),rep(3,n3))
  #pairs(X[,1:pnum],col=l) 
  Xkp=data.frame(X)#####convert categorical into factors
  for(i in 1:pcat){
    Xkp[,(pnum+i)]=as.factor(X[,(pnum+i)])
  }
  dfc = Xkp |> as_tibble()
  
  myk=2
  
  
  t1=proc.time()
  D_int = distance_by_method(dfc,method="u_dep",interaction=TRUE)
  t2=proc.time()
  t_N2_NC[ite,1]=t2[3]-t1[3]
  t1=proc.time()
  D_no_int = distance_by_method(dfc,method="u_dep",interaction =FALSE)
  t2=proc.time()
  t_N2_NC[ite,2]=t2[3]-t1[3]
  t1=proc.time()
  # D_gudmm = distance_by_method(dfc,method="gudmm",interaction=FALSE) |> as.matrix()
  D_naive= distance_by_method(dfc,method="naive",interaction=FALSE)
  t2=proc.time()
  t_N2_NC[ite,3]=t2[3]-t1[3]
  t1=proc.time()
  D_gow = distance_by_method(dfc,method="gower",interaction=FALSE)
  t2=proc.time()
  t_N2_NC[ite,4]=t2[3]-t1[3]
  
  t1=proc.time()
  D_mod_gow = distance_by_method(dfc,method="mod_gower",interaction=FALSE) |> as.matrix()
  t2=proc.time()
  t_N2_NC[ite,5]=t2[3]-t1[3]
  
  t1=proc.time()
  D_dkss = distance_by_method(dfc,method="dkss",interaction=FALSE)
  t2=proc.time()
  t_N2_NC[ite,6]=t2[3]-t1[3]
  
  
  
  
  cl_int = spectral_from_dist(D_int, k=myk,affinity_method = "gaussian")
  cl_no_int = spectral_from_dist(D_no_int, k=myk,affinity_method = "gaussian")
  # cl_gudmm = spectral_from_dist(D_gudmm, k=myk,affinity_method = "gaussian")
  cl_naive = spectral_from_dist(D_naive, k=myk,affinity_method = "gaussian")
  cl_gow = spectral_from_dist(D_gow, k=myk,affinity_method = "gaussian")
  cl_mod_gow = spectral_from_dist(D_mod_gow, k=myk,affinity_method = "gaussian")
  cl_dkss = spectral_from_dist(D_dkss, k=myk)
  
  
  
  ARIM_N2_NC[ite,1]=adjustedRandIndex(cl_int, truth)
  ARIM_N2_NC[ite,2]=adjustedRandIndex(cl_no_int, truth)
  #adjustedRandIndex(cl_gudmm, truth)
  ARIM_N2_NC[ite,3]=adjustedRandIndex(cl_gow, truth)
  ARIM_N2_NC[ite,4]=adjustedRandIndex(cl_mod_gow, truth)
  ARIM_N2_NC[ite,5]=adjustedRandIndex(cl_naive, truth)
  ARIM_N2_NC[ite,6]= adjustedRandIndex(cl_dkss, truth)
  
  colnames(t_N2_NC)=names
  colnames(ARIM_N2_NC)=names
  print(ite)
}

save(t_N2_NC,file='t_N2_NC2.rdata')
save(ARIM_N2_NC,file='ARIM_N2_NC2.rdata')



names=c('Udep Int.','Udep','naive','Gower','mod G.', 'dkss' )
pnum=4
pcat=2
p=pnum+pcat

niter=50
t_N2_NC_42=matrix(0,niter,6)
ARIM_N2_NC_42=matrix(0,niter,6)

for(ite in 1:niter){
  
  data2 <- mlbench.spirals(n, 2, 0.05)  # alternative: spirals
  d <- data2$x[order(data2$classes), ]
  moon <- shapes.two.moon(n/2)
  # pnumN=pnum-4
  # sig=matrix(0.5,pnumN,pnumN) ####correlation for the numerical
  # diag(sig)=1
  # nor1=rmvnorm(n1,rep(0,pnumN),sigma=sig) ###usi
  # sig=matrix(-0.5,pnumN,pnumN) ####correlation for the numerical
  # diag(sig)=1
  # nor2=rmvnorm(n2,rep(4,pnumN),sigma=sig) ###using mean 0
  data=cbind(d,moon$data)#,rbind(nor1,nor2))
  truth=moon$clusters
  pairs(data, col =truth, pch = 16,cex=0.8)
  ##Choice of number of binary variables 
  pbin=pcat #all categoricals binary for now
  #betaBin=matrix(c(1,rep(2,p)),pbin,(p+1),byrow = TRUE)## the beta0 is 1 all the others 2
  Xn=data[which(truth==1),]
  #XX=cbind(1,Xn) ##for intercept
  Xcat=matrix(0,n1,pcat)
  #XX=cbind(XX,Xcat)
  for(i in 1:n1){
    Xcat[i,]=rbinom(pbin,1,p=rep(0.85,pbin))
    # for(j in 1:pbin){
    #   Xcat[i,j]=rbinom(1,1,p=1/(1+exp(betaBin[j,1:(j+pbin-1)]%*%XX[i,1:(j+pbin-1)])))
    # }
  }
  #  boxplot(Xn[,1]~Xcat[,1])## just to visualize one association
  # table(Xcat[,1],Xcat[,2])## just to visualize one association
  #  table(Xcat[,1],Xcat[,3])## just to visualize one association
  X1=cbind(Xn,Xcat)
  
  ##Choice of number of binary variables 
  #betaBin=matrix(c(1,rep(-2,p)),pbin,(p+1),byrow = TRUE)## the beta0 is 1 all the others 2
  Xn=data[which(truth==2),]
  #XX=cbind(1,Xn) ##for intercept
  Xcat=matrix(0,n2,pcat)
  #XX=cbind(XX,Xcat)
  
  for(i in 1:n2){
    Xcat[i,]=rbinom(pbin,1,p=rep(0.15,pbin))
    # for(j in 1:pbin){
    #   
    #   Xcat[i,j]=rbinom(1,1,p=1/(1+exp(betaBin[j,1:(j+pnum-1)]%*%XX[i,1:(j+pnum-1)])))
    # }
  }
  
  #  boxplot(Xn[,1]~Xcat[,1])## just to visualize one association
  #  table(Xcat[,1],Xcat[,2])## just to visualize one association
  # table(Xcat[,1],Xcat[,3])## just to visualize one association
  X2=cbind(Xn,Xcat)
  
  
  X=rbind(X1,X2)
  # l=c(rep(1,n1),rep(2,n2),rep(3,n3))
  #pairs(X[,1:pnum],col=l) 
  Xkp=data.frame(X)#####convert categorical into factors
  for(i in 1:pcat){
    Xkp[,(pnum+i)]=as.factor(X[,(pnum+i)])
  }
  dfc = Xkp |> as_tibble()
  
  myk=2
  
  
  t1=proc.time()
  D_int = distance_by_method(dfc,method="u_dep",interaction=TRUE)
  t2=proc.time()
  t_N2_NC_42[ite,1]=t2[3]-t1[3]
  t1=proc.time()
  D_no_int = distance_by_method(dfc,method="u_dep",interaction =FALSE)
  t2=proc.time()
  t_N2_NC_42[ite,2]=t2[3]-t1[3]
  t1=proc.time()
  # D_gudmm = distance_by_method(dfc,method="gudmm",interaction=FALSE) |> as.matrix()
  D_naive= distance_by_method(dfc,method="naive",interaction=FALSE)
  t2=proc.time()
  t_N2_NC_42[ite,3]=t2[3]-t1[3]
  t1=proc.time()
  D_gow = distance_by_method(dfc,method="gower",interaction=FALSE)
  t2=proc.time()
  t_N2_NC_42[ite,4]=t2[3]-t1[3]
  
  t1=proc.time()
  D_mod_gow = distance_by_method(dfc,method="mod_gower",interaction=FALSE) |> as.matrix()
  t2=proc.time()
  t_N2_NC_42[ite,5]=t2[3]-t1[3]
  
  t1=proc.time()
  D_dkss = distance_by_method(dfc,method="dkss",interaction=FALSE)
  t2=proc.time()
  t_N2_NC_42[ite,6]=t2[3]-t1[3]
  
  
  
  
  cl_int = spectral_from_dist(D_int, k=myk,affinity_method = "gaussian")
  cl_no_int = spectral_from_dist(D_no_int, k=myk,affinity_method = "gaussian")
  # cl_gudmm = spectral_from_dist(D_gudmm, k=myk,affinity_method = "gaussian")
  cl_naive = spectral_from_dist(D_naive, k=myk,affinity_method = "gaussian")
  cl_gow = spectral_from_dist(D_gow, k=myk,affinity_method = "gaussian")
  cl_mod_gow = spectral_from_dist(D_mod_gow, k=myk,affinity_method = "gaussian")
  cl_dkss = spectral_from_dist(D_dkss, k=myk)
  
  
  
  ARIM_N2_NC_42[ite,1]=adjustedRandIndex(cl_int, truth)
  ARIM_N2_NC_42[ite,2]=adjustedRandIndex(cl_no_int, truth)
  #adjustedRandIndex(cl_gudmm, truth)
  ARIM_N2_NC_42[ite,3]=adjustedRandIndex(cl_gow, truth)
  ARIM_N2_NC_42[ite,4]=adjustedRandIndex(cl_mod_gow, truth)
  ARIM_N2_NC_42[ite,5]=adjustedRandIndex(cl_naive, truth)
  ARIM_N2_NC_42[ite,6]= adjustedRandIndex(cl_dkss, truth)
  
  colnames(t_N2_NC_42)=names
  colnames(ARIM_N2_NC_42)=names
  print(ite)
}

save(t_N2_NC_42,file='t_N2_NC_42.rdata')
save(ARIM_N2_NC_42,file='ARIM_N2_NC_42.rdata')



boxplot(ARIM_N2_NC_42, main = "", ylab='ARI',
        las = 2,                # Rotate axis labels for readability
        col = "lightblue",      # Fill color
        border = "darkblue",    # Border color
        outline = TRUE)         # Show outliers

names=c('Udep Int.','Udep','naive','Gower','mod G.', 'dkss' )
pnum=4
pcat=8
p=pnum+pcat

niter=50
t_N2_NC_48=matrix(0,niter,6)
ARIM_N2_NC_48=matrix(0,niter,6)

for(ite in 1:niter){
  
  data2 <- mlbench.spirals(n, 2, 0.05)  # alternative: spirals
  d <- data2$x[order(data2$classes), ]
  moon <- shapes.two.moon(n/2)
  # pnumN=pnum-4
  # sig=matrix(0.5,pnumN,pnumN) ####correlation for the numerical
  # diag(sig)=1
  # nor1=rmvnorm(n1,rep(0,pnumN),sigma=sig) ###usi
  # sig=matrix(-0.5,pnumN,pnumN) ####correlation for the numerical
  # diag(sig)=1
  # nor2=rmvnorm(n2,rep(4,pnumN),sigma=sig) ###using mean 0
  data=cbind(d,moon$data)#,rbind(nor1,nor2))
  truth=moon$clusters
  pairs(data, col =truth, pch = 16,cex=0.8)
  ##Choice of number of binary variables 
  pbin=pcat #all categoricals binary for now
  #betaBin=matrix(c(1,rep(2,p)),pbin,(p+1),byrow = TRUE)## the beta0 is 1 all the others 2
  Xn=data[which(truth==1),]
  #XX=cbind(1,Xn) ##for intercept
  Xcat=matrix(0,n1,pcat)
  #XX=cbind(XX,Xcat)
  for(i in 1:n1){
    Xcat[i,]=rbinom(pbin,1,p=rep(0.85,pbin))
    # for(j in 1:pbin){
    #   Xcat[i,j]=rbinom(1,1,p=1/(1+exp(betaBin[j,1:(j+pbin-1)]%*%XX[i,1:(j+pbin-1)])))
    # }
  }
  #  boxplot(Xn[,1]~Xcat[,1])## just to visualize one association
  # table(Xcat[,1],Xcat[,2])## just to visualize one association
  #  table(Xcat[,1],Xcat[,3])## just to visualize one association
  X1=cbind(Xn,Xcat)
  
  ##Choice of number of binary variables 
  #betaBin=matrix(c(1,rep(-2,p)),pbin,(p+1),byrow = TRUE)## the beta0 is 1 all the others 2
  Xn=data[which(truth==2),]
  #XX=cbind(1,Xn) ##for intercept
  Xcat=matrix(0,n2,pcat)
  #XX=cbind(XX,Xcat)
  
  for(i in 1:n2){
    Xcat[i,]=rbinom(pbin,1,p=rep(0.15,pbin))
    # for(j in 1:pbin){
    #   
    #   Xcat[i,j]=rbinom(1,1,p=1/(1+exp(betaBin[j,1:(j+pnum-1)]%*%XX[i,1:(j+pnum-1)])))
    # }
  }
  
  #  boxplot(Xn[,1]~Xcat[,1])## just to visualize one association
  #  table(Xcat[,1],Xcat[,2])## just to visualize one association
  # table(Xcat[,1],Xcat[,3])## just to visualize one association
  X2=cbind(Xn,Xcat)
  
  
  X=rbind(X1,X2)
  # l=c(rep(1,n1),rep(2,n2),rep(3,n3))
  #pairs(X[,1:pnum],col=l) 
  Xkp=data.frame(X)#####convert categorical into factors
  for(i in 1:pcat){
    Xkp[,(pnum+i)]=as.factor(X[,(pnum+i)])
  }
  dfc = Xkp |> as_tibble()
  
  myk=2
  
  
  t1=proc.time()
  D_int = distance_by_method(dfc,method="u_dep",interaction=TRUE)
  t2=proc.time()
  t_N2_NC_48[ite,1]=t2[3]-t1[3]
  t1=proc.time()
  D_no_int = distance_by_method(dfc,method="u_dep",interaction =FALSE)
  t2=proc.time()
  t_N2_NC_48[ite,2]=t2[3]-t1[3]
  t1=proc.time()
  # D_gudmm = distance_by_method(dfc,method="gudmm",interaction=FALSE) |> as.matrix()
  D_naive= distance_by_method(dfc,method="naive",interaction=FALSE)
  t2=proc.time()
  t_N2_NC_48[ite,3]=t2[3]-t1[3]
  t1=proc.time()
  D_gow = distance_by_method(dfc,method="gower",interaction=FALSE)
  t2=proc.time()
  t_N2_NC_48[ite,4]=t2[3]-t1[3]
  
  t1=proc.time()
  D_mod_gow = distance_by_method(dfc,method="mod_gower",interaction=FALSE) |> as.matrix()
  t2=proc.time()
  t_N2_NC_48[ite,5]=t2[3]-t1[3]
  
  t1=proc.time()
  D_dkss = distance_by_method(dfc,method="dkss",interaction=FALSE)
  t2=proc.time()
  t_N2_NC_48[ite,6]=t2[3]-t1[3]
  
  
  
  
  cl_int = spectral_from_dist(D_int, k=myk,affinity_method = "gaussian")
  cl_no_int = spectral_from_dist(D_no_int, k=myk,affinity_method = "gaussian")
  # cl_gudmm = spectral_from_dist(D_gudmm, k=myk,affinity_method = "gaussian")
  cl_naive = spectral_from_dist(D_naive, k=myk,affinity_method = "gaussian")
  cl_gow = spectral_from_dist(D_gow, k=myk,affinity_method = "gaussian")
  cl_mod_gow = spectral_from_dist(D_mod_gow, k=myk,affinity_method = "gaussian")
  cl_dkss = spectral_from_dist(D_dkss, k=myk)
  
  
  
  ARIM_N2_NC_48[ite,1]=adjustedRandIndex(cl_int, truth)
  ARIM_N2_NC_48[ite,2]=adjustedRandIndex(cl_no_int, truth)
  #adjustedRandIndex(cl_gudmm, truth)
  ARIM_N2_NC_48[ite,3]=adjustedRandIndex(cl_gow, truth)
  ARIM_N2_NC_48[ite,4]=adjustedRandIndex(cl_mod_gow, truth)
  ARIM_N2_NC_48[ite,5]=adjustedRandIndex(cl_naive, truth)
  ARIM_N2_NC_48[ite,6]= adjustedRandIndex(cl_dkss, truth)
  
  colnames(t_N2_NC_48)=names
  colnames(ARIM_N2_NC_48)=names
  print(ite)
}

save(t_N2_NC_48,file='t_N2_NC_48.rdata')
save(ARIM_N2_NC_48,file='ARIM_N2_NC_48.rdata')



boxplot(ARIM_N2_NC_48, main = "", ylab='ARI',
        las = 2,                # Rotate axis labels for readability
        col = "lightblue",      # Fill color
        border = "darkblue",    # Border color
        outline = TRUE)         # Show outliers







