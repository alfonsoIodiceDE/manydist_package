library(doParallel)
parallel::detectCores()                        #How many cores are available to the computer?
n.cores = 8       
my.cluster = parallel::makeCluster(            #Set up the parallel backend
  n.cores, 
  type = "PSOCK"
)
print(my.cluster)                               #Check if its successful
doParallel::registerDoParallel(cl = my.cluster) #Register it to execute loops in parallel.
foreach::getDoParRegistered()                   #Check if the registration is successful.
foreach::getDoParWorkers()                      #Check how many workers are available.


### Generate data from normal plus glm 3
library(mvtnorm)
library(Rmixmod)
names=c('Udep Int.','Udep','naive','Gower','mod G.', 'dkss' )
pnum=15
pcat=15
p=pnum+pcat
n=500

niter=50
t_N1E=matrix(0,niter,6)
ARIM_N1E=matrix(0,niter,6)

for(ite in 1:niter){
  ####### Cluster 1 
  n1=round(0.2*n)
  sig=matrix(0.5,pnum,pnum) ####correlation for the numerical
  diag(sig)=1
  Xn=rmvnorm(n1,rep(0,pnum),sigma=sig) ###using mean 0
  
  ##Choice of number of binary variables 
  pbin=pcat #all categoricals binary for now
  betaBin=matrix(c(1,rep(2,p)),pbin,(p+1),byrow = TRUE)## the beta0 is 1 all the others 2
  
  XX=cbind(1,Xn) ##for intercept
  Xcat=matrix(0,n1,pcat)
  XX=cbind(XX,Xcat)
  for(i in 1:n1){
    for(j in 1:pbin){
      Xcat[i,j]=rbinom(1,1,p=1/(1+exp(betaBin[j,1:(j+pnum-1)]%*%XX[i,1:(j+pnum-1)])))
    }}
  
#  boxplot(Xn[,1]~Xcat[,1])## just to visualize one association
 # table(Xcat[,1],Xcat[,2])## just to visualize one association
#  table(Xcat[,1],Xcat[,3])## just to visualize one association
  X1=cbind(Xn,Xcat)
  
  ####### Cluster 2 
  n2=round(0.3*n)
  sig=matrix(0.8,pnum,pnum) ####correlation for the numerical
  diag(sig)=1
  Xn=rmvnorm(n2,rep(0,pnum),sigma=sig) ###using mean 0
  
  ##Choice of number of binary variables 
  betaBin=matrix(c(1,rep(-2,p)),pbin,(p+1),byrow = TRUE)## the beta0 is 1 all the others 2
  
  XX=cbind(1,Xn) ##for intercept
  Xcat=matrix(0,n2,pcat)
  XX=cbind(XX,Xcat)

    for(i in 1:n2){
      for(j in 1:pbin){
        Xcat[i,j]=rbinom(1,1,p=1/(1+exp(betaBin[j,1:(j+pnum-1)]%*%XX[i,1:(j+pnum-1)])))
      }}

#  boxplot(Xn[,1]~Xcat[,1])## just to visualize one association
#  table(Xcat[,1],Xcat[,2])## just to visualize one association
 # table(Xcat[,1],Xcat[,3])## just to visualize one association
  X2=cbind(Xn+4,Xcat)
  
  ####### Cluster 3 
  n3=round(n*0.5)
  sig=matrix(-0.5,pnum,pnum) ####correlation for the numerical
  diag(sig)=1
  Xn=rmvnorm(n3,rep(0,pnum),sigma=sig) ###using mean 0
  
  ##Choice of number of binary variables 
  
  betaBin=matrix(c(1,rep(2,p)),pbin,(p+1),byrow = TRUE)## the beta0 is 1 all the others 2
  
  XX=cbind(1,Xn) ##for intercept
  Xcat=matrix(0,n3,pcat)
  XX=cbind(XX,Xcat)
  for(i in 1:n3){
    for(j in 1:pbin){
      Xcat[i,j]=rbinom(1,1,p=1/(1+exp(betaBin[j,1:(j+pnum-1)]%*%XX[i,1:(j+pnum-1)])))
    }}
  
# boxplot(Xn[,1]~Xcat[,1])## just to visualize one association
 # table(Xcat[,1],Xcat[,2])## just to visualize one association
  #table(Xcat[,1],Xcat[,3])## just to visualize one association
  Xn=sweep(Xn,2,rep(c(-1,5),length=pnum),FUN="+")
  X3=cbind(Xn,Xcat)
  
  X=rbind(X1,X2,X3)
 # l=c(rep(1,n1),rep(2,n2),rep(3,n3))
  #pairs(X[,1:pnum],col=l) 
  Xkp=data.frame(X)#####convert categorical into factors
  for(i in 1:pcat){
    Xkp[,(pnum+i)]=as.factor(X[,(pnum+i)])
  }
  dfc = Xkp |> as_tibble()
  truth =c(rep(1,n1),rep(2,n2),rep(3,n3))
  myk=3
  
  t1=proc.time()
  D_int = distance_by_method(dfc,method="u_dep",interaction=TRUE)
  t2=proc.time()
  t_N1E[ite,1]=t2[3]-t1[3]
  t1=proc.time()
  D_no_int = distance_by_method(dfc,method="u_dep",interaction =FALSE)
  t2=proc.time()
  t_N1E[ite,2]=t2[3]-t1[3]
  t1=proc.time()
  # D_gudmm = distance_by_method(dfc,method="gudmm",interaction=FALSE) |> as.matrix()
  D_naive= distance_by_method(dfc,method="naive",interaction=FALSE)
  t2=proc.time()
  t_N1E[ite,3]=t2[3]-t1[3]
  t1=proc.time()
  D_gow = distance_by_method(dfc,method="gower",interaction=FALSE)
  t2=proc.time()
  t_N1E[ite,4]=t2[3]-t1[3]
  
  t1=proc.time()
  D_mod_gow = distance_by_method(dfc,method="mod_gower",interaction=FALSE) |> as.matrix()
  t2=proc.time()
  t_N1E[ite,5]=t2[3]-t1[3]
  
  t1=proc.time()
  D_dkss = distance_by_method(dfc,method="dkss",interaction=FALSE)
  t2=proc.time()
  t_N1E[ite,6]=t2[3]-t1[3]
  
  
  
  
  cl_int = spectral_from_dist(D_int, k=myk,affinity_method = "gaussian")
  cl_no_int = spectral_from_dist(D_no_int, k=myk,affinity_method = "gaussian")
  # cl_gudmm = spectral_from_dist(D_gudmm, k=myk,affinity_method = "gaussian")
  cl_naive = spectral_from_dist(D_naive, k=myk,affinity_method = "gaussian")
  cl_gow = spectral_from_dist(D_gow, k=myk,affinity_method = "gaussian")
  cl_mod_gow = spectral_from_dist(D_mod_gow, k=myk,affinity_method = "gaussian")
  cl_dkss = spectral_from_dist(D_dkss, k=myk)
  
  
  
  ARIM_N1E[ite,1]=adjustedRandIndex(cl_int, truth)
  ARIM_N1E[ite,2]=adjustedRandIndex(cl_no_int, truth)
  #adjustedRandIndex(cl_gudmm, truth)
  ARIM_N1E[ite,3]=adjustedRandIndex(cl_gow, truth)
  ARIM_N1E[ite,4]=adjustedRandIndex(cl_mod_gow, truth)
  ARIM_N1E[ite,5]=adjustedRandIndex(cl_naive, truth)
  ARIM_N1E[ite,6]= adjustedRandIndex(cl_dkss, truth)
  
colnames(t_N1E)=names
colnames(ARIM_N1E)=names
print(ite)
}

save(t_N1E,file='t_N1E.rdata')
save(ARIM_N1E,file='ARIM_N1E.rdata')



boxplot(ARIM_N1E, main = "", ylab='ARI',
            las = 2,                # Rotate axis labels for readability
               col = "lightblue",      # Fill color
                border = "darkblue",    # Border color
                 outline = TRUE)         # Show outliers


