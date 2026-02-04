library(tidyverse)
library(sampling)
library(SimDesign)
library(mclust)
conflicted::conflicts_prefer(SimDesign::rmvnorm)
devtools::load_all()
maxite=2
ARIM_x_N2=matrix(0,maxite,5)
t_x_N2=matrix(0,maxite,5)
n=1000
pnum=6
for(ite in 1:maxite){

  c1=rmvnorm(n,rep(0,pnum))

  cat1=rep(1,n)
  sel=which(srswor(n*0.6,n)==1)
  cat1[sel[1:n*0.2]]=2
  cat1[sel[(n*0.2+1):n*0.6]]=3
  cat1=factor(cat1)

  cat2=rep(1,n)
  cat2[which(srswor(n*0.2,n)==1)]=2
  cat2=factor(cat2)

  cat4=rep(1,n)
  sel=which(srswor(n*0.6,n)==1)
  cat4[sel[1:n*0.3]]=2
  cat4[sel[(n*0.3+1):n*0.6]]=3
  cat4=factor(cat4)

  cat3=rep(1,n)
  cat3[which(srswor(n*0.4,n)==1)]=2
  cat3=factor(cat3)


  pcat=4

  beta=matrix(c(0.8,0.5,-0.4,0.8,-0.4,0.5),3,2,1)
  truth=rep(3,nrow(c1))
  for(i in 1:3){
    for(j in 1:2){
      for(h in (pnum/2):1){
        c1[which(cat1==i&cat2==j),h]=apply(c1[which(cat1==i&cat2==j),(h+1):pnum]*beta[i,j],1,sum) #+rnorm(length(which(cat1==i&cat2==j)))*0.02
        if(beta[i,j]==0.5){c1[which(cat1==i&cat2==j),h]=c1[which(cat1==i&cat2==j),h]+8}
      }
      if(beta[i,j]==0.8){truth[which(cat1==i&cat2==j)]=1
      }else if(beta[i,j]==-0.4){truth[which(cat1==i&cat2==j)]=2
      }}
  }

  test=cbind(c1,cat1,cat2,cat3,cat4)

  Xkp=data.frame(test)
  nn=paste('V',1:pnum)
  nncat=paste('C',1:pcat)
  names(Xkp)=c(nn,nncat)
  for(i in 1:pnum){
    Xkp[,i]=as.vector( Xkp[,i])

  }



  #Xkp[which(Xkp[,15]>0),15]=2
  #Xkp[which(Xkp[,15]<0),15]=1


  for(i in (pnum+1):(pnum+pcat)){
    Xkp[,i]=factor(unlist(Xkp[,i]))
  }


  #dd=mice(Xkp,m=10)
  #Comp=complete(dd)


  maxite=1
  #ARIM_Real=matrix(0,maxite,8)
  #for(ite in 1:maxite){

  # for(ite in 1:10){
  # Xkp=complete(dd,ite)
  myk=3
  dfc = Xkp |> as_tibble()



  t1=proc.time()
  udep_int = mdist(x=dfc,preset="custom", distance_cont = "euclidean",
                   distance_cat = "tot_var_dist",
                   commensurable = FALSE,
                   scaling_cont = 'pc_scores',interaction = TRUE,prop_nn = .05,score="logloss")
  D_int = udep_int$distance |> as.matrix()
  t2=proc.time()
  t_x_N2[ite,1]=t2[3]-t1[3]
  t1=proc.time()
  udep_no_int = mdist(x=dfc,preset="custom", distance_cont = "manhattan",
                      distance_cat = "tot_var_dist",
                      commensurable = TRUE,
                      scaling_cont = 'std',interaction = FALSE)
  D_no_int = udep_no_int$distance |> as.matrix()
  t2=proc.time()
  t_x_N2[ite,2]=t2[3]-t1[3]
  t1=proc.time()
  # D_gudmm = distance_by_method(dfc,method="gudmm",interaction=FALSE) |> as.matrix()
  naive = mdist(x=dfc,preset="euclidean_onehot")
  D_naive = naive$distance |> as.matrix()
  t2=proc.time()
  t_x_N2[ite,3]=t2[3]-t1[3]
  t1=proc.time()
  gow = mdist(x=dfc,preset="gower")
  D_gow = gow$distance |> as.matrix()
  t2=proc.time()
  t_x_N2[ite,4]=t2[3]-t1[3]

  t1=proc.time()
  mod_gow = mdist(x=dfc,preset="mod_gower")
  D_mod_gow = mod_gow$distance |> as.matrix()
  t2=proc.time()
  t_x_N2[ite,5]=t2[3]-t1[3]




  cl_int = spectral_from_dist(D_int, k=myk,affinity_method = "gaussian")
  cl_no_int = spectral_from_dist(D_no_int, k=myk,affinity_method = "gaussian")
  # cl_gudmm = spectral_from_dist(D_gudmm, k=myk,affinity_method = "gaussian")
  cl_naive = spectral_from_dist(D_naive, k=myk,affinity_method = "gaussian")
  cl_gow = spectral_from_dist(D_gow, k=myk,affinity_method = "gaussian")
  cl_mod_gow = spectral_from_dist(D_mod_gow, k=myk,affinity_method = "gaussian")

  #res <- clustMD(X = dfc, G = myk, CnsIndx = 15, OrdIndx = 23, model = "EVI")

  #ite=1
  ARIM_x_N2[ite,1]=adjustedRandIndex(cl_int, truth)
  ARIM_x_N2[ite,2]=adjustedRandIndex(cl_no_int, truth)
  #adjustedRandIndex(cl_gudmm, truth)
  ARIM_x_N2[ite,3]=adjustedRandIndex(cl_gow, truth)
  ARIM_x_N2[ite,4]=adjustedRandIndex(cl_mod_gow, truth)
  ARIM_x_N2[ite,5]=adjustedRandIndex(cl_naive, truth)
  print(ARIM_x_N2[ite,])
}
ARIM_x_N2

names=c('Udep Int.','Udep','Gower','mod G.','naive')
colnames(ARIM_x_N2)=names


save(t_x_N2,file='t_x_N2.rdata')
save(ARIM_x_N2,file='ARIM_x_N2.rdata')



boxplot(ARIM_x_N2, main = "", ylab='ARI',
        las = 2,                # Rotate axis labels for readability
        col = "lightblue",      # Fill color
        border = "darkblue",    # Border color
        outline = TRUE)         # Show outliers
