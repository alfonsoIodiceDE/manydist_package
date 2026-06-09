library(archdata)
data(DartPoints)
library(manydist)
#library(cluster)
#library(FactoMineR)
#library(aricode)
#library(tidymodels)
#library(manydist)
#library(smacof)

df<-DartPoints

# The B.Width (basal width) variable has 22 missing values
# out of 91. The remaining six categorical shape variables have at most two missing values each. Carlson (2017, p. 148) notes that excluding B.Width and the three cases with other missing values yields 88 complete observations.
df<-df[,-c(2,3,4)]
df<-na.omit(df) # If all variables used, this results in 68 observations
# Alternatively: remove B.width, col 8:
df<-DartPoints
df<-df[,-c(2,3,4,8)]
df<-na.omit(df)  # 88 obs

### LOVOs for supervised methods:
distcustoms<-c("tvd","gifi_chi2","le_and_ho","kumar-johnson")
resp<-TRUE
comm<-TRUE

meths<- list(tvd_sup = lovo_method_spec(
  response = Name,
  distance_cat = "tvd",
  response_used = resp,
  commensurable = comm
),
gifi_sup = lovo_method_spec(
  response = Name,
  distance_cat = "gifi_chi2",
  response_used = resp,
  commensurable = comm
),
leho_sup = lovo_method_spec(
  response = Name,
  distance_cat = "le_and_ho",
  response_used = resp,
  commensurable = comm
),
kj_sup = lovo_method_spec(
  response = Name,
  distance_cat = "kumar-johnson",
  response_used = resp,
  commensurable = comm
))

lovo_cmp <- compare_lovo_mdist(df,methods = meths)#,cluster_k=5)

autoplot(lovo_cmp)
autoplot(lovo_cmp,"mad_importance")
