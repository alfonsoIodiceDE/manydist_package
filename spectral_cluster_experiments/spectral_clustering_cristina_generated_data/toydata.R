library("mvtnorm")
library("reshape2")
library("tidyverse")

pnum=4
pcat=2
p=pnum+pcat
n=500



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
  truth =c(rep(1,n1),rep(2,n2),rep(3,n3))

  colnames(Xkp)=paste0("V", 1:6)
  pairs(Xkp[,1:4],col=truth,pch=truth)
  library(xtable)
  xtable(cbind(table(Xkp[1:n1,5],Xkp[1:n1,6]),table(Xkp[(n1+1):n2,5],Xkp[(n1+1):n2,6]),  table(Xkp[(n2+1):n3,5],Xkp[(n2+1):n3,6])))

  library(tidyverse)

  # add cluster label
  Xkp$cluster <- factor(truth, labels = c("C1","C2","C3"))

  # function to compute conditional means
  compute_profiles <- function(data, cat_var) {
    data %>%
      group_by(cluster, !!sym(cat_var)) %>%
      summarise(across(V1:V4, mean), .groups = "drop") %>%
      pivot_longer(V1:V4, names_to = "variable", values_to = "value") %>%
      mutate(category = as.factor(.data[[cat_var]]),
             cat_var = cat_var)
  }

  # profiles for V5 and V6
  prof_V5 <- compute_profiles(Xkp, "V5")
  prof_V6 <- compute_profiles(Xkp, "V6")

  profiles <- bind_rows(prof_V5, prof_V6)

  diff_df <- profiles %>%
    mutate(category = as.character(category)) %>%
    select(cluster, cat_var, variable, category, value) %>%
    pivot_wider(
      names_from = category,
      values_from = value
    ) %>%
    mutate(diff = `1` - `0`) %>%
    filter(!is.na(diff))

  # plot
  fig2plot = ggplot(diff_df, aes(x = variable, y = diff, fill = diff > 0)) +
    geom_col(width = 0.5, show.legend = FALSE, alpha = .75) +
    scale_fill_manual(values = c("indianred", "dodgerblue")) +
    facet_grid(cluster ~ cat_var) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
      x = "Continuous variables",
      y = "Conditional mean difference (1 - 0)"
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "grey95"),
      strip.text = element_text(face = "bold")
    )


  ggsave(
    filename = "figure2.png",
    plot = fig2plot,
    width = 7,
    height = 5,
    dpi = 300,
    device = ragg::agg_png
  )

  ggsave(
    "figure2.pdf",
    plot = fig2plot,
    width = 7,
    height = 5
  )

  #  ##### C 1
  # data=Xkp[1:n1,1:4]
  # group =Xkp[1:n1,5]
  #
  # df1 <- as.data.frame(data[which(group==0),])
  # df2 <- as.data.frame(data[which( group==1),])
  #
  # df1$Group <- "0"
  # df2$Group <- "1"
  #
  # df1$Obs <- 1:nrow(df1)
  # df2$Obs <- 1:nrow(df2)
  #
  # # Combine
  # df <- rbind(df1, df2)
  #
  # # Melt to long format
  # df_long <- reshape2::melt(df, id.vars = c("Group", "Obs"),
  #                 variable.name = "Method", value.name = "Value")
  #
  # # Plot
  # ggplot(df_long, aes(x = interaction(Group, Method), y = Value, fill = Group)) +
  #   geom_boxplot() +
  #   scale_x_discrete(labels = rep(paste0("V", 1:4), each = 2)) +
  #   labs(x = "Variables", y = "", fill = 'V5',
  #        title = "Cluster 1") +
  #   theme_bw() +
  #   theme(
  #     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
  #     plot.title = element_text(hjust = 0.5)
  #   )
  #
  #
  # #C1_2
  # ##### C 1
  # data=Xkp[1:n1,1:4]
  # group =Xkp[1:n1,6]
  #
  # df1 <- as.data.frame(data[which(group==0),])
  # df2 <- as.data.frame(data[which( group==1),])
  #
  # df1$Group <- "0"
  # df2$Group <- "1"
  #
  # df1$Obs <- 1:nrow(df1)
  # df2$Obs <- 1:nrow(df2)
  #
  # # Combine
  # df <- rbind(df1, df2)
  #
  # # Melt to long format
  # df_long <- melt(df, id.vars = c("Group", "Obs"),
  #                 variable.name = "Method", value.name = "Value")
  #
  # # Plot
  # ggplot(df_long, aes(x = interaction(Group, Method), y = Value, fill = Group)) +
  #   geom_boxplot() +
  #   scale_x_discrete(labels = rep(paste0("V", 1:4), each = 2)) +
  #   labs(x = "Variables", y = "", fill = 'V6',
  #        title = "Cluster 1 ") +
  #   theme_bw() +
  #   theme(
  #     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
  #     plot.title = element_text(hjust = 0.5)
  #   )
  #
  #
  #
  #
  #
  #
  # ##############C2
  #
  #
  #
  # ##### C 1
  # data=Xkp[(n1+1):n2,1:4]
  # group =Xkp[(n1+1):n2,5]
  #
  # df1 <- as.data.frame(data[which(group==0),])
  # df2 <- as.data.frame(data[which( group==1),])
  #
  # df1$Group <- "0"
  # df2$Group <- "1"
  #
  # df1$Obs <- 1:nrow(df1)
  # df2$Obs <- 1:nrow(df2)
  #
  # # Combine
  # df <- rbind(df1, df2)
  #
  # # Melt to long format
  # df_long <- melt(df, id.vars = c("Group", "Obs"),
  #                 variable.name = "Method", value.name = "Value")
  #
  # # Plot
  # ggplot(df_long, aes(x = interaction(Group, Method), y = Value, fill = Group)) +
  #   geom_boxplot() +
  #   scale_x_discrete(labels = rep(paste0("V", 1:4), each = 2)) +
  #   labs(x = "Variables", y = "", fill = 'V5',
  #        title = "Cluster 2") +
  #   theme_bw() +
  #   theme(
  #     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
  #     plot.title = element_text(hjust = 0.5)
  #   )
  #
  #
  # #C1_2
  # ##### C 1
  # data=Xkp[(n1+1):n2,1:4]
  # group =Xkp[(n1+1):n2,6]
  #
  # df1 <- as.data.frame(data[which(group==0),])
  # df2 <- as.data.frame(data[which( group==1),])
  #
  # df1$Group <- "0"
  # df2$Group <- "1"
  #
  # df1$Obs <- 1:nrow(df1)
  # df2$Obs <- 1:nrow(df2)
  #
  # # Combine
  # df <- rbind(df1, df2)
  #
  # # Melt to long format
  # df_long <- melt(df, id.vars = c("Group", "Obs"),
  #                 variable.name = "Method", value.name = "Value")
  #
  # # Plot
  # ggplot(df_long, aes(x = interaction(Group, Method), y = Value, fill = Group)) +
  #   geom_boxplot() +
  #   scale_x_discrete(labels = rep(paste0("V", 1:4), each = 2)) +
  #   labs(x = "Variable", y = "", fill = 'V6',
  #        title = "Cluster 2") +
  #   theme_bw() +
  #   theme(
  #     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
  #     plot.title = element_text(hjust = 0.5)
  #   )
  #
  #
  #
  #
  #
  #
  # ########C3
  #
  #
  #
  # ##############C2
  #
  #
  #
  # ##### C 1
  # data=Xkp[(n2+1):n3,1:4]
  # group =Xkp[(n2+1):n3,5]
  #
  # df1 <- as.data.frame(data[which(group==0),])
  # df2 <- as.data.frame(data[which( group==1),])
  #
  # df1$Group <- "0"
  # df2$Group <- "1"
  #
  # df1$Obs <- 1:nrow(df1)
  # df2$Obs <- 1:nrow(df2)
  #
  # # Combine
  # df <- rbind(df1, df2)
  #
  # # Melt to long format
  # df_long <- melt(df, id.vars = c("Group", "Obs"),
  #                 variable.name = "Method", value.name = "Value")
  #
  # # Plot
  # ggplot(df_long, aes(x = interaction(Group, Method), y = Value, fill = Group)) +
  #   geom_boxplot() +
  #   scale_x_discrete(labels = rep(paste0("V", 1:4), each = 2)) +
  #   labs(x = "Variables", y = "", fill = 'V5',
  #        title = "Cluster 3") +
  #   theme_bw() +
  #   theme(
  #     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
  #     plot.title = element_text(hjust = 0.5)
  #   )
  #
  #
  # #C1_2
  # ##### C 1
  # data=Xkp[(n2+1):n3,1:4]
  # group =Xkp[(n2+1):n3,6]
  #
  # df1 <- as.data.frame(data[which(group==0),])
  # df2 <- as.data.frame(data[which( group==1),])
  #
  # df1$Group <- "0"
  # df2$Group <- "1"
  #
  # df1$Obs <- 1:nrow(df1)
  # df2$Obs <- 1:nrow(df2)
  #
  # # Combine
  # df <- rbind(df1, df2)
  #
  # # Melt to long format
  # df_long <- melt(df, id.vars = c("Group", "Obs"),
  #                 variable.name = "Method", value.name = "Value")
  #
  # # Plot
  # ggplot(df_long, aes(x = interaction(Group, Method), y = Value, fill = Group)) +
  #   geom_boxplot() +
  #   scale_x_discrete(labels = rep(paste0("V", 1:4), each = 2)) +
  #   labs(x = "Variable", y = " ", fill = 'V6',
  #        title = "Cluster 3") +
  #   theme_bw() +
  #   theme(
  #     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
  #     plot.title = element_text(hjust = 0.5)
  #   )
  #
  #
  #
