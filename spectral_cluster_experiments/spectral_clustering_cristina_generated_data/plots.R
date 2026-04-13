library(ggplot2)
library(reshape2)
names=c('Udep Int.','Udep','Gower','mod G.','naive', 'dkss' )
# Convert to data frame with labels
df1 <- as.data.frame(ARIM_N1E)
df2 <- as.data.frame(ARIM_N2E)

df1$Group <- "500"
df2$Group <- "1000"

df1$Obs <- 1:nrow(df1)
df2$Obs <- 1:nrow(df2)

# Combine
df <- rbind(df1, df2)

# Melt to long format
df_long <- melt(df, id.vars = c("Group", "Obs"),
                variable.name = "Method", value.name = "ARI")

# Plot
ggplot(df_long, aes(x = interaction(Group, Method), y = ARI, fill = Group)) +
  geom_boxplot() +
  scale_x_discrete(labels = rep(names, each = 2)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Method", y = "ARI", fill = 'n',
       title = "15 continuous and 15 categorical variables") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(hjust = 0.5)
  )
####### plot 2
# Convert to data frame with labels
df1 <- as.data.frame(ARIM_N1_1020)
df2 <- as.data.frame(ARIM_N2_1020)

df1$Group <- "500"
df2$Group <- "1000"

df1$Obs <- 1:nrow(df1)
df2$Obs <- 1:nrow(df2)

# Combine
df <- rbind(df1, df2)

# Melt to long format
df_long <- melt(df, id.vars = c("Group", "Obs"),
                variable.name = "Method", value.name = "ARI")

# Plot
ggplot(df_long, aes(x = interaction(Group, Method), y = ARI, fill = Group)) +
  geom_boxplot() +
  scale_x_discrete(labels = rep(names, each = 2)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Method", y = "ARI", fill = 'n',
       title = "10 continuous and 20 categorical variables") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

####### plot 3
# Convert to data frame with labels
df1 <- as.data.frame(ARIM_N1_2010)
df2 <- as.data.frame(ARIM_N2_2010)

df1$Group <- "500"
df2$Group <- "1000"

df1$Obs <- 1:nrow(df1)
df2$Obs <- 1:nrow(df2)

# Combine
df <- rbind(df1, df2)

# Melt to long format
df_long <- melt(df, id.vars = c("Group", "Obs"),
                variable.name = "Method", value.name = "ARI")

# Plot
ggplot(df_long, aes(x = interaction(Group, Method), y = ARI, fill = Group)) +
  geom_boxplot() +
  scale_x_discrete(labels = rep(names, each = 2)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Method", y = "ARI", fill = 'n',
       title = "20 continuous and 10 categorical variables") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(hjust = 0.5)
  )


#######Plot 4
df1 <- as.data.frame(ARIM_noise)
df2 <- as.data.frame(ARIM_noise_N2)

df1$Group <- "500"
df2$Group <- "1000"

df1$Obs <- 1:nrow(df1)
df2$Obs <- 1:nrow(df2)

# Combine
df <- rbind(df1, df2)

# Melt to long format
df_long <- melt(df, id.vars = c("Group", "Obs"),
                variable.name = "Method", value.name = "ARI")

# Plot
ggplot(df_long, aes(x = interaction(Group, Method), y = ARI, fill = Group)) +
  geom_boxplot() +
  scale_x_discrete(labels = rep(names, each = 2)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Method", y = "ARI", fill = 'n',
    title = "10 continuous, 10 categorical variables, and 10 noise"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(hjust = 0.5)
  )
###### tables

Aris=cbind(apply(ARIM_N1E,2,mean),apply(ARIM_N2E,2,mean),apply(ARIM_N1_1020,2,mean),apply(ARIM_N2_1020,2,mean),apply(ARIM_N1_2010,2,mean),apply(ARIM_N2_2010,2,mean),apply(ARIM_noise,2,mean),apply(ARIM_noise_N2,2,mean))
timess=cbind(apply(t_N1E,2,mean),apply(t_N2E,2,mean),apply(t_N1_1020,2,mean),apply(t_N2_1020,2,mean),apply(t_N1_2010,2,mean),apply(t_N2_2010,2,mean),apply(t_noise,2,mean),apply(t_noise_N2,2,mean))


#######Plot 6
df1 <- as.data.frame(ARIM_N1_NC_48)
df2 <- as.data.frame(ARIM_N2_NC_48)

df1$Group <- "500"
df2$Group <- "1000"

df1$Obs <- 1:nrow(df1)
df2$Obs <- 1:nrow(df2)

# Combine
df <- rbind(df1, df2)

# Melt to long format
df_long <- melt(df, id.vars = c("Group", "Obs"),
                variable.name = "Method", value.name = "ARI")

# Plot
ggplot(df_long, aes(x = interaction(Group, Method), y = ARI, fill = Group)) +
  geom_boxplot() +
  scale_x_discrete(labels = rep(names, each = 2)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Method", y = "ARI", fill = 'n',
    title = "4 continuous and 8 categorical variables"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

df1 <- as.data.frame(ARIM_N1_NC)
df2 <- as.data.frame(ARIM_N2_NC)

df1$Group <- "500"
df2$Group <- "1000"

df1$Obs <- 1:nrow(df1)
df2$Obs <- 1:nrow(df2)

# Combine
df <- rbind(df1, df2)

# Melt to long format
df_long <- melt(df, id.vars = c("Group", "Obs"),
                variable.name = "Method", value.name = "ARI")

# Plot
ggplot(df_long, aes(x = interaction(Group, Method), y = ARI, fill = Group)) +
  geom_boxplot() +
  scale_x_discrete(labels = rep(names, each = 2)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Method", y = "ARI", fill = 'n',
    title = "4 continuous and 4 categorical variables"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(hjust = 0.5)
  )


#######Plot 7
df1 <- as.data.frame(ARIM_N1_NC_42)
df2 <- as.data.frame(ARIM_N2_NC_42)

df1$Group <- "500"
df2$Group <- "1000"

df1$Obs <- 1:nrow(df1)
df2$Obs <- 1:nrow(df2)

# Combine
df <- rbind(df1, df2)

# Melt to long format
df_long <- melt(df, id.vars = c("Group", "Obs"),
                variable.name = "Method", value.name = "ARI")

# Plot
ggplot(df_long, aes(x = interaction(Group, Method), y = ARI, fill = Group)) +
  geom_boxplot() +
  scale_x_discrete(labels = rep(names, each = 2)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Method", y = "ARI", fill = 'n',
    title = "4 continuous and 8 categorical variables"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

df1 <- as.data.frame(ARIM_N1_NC)
df2 <- as.data.frame(ARIM_N2_NC)

df1$Group <- "500"
df2$Group <- "1000"

df1$Obs <- 1:nrow(df1)
df2$Obs <- 1:nrow(df2)

# Combine
df <- rbind(df1, df2)

# Melt to long format
df_long <- melt(df, id.vars = c("Group", "Obs"),
                variable.name = "Method", value.name = "ARI")

# Plot
ggplot(df_long, aes(x = interaction(Group, Method), y = ARI, fill = Group)) +
  geom_boxplot() +
  scale_x_discrete(labels = rep(names, each = 2)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Method", y = "ARI", fill = 'n',
    title = "4 continuous and 2 categorical variables"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

