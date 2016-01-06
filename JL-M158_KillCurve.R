df1 <- read.csv("~/R-tables/JL-M158_Luciferase_killcurve.csv", check.names = FALSE)

#############################################################################################################################
## Subsets into two tables - 152 CTRL and 152 CD47 KO
library(reshape2)

df1 <- melt(df1, id=c("Tx","Ratio","Target"))
colnames(df1) <- c("Tx","Ratio","Target", "Rep", "Lum")

#Re-orders table
attach(df1)
df1 <- df1[order(Target, Tx, Ratio, Rep),]
detach(df1)

#############################################################################################################################
## Calculates average of blanks and subtracts
library(plyr)

df1_summary <- ddply(df1, .(Tx, Ratio, Target), summarise, mean=mean(Lum))
df1$Lum <- df1$Lum - df1_summary[df1_summary$Tx %in% "PBS Blank",4]

#############################################################################################################################
##Splits by treatment groups and calculates percent killing
df1_CTRL <- df1[df1$Target %in% "152 CTRL",]
df1_KO <- df1[df1$Target %in% "152 CD47 KO",]

df1_CTRL$Kill <- df1_CTRL$Lum*100 / df1_summary[(df1_summary$Tx == "Tumor Alone") & (df1_summary$Target == "152 CTRL"),4]
df1_KO$Kill <- df1_KO$Lum*100 / df1_summary[(df1_summary$Tx == "Tumor Alone") & (df1_summary$Target == "152 CD47 KO"),4]

#Merges treatment groups together again & removes tumor alone
df2 <- rbind.data.frame(df1_CTRL, df1_KO)
df2 <- df2[df2$Tx != "Tumor Alone",]

#Calculates summary of df2 for plotting
df2_summary <- ddply(df2, .(Tx, Ratio, Target), summarise, mean=mean(Kill), sd=sd(Kill))

#############################################################################################################################
##Calculates t-test p-value for each group

library(plyr)
ttest <- ddply(df2,c("Ratio", "Tx"),
               function(x) {
                 w <- t.test(Kill~Target,data=x)
                 with(w,data.frame(statistic,p.value))
               })

#Reorders ttest table
ttest$Ratio <- ordered(ttest$Ratio, levels = c("1:1", "5:1", "10:1"))
ttest$Tx <- ordered(ttest$Tx, levels = c("Mock", "CpG", "CpG + Etox"))
ttest <- ttest[with(ttest, order(Tx, Ratio)), ]

#Adds significance symbol
ttest$sym <- "ns"
ttest$sym[ttest$p.value <0.05] <- "*"
ttest$sym[ttest$p.value <0.01] <- "**"
ttest$sym[ttest$p.value <0.001] <- "****"

#############################################################################################################################
## Makes faceted line plots of surviving cell percentages (facet by Treatment)
library(ggplot2)

#Sets order of points and faceting
df2_summary$Tx <- factor(df2_summary$Tx, levels = c("Mock", "CpG", "CpG + Etox"))
df2_summary$Target <- factor(df2_summary$Target, levels = c("152 CTRL", "152 CD47 KO"))
df2_summary$Ratio <- factor(df2_summary$Ratio, levels = c("1:1", "5:1", "10:1"))

p <- ggplot(df2_summary, aes(x=factor(Ratio), y=mean, group=Target, colour=Target)) +
  geom_line(size=1) +
  geom_point() +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1, size=1) + 
  scale_colour_manual(
    name  = 'Target Cell',
    values = c("black", "red")
  ) +
  #scale_y_log10(breaks = c(1, 10, 100)) +
  xlab("E:T Ratio") +
  ylab("Surviving cells (%)") +
  facet_grid(. ~Tx) +
  theme(
    panel.border = element_rect(colour = "black", fill=NA),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(size=.01, color="gray"), 
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 18, color = "black", face = "bold"),
    axis.title.y = element_text(size = 18, color = "black", vjust=1.4, face = "bold"),
    strip.text = element_text(size = 18, color = "black", face = "bold"),
    legend.background = element_rect(color = "black"),
    legend.title = element_text(size = 18, color = "black", face = "bold"), 
    legend.text = element_text(size = 14, color = "black"),
    legend.key = element_blank(), 
    legend.position=c(.1,.1)
  )
#############################################################################################################################
## Makes faceted line plots of surviving cell percentages (facet by Cell line)
library(ggplot2)

#Sets order of points and faceting
df2_summary$Tx <- factor(df2_summary$Tx, levels = c("Mock", "CpG", "CpG + Etox"))
df2_summary$Target <- factor(df2_summary$Target, levels = c("152 CTRL", "152 CD47 KO"))
df2_summary$Ratio <- factor(df2_summary$Ratio, levels = c("1:1", "5:1", "10:1"))

q <- ggplot(df2_summary, aes(x=factor(Ratio), y=mean, group=Tx, colour=Tx)) +
  geom_line(size=1) +
  geom_point() +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1, size=1) + 
  scale_colour_manual(
    name  = 'Treatment',
    values = c("black", "red", "blue")
  ) +
  #scale_y_log10(breaks = c(1, 10, 100)) +
  xlab("E:T Ratio") +
  ylab("Surviving cells (%)") +
  facet_grid(. ~Target) +
  theme(
    panel.border = element_rect(colour = "black", fill=NA),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(size=.01, color="gray"), 
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 18, color = "black", face = "bold"),
    axis.title.y = element_text(size = 18, color = "black", vjust=1.4, face = "bold"),
    strip.text = element_text(size = 18, color = "black", face = "bold"),
    legend.background = element_rect(color = "black"),
    legend.title = element_text(size = 18, color = "black", face = "bold"), 
    legend.text = element_text(size = 14, color = "black"),
    legend.key = element_blank(), 
    legend.position=c(.09, 0.15)
  )

#############################################################################################################################
## Plots symbols with graph
#Subsets CTRL data and rearranges for coordinates
df3 <- df2_summary[df2_summary$Target %in% "152 CTRL",]
df3$Ratio <- ordered(df3$Ratio, levels = c("1:1", "5:1", "10:1"))
df3$Tx <- ordered(df3$Tx, levels = c("Mock", "CpG", "CpG + Etox"))
df3 <- df3[with(df3, order(Tx, Ratio)), ]

#Attaches statistic values to coordinate table (df3)
df4 <- cbind.data.frame(df3, ttest[,3:5])
p + geom_text(aes(factor(Ratio), mean + sd + 10, label=sym, group=NULL), size=6, show_guide  = F, data=df4)
