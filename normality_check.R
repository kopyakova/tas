library(ggplot2)
library(normtest)

age <- Eredivisie28[,1]
market_value <- Eredivisie28[,2]

#Check normality of the variables
jb.norm.test(age)
jb.norm.test(market_value)

jb.norm.test(log(age))
jb.norm.test(log(market_value))


# Histogram with density plot
ggplot(Eredivisie28, aes(x=Age)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") 

ggplot(Eredivisie28, aes(x=MarketValue)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") 

ggplot(log(Eredivisie28), aes(x=Age)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") 

ggplot(log(Eredivisie28), aes(x=MarketValue)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") 

# Basic scatter plot
gg1 = ggplot(Eredivisie28, aes(x = Age, y=MarketValue)) + geom_point(colour = "dodgerblue4")
gg1
#Scatter plot of log transformed data
log_Eredivisie28 = log(Eredivisie28)
detMCD = covDetMCD(log_Eredivisie28)

log_Eredivisie28 <- cbind(log_Eredivisie28, detMCD$weights)
colnames(log_Eredivisie28) = c(paste0("Log", colnames(Eredivisie28)), "weights")
theme_set(theme_bw())
gg = ggplot(log_Eredivisie28, aes(x = LogAge, y=LogMarketValue)) + geom_point(aes(color= as.factor(weights))) + 
    stat_ellipse(data = log_Eredivisie28[detMCD$weights > 0,]) + stat_ellipse(data = log_Eredivisie28, linetype = 2) +
    xlab("Logarithm of Age") + ylab("Logarithm of Market Value") + scale_color_manual(name="Observations",labels=c("Outlier", "Not outlier"), values = c("dodgerblue4", "firebrick2"))
gg
