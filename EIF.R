#only run the first line 1 time
Eredivisie28=log(Eredivisie28)


df = (Eredivisie28)
# Choose type of outlier
age=FALSE
marketvalue=TRUE
mvage= FALSE

if (age==TRUE){
  #EIF AGE
  df$Age[sample(nrow(df),1)] <- 3.5
} else if (marketvalue==TRUE){
  #EIF MARKETVALUE
  df$MarketValue[sample(nrow(df),1)] <- 18
} else if (mvage==TRUE){
  #EIF MARKETVALUE & AGE
  df[sample(nrow(df),1),] <- c(3.5, 18)
}


# deMCD, Age on MarketValue
coef_MCD <-lmDetMCD(data.frame(Eredivisie28$Age), data.frame(Eredivisie28$MarketValue))[[1]]
coef_EIF_MCD <- lmDetMCD(data.frame(df$Age), data.frame(df$MarketValue))[[1]]
# Emperical Influence Function
EIF_MCD = nrow(df)*(coef_EIF_MCD - coef_MCD)

# least-trimmed squares, Age on MarketValue
coef_lts <- ltsReg(Eredivisie28$Age, Eredivisie28$MarketValue)$coefficients
coef_EIF_lts <- ltsReg(df$Age, df$MarketValue)$coefficients
# emperical influence function
EIF_lts = nrow(df)*(coef_EIF_lts - coef_lts)

# OLS, Age on MarketValue
coef_ols = lm(Eredivisie28$MarketValue ~ Eredivisie28$Age)$coefficients
coef_EIF_ols = lm(df$MarketValue ~ df$Age)$coefficients
# emperical influence function
EIF_ols = nrow(df)*(coef_EIF_ols - coef_ols)

plots = rbind(coef_lts, coef_ols, coef_MCD[,])
ggplot(data=Eredivisie28, aes(y=MarketValue,x=Age))+geom_point()+
  geom_abline(aes(intercept=plots[1,1], slope=plots[1,2], color='LTS'))+
  geom_abline(aes(intercept=plots[2,1], slope=plots[2,2], color='OLS'))+
  geom_abline(aes(intercept=plots[3,1], slope=plots[3,2], color='detMCD'))+
  theme(legend.title = element_blank())+
  labs(y='log(MarketValue)', x='log(Age)')+
  scale_colour_brewer(palette="Set1")


plotsEIF = rbind(coef_EIF_lts, coef_EIF_ols, coef_EIF_MCD[,])
ggplot(data=df, aes(y=MarketValue,x=Age))+geom_point()+
  geom_abline(aes(intercept=plotsEIF[1,1], slope=plotsEIF[1,2], color='LTS'))+
  geom_abline(aes(intercept=plotsEIF[2,1], slope=plotsEIF[2,2], color='OLS'))+
  geom_abline(aes(intercept=plotsEIF[3,1], slope=plotsEIF[3,2], color='detMCD'))+
  theme(legend.title = element_blank())+
  labs(y='log(MarketValue)', x='log(Age)')+
  scale_colour_brewer(palette="Set1")
  

