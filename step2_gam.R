library(mgcv)
library(fmsb)
library(Hmisc)
library(car)
library(tidyverse)
require(magrittr)

rm(list = ls())

# setting the working directory and select the datasets
#### Settings
input_dir="~/CNR/MSFD/github/release"

sspp= "RJC" ## set species 3 alpha code

gsa="17" # Assign gsa code to output file name (i.e. 9_11_ for two GSAs) 

sel="sizethreshold" ## set kind of input file. Alternatives: "selection" ; "total; sizethreshold"

indicator="L95" # alternatives= L95 ; Pmega

st_month=6 ## Select month for standardization

### set WD #####
dir_t=paste0(input_dir,"/","GSA",gsa,"/", sspp,"/", sep="")
dir.create(paste0(dir_t, "/gam"))

## Import data
pop=ifelse(sel=="total", "_whole_", "_mat_")

data= read_csv(paste0(dir_t, sspp,"_GSA",gsa,"_Indicators_", sel,".csv"))%>%
  dplyr::mutate(month=lubridate::month(as.Date(mean_doy, origin=paste0(as.numeric(year), "-01-01"))))%>%
  dplyr::select(year, nhauls, month, Indicator=indicator)


#### PAIR PLOT FOR COLLINEARITY CHECK
panel.smooth2<-function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
                         cex = 1, col.smooth = "red", span = 2/3, iter = 3, ...) 
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
          col = 1, ...)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor ,...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  #r <- abs(cor(x, y))
  r <- (cor(x, y))
  #txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- format(c(r, 0.123456789), digits = 1)[1]
  txt <- paste(prefix, txt, sep = "")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  #text(0.5, 0.5, txt, cex = cex.cor * r)
  text(0.5, 0.5, txt, cex = 3*abs(r))
}

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "white", ...)
}

png(paste0(dir_t, paste0("/gam/",indicator,"pairplot.png")), height = 20, width = 20, units = "cm", res = 600) 
#op<-par(mfrow=c(1,1), mar=c(4,4,2,2))
Z<-data%>%dplyr::select(nhauls, month,year)
pairs(Z,las=1 ,lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist, cex=0.9, cex.labels=1.0)
# box()
dev.off()


# VIF analysis
mod_lm <- lm(Indicator ~ nhauls+ month+year, data=data)
vif<-vif(lm(mod_lm, data = data)); vif

# correlation matrix
rcorr(as.matrix(data)) 

# Stepwise forward inclusion for GAM modeling
mod1 <- gam(Indicator ~ s(year)+0, family= gaussian (link="identity"),data=data, select=T) 
mod2 <- gam(Indicator ~ factor( month)+0, family= gaussian (link="identity"),data=data, select=T)
mod3 <- gam(Indicator ~ s(nhauls)+0, family= gaussian (link="identity"),data=data, select=T)


summary(mod1) 
summary(mod2) 
summary(mod3) 
AIC(mod1,mod2,mod3)

mod11 <- gam(Indicator ~ s(year)+ factor(month), family= gaussian (link="identity"),data=data, select=T)
mod12 <- gam(Indicator ~ s(year)+ s(nhauls,k=3), family= gaussian (link="identity"),data=data, select=T)

summary(mod11) 
summary(mod12) 
AIC(mod11,mod12)

mod21 <- gam(Indicator ~ s(year)+ s(nhauls) + s(month,k=3), family= gaussian (link="identity"),data=data, select=T)
summary(mod21)
AIC(mod21)

## GAM.checks
par(mfrow=c(2,2))
gam.check(mod11)
gam.check(mod12)
gam.check(mod21)

## factors effect
par(mfrow=c(1,2))
plot(mod11)
termplot(mod11)
plot(mod12)
par(mfrow=c(1,3))
plot(mod21)


##### Model selection
mod=mod11

#residui pearson
png(paste0(dir_t, "/gam/", indicator,pop,"pearson_res.png"), height = 20, width = 20, units = "cm",  res = 600) 
par(mfrow=c(1,1))
E.m4<- resid(mod, type = "pearson")
Fit.m4 <- fitted(mod)
#e4<-resid(M3)
plot(x = Fit.m4, y = E.m4, xlab = "Fitted values",
     ylab = "Residuals")
#plot(x = Fit.m4, y = E.m4, xlab = "Fitted values",
#    ylab = "Residuals", main = "logKG.KM2")
abline(0,0)
tmp <- loess(E.m4 ~Fit.m4 ,span=0.75)
tmp2 <- predict(tmp,se=T)
I1 <- order(Fit.m4)
lines(Fit.m4[I1], tmp2$fit[I1], lty=1)
lines(Fit.m4[I1], tmp2$fit[I1] + 2*tmp2$se.fit[I1], lty = 2)
lines(Fit.m4[I1], tmp2$fit[I1] - 2*tmp2$se.fit[I1], lty = 2)
dev.off()


# ASSUNZIONE SU INDIPENDENZA, SI PLOTTANO I RESIDUI CONTRO LA VARIABILE INDIPENDENTE:

# ASSUNZIONE SU INDIPENDENZA, SI PLOTTANO I RESIDUI CONTRO LA VARIABILE INDIPENDENTE:
png(paste0(dir_t, "/gam/", indicator, pop, "autocorrelation.png"), height = 20, width = 30, units = "cm",  res = 600) 
par(mfrow = c(1,1))
acf(residuals(mod))
dev.off()


plot_variabili=function(beta, teta){
  var1 <- beta
  plot(x = var1, y = E.m4, xlab = teta,
       ylab = "Residuals")
  abline(0,0)
  tmp <- loess(E.m4 ~var1 ,span=0.8)
  tmp2 <- predict(tmp,se=T)
  I1 <- order(var1)
  lines(var1[I1], tmp2$fit[I1], lty=1)
  lines(var1[I1], tmp2$fit[I1] + 2*tmp2$se.fit[I1], lty = 2)
  lines(var1[I1], tmp2$fit[I1] - 2*tmp2$se.fit[I1], lty = 2)
  
}

png(paste0(dir_t, "/gam/", indicator, pop, "res_by_var.png"), height = 20, width = 30, units = "cm",  res = 600) 
par(mfrow=c(1,3))
plot_variabili(data$year, "Year")
plot_variabili(data$month, "Month")
plot_variabili(data$nhauls, "nhauls")
dev.off()

png(paste0(dir_t, "/gam/", indicator,pop,  "res_check.png"), height = 20, width = 30, units = "cm", res = 600) 
par(mfrow=c(1,3))
#Normality
par(mar=c(5,5,2,2))
hist(E.m4, ylab = "Frequency", xlab = "Residuals", las=1,breaks=16, cex.lab=1.1, cex.axis=1.1, main=NULL)

#Or qq-plot
par(mar=c(5,5,2,2))
qqnorm(E.m4, main=NULL, lwd=1.5,cex.lab=1.1, las=1,cex.axis=1, bty="l", xlab="Theoretical quantiles", ylab="Std. Dev. quantiles")
abline(0.0,1.3, lty=1, lwd=1.5)


res <- residuals(mod)
plot(res, ylim=c(-5.5,5.5))
abline(0,0, col="red")
dev.off()

png(paste0(dir_t, "/gam/", indicator, "cook.png"), height = 20, width = 30, units = "cm",  res = 600) 
# cOOK'S DISTANCE
par(mfrow=c(1,1))
plot(cooks.distance(mod), ylab="Cook's distance", type = "h", ylim=c(0,1))
abline(h=1, col=1,lwd=2)
dev.off()

# inspection of residuals and diagnostic plots



png(paste0(dir_t, "/gam/", indicator,pop,  "gam_check.png"), height = 20, width = 30, units = "cm",  res = 600) 
par(mfrow=c(2,2))
gam.check(mod)
dev.off()

png(paste0(dir_t, "/gam/", indicator,pop,  "factor_effects.png"), height = 20, width = 30, units = "cm",  res = 600) 
par(mfrow=c(1,2))
plot(mod,select =1)
plot(mod, select =2)
plot(mod, select =3)
plot(mod, select =4)
termplot(mod, se=T)
dev.off()


################################# Prediction anbd standardization
xx <- predict(mod, type = "response", se.fit = TRUE)
data$pred <- xx$fit
data$pred.se <- xx$se.fit

stand <- data %>% dplyr::mutate(month = st_month, nhauls = 120)
yy <- predict(mod, newdata=stand, type = "response", se.fit = TRUE)
data$stand <- yy$fit
data$stand.se <- yy$se.fit

data %<>% mutate(pred.lwr = pred - 1.96*pred.se,
                 pred.upr = pred + 1.96*pred.se,
                 stand.lwr = stand - 1.96*stand.se,
                 stand.upr = stand + 1.96*stand.se) 


gam1 <- ggplot()+
  geom_point(aes(x = year, y = Indicator), color = "black", data = data) +
  geom_smooth(data = data, aes(x = year, y = stand, 
                               ymin = stand -1.96 * stand.se,
                               ymax = stand + 1.96 * stand.se, 
                               color = "Standardized on month 6",
                               fill = "Standardized on month 6"), 
              stat = "identity") +
  geom_smooth(data = data, aes(x = year, y = pred, 
                               ymin = pred -1.96 * pred.se,
                               ymax = pred + 1.96 * pred.se, 
                               color = "Predicted",
                               fill = "Predicted"), 
              stat = "identity") +
  
  theme_bw()+
  ylab("L95 (mm)") +
  xlab("Year") +
  scale_x_continuous(breaks = seq(1994,2020,2)) +
  scale_color_discrete(name = "Source", 
                       labels = c("Predicted", paste("Standardized on month",st_month)))+
  scale_fill_discrete(name = "Source", 
                      labels = c("Predicted", paste("Standardized on month",st_month)))+
  ggtitle(paste(sspp,"GSA",gsa, indicator, pop)) + 
  theme(legend.position="bottom")

ggsave(paste0(dir_t, "/gam/", indicator, pop, "indicator_trend.png"), plot = gam1, width = 14, height = 8.5)

write.csv(data, paste0(dir_t, "gam/",indicator,pop,  sspp,"_standardized_", st_month, ".csv"), row.names = FALSE)


