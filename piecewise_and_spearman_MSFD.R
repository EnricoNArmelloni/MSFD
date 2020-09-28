library(tidyverse)
library(data.table)
library(writexl)
library(segmented)
library(MASS)
library(caret)
rm(list=ls()) 
set.seed(123)


###--- Set parameters
mydir_input="~/CNR/MSFD/github/MSFD/input"
mydir_out="~/CNR/MSFD/github/MSFD/output/models_trial" 

area="GSA17_18"
sel_ind<-"L95"
timeframe="year" ## alternatives: "year", "haul"

species_list=c("MUT", "HKE")
depth_list=c(200, 600)

variables_haul=c("year", "month", "meandepth", "hour", "Lat", "Lon", "vessel") ### Data_by_haul
variables_year=c("year", "month") ### Data_by_year

filename_haul="_L95_haul.csv"
filename_year="_L95_plot.csv"



for(i in 1:length(species_list)){
  
species=species_list[i]
max_depth=depth_list[i]

###  Run ####
if(dir.exists(paste0(mydir_out,"/",area, "_", species))==TRUE){
  NA
}else{
  dir.create(paste0(mydir_out,"/",area, "_",species))
}
mydir_out=paste0(mydir_out,"/",area, "_",species)


if(timeframe=="year"){
  filename=filename_year
  variables=variables_year
}else{
  filename=filename_haul
  variables=variables_haul
}

L95= read_csv(paste0("~/CNR/MSFD/github/MSFD/input/",area, "/", species,"/",species,"_", area, filename))%>%
  dplyr::select(Indicator, variables)

if(timeframe=="haul"){
  data=L95%>%
    dplyr::filter(meandepth <=max_depth)
}else{
  data=L95
}


## Piecewise regression ####
# Model selection
if(timeframe=="haul"){
  train.control <- trainControl(method = "cv", number = 10) # Set up repeated k-fold cross-validation
  
  step.model <- train(Indicator ~., data = data,
                      method = "leapBackward", 
                      tuneGrid = data.frame(nvmax = 1:5),
                      trControl = train.control)
  
  best_model=step.model[["finalModel"]][["xnames"]][-1]
  
  if(length(best_model[str_detect(best_model, "vessel")==TRUE])>0){
    
    best_model=c(best_model[str_detect(best_model, "vessel")==FALSE], "vessel")
    
  }else{
    NA
  }
  
  form=as.formula(paste("Indicator ~ ",paste(best_model,collapse="+"),sep = ""))
  
}else{
  train.control <- trainControl(method = "cv", number = 3) # Set up repeated k-fold cross-validation
  
  step.model <- train(Indicator ~., data = data,
                      method = "leapBackward", 
                      tuneGrid = data.frame(nvmax = 1:3),
                      trControl = train.control)
  
  best_model=step.model[["finalModel"]][["xnames"]][-1]
  
  form=as.formula(paste("Indicator ~ ",paste(best_model,collapse="+"),sep = ""))
}



lm0=glm(form, data=data)
summary(lm0)
png(filename=paste0(mydir_out, "/",  sel_ind, species,"_", timeframe, "_summary_model.png"))
par(mfrow=c(2,2))
plot(lm0)
dev.off()

# regression
my.seg <- segmented(lm0, seg.Z = ~ year)

png(filename=paste0(mydir_out, "/",  sel_ind, species, "_", timeframe, "_segmented.png"))
par(mfrow=c(1,1))
plot(my.seg)
dev.off()
plot(my.seg$y)
a=my.seg
###### Test
bp=my.seg$psi
AP=data%>%dplyr::filter(year > bp[2])
BRP=data%>%dplyr::filter(year <= bp[2])

variance_sign<-var.test(BRP$Indicator, AP$Indicator)[["p.value"]]

if(variance_sign >= 0.05){
  ttest<-t.test(BRP$Indicator , AP$Indicator,  var.equal = TRUE)
}else{
  ttest<-t.test(BRP$Indicator , AP$Indicator,  var.equal = FALSE)
}
regr.p=ttest$p.value  



## Spearman analysis ####
sprmn=cor.test( ~ Indicator + year, 
          data=data,
          method = "spearman",
          continuity = FALSE,
          conf.level = 0.95)
sprmn.p= sprmn$p.value
rho=     sprmn$estimate[[1]]

#### Result table ####
dev_expl=1-(summary(lm0)[["deviance"]]/summary(lm0)[["null.deviance"]])

res=tibble(species=paste(species, area, sep="_"),formula= paste(best_model,collapse="+"), bp.year=bp[2], slope.1= slope(my.seg)$year[1,1], slope.2=slope(my.seg)$year[2,1], bp.fitted=dev_expl, bp.test= regr.p , spearman.p= sprmn.p, rho.spearman=rho)

### Final plot ####
newdat=bind_cols(data, fit=my.seg$fitted.values)
dataplot=bind_rows(newdat%>%dplyr::group_by(year)%>%dplyr::summarize(L95=mean(Indicator))%>%dplyr::mutate(source="observed"),
          newdat%>%dplyr::group_by(year)%>%dplyr::summarize(L95=mean(fit))%>%dplyr::mutate(source="fitted"))
ggplot(data=dataplot)+geom_line(aes(x=year, y=L95, color=source))+geom_vline(aes(xintercept=bp[2]), linetype="dashed")
fitted(my.seg)
### Save ####

writexl::write_xlsx(res,paste0(mydir_out, "/",  sel_ind, species,"_", timeframe,  "_results.xlsx") )
saveRDS(list(lm0, my.seg), file=paste0(mydir_out, "/",  sel_ind, species, "_", timeframe,  "_models.RData"))
ggsave(file=paste0(mydir_out, "/",  sel_ind, species, "_", timeframe, "_fitted_vs_observed.png"), width=12)
}
