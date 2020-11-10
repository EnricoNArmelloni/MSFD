library(tidyverse)
library(data.table)
library(writexl)
library(segmented)
library(MASS)
library(caret)
rm(list=ls()) 
set.seed(123)

###--- Set parameters
mydir_input="~/CNR/MSFD/github/MSFD/input" ## This folder must contain the "settings" file and subfolder paths as "/GSA/Species"

mydir_out_base="~/CNR/MSFD/github/MSFD/output/models_trial" 

area="GSA17_18"

sel_ind<-"L95"

species_list=c("MUT", "HKE")

filename="_L95_plot.csv"

settings <- read_csv(paste0(mydir_input, "/settings.csv"))%>%
  gather(key="var_id", value="variable", -species_id)%>%
  dplyr::filter(!is.na(variable))%>%
  dplyr::select(-var_id)

for(i in 1:length(species_list)){

  species=species_list[i]
  
  ###  Run ####
  if(dir.exists(paste0(mydir_out_base,"/",area, "_", species))==TRUE){
    NA
  }else{
    dir.create(paste0(mydir_out_base,"/",area, "_",species))
  }
  mydir_out=paste0(mydir_out_base,"/",area, "_",species)
  
  sets=settings%>%dplyr::filter(species_id == species)
  
  variables=unique(sets$variable)
  
  file_import=paste0(species,"_", area, filename)
  
    data= read_csv(paste0( mydir_input, "/",area, "/", species,"/", file_import))%>%
      dplyr::select(Indicator, variables)
 
    
  
  ## Piecewise regression ####
  # Model selection
    if(length(variables)>1){
      train.control <- trainControl(method = "cv", number = 3) # Set up repeated k-fold cross-validation
      
      step.model <- train(Indicator ~., data = data,
                          method = "leapBackward", 
                          tuneGrid = data.frame(nvmax = 1:3),
                          trControl = train.control)
      
      best_model=step.model[["finalModel"]][["xnames"]][-1]
      
      form=as.formula(paste("Indicator ~ ",paste(best_model,collapse="+"),sep = ""))
      
      
      lm0=glm(form, data=data)
      
    }else{
      
      lm0=glm(Indicator ~ year, data=data)
      
    }
      
      summary(lm0)
      png(filename=paste0(mydir_out, "/",  sel_ind, species, "_summary_model.png"))
      par(mfrow=c(2,2))
      plot(lm0)
      dev.off()
   
   

  
  # regression
  my.seg <- segmented(lm0, seg.Z = ~ year)
  
  png(filename=paste0(mydir_out, "/",  sel_ind, species,  "_segmented.png"))
  par(mfrow=c(1,1))
  plot(my.seg)
  dev.off()
  
  plot(my.seg$y)
  
  
  ###### Test
  bp=my.seg$psi
  
  AP=data%>%
    dplyr::filter(year > round(as.numeric(bp[2])))
  BRP=data%>%
    dplyr::filter(year <= round(as.numeric(bp[2])))
  
  variance_sign<-var.test(BRP$Indicator, AP$Indicator)[["p.value"]]
  
  if(variance_sign >= 0.05){
    ttest<-t.test(BRP$Indicator , AP$Indicator,  var.equal = TRUE)
  }else{
    ttest<-t.test(BRP$Indicator , AP$Indicator,  var.equal = FALSE)
  }
  meanBRP=mean(BRP$Indicator)
  meanAP=mean(AP$Indicator)
  regr.p=ttest$p.value  
  AP$Val=AP$Indicator
  BRP$Val=BRP$Indicator
  
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
  
  res=tibble(species=paste(species, area, sep="_"),formula= paste(best_model,collapse="+"), bp.year=bp[2], slope.1= slope(my.seg)$year[1,1], slope.2=slope(my.seg)$year[2,1], bp.fitted=dev_expl,mean.BRP= meanBRP, meanAP= meanAP, bp.test= regr.p , spearman.p= sprmn.p, rho.spearman=rho)
  
  ### Final plot ####
  newdat=bind_cols(data, fit=my.seg$fitted.values)
  
  dataplot=bind_rows(newdat%>%dplyr::group_by(year)%>%dplyr::summarize(L95=mean(Indicator))%>%dplyr::mutate(source="observed"),
                     newdat%>%dplyr::group_by(year)%>%dplyr::summarize(L95=mean(fit))%>%dplyr::mutate(source="fitted"))
  
  
  if(nrow(dataplot%>%distinct(year)) > 20 ){
    spaz=5
  }else{
    spaz=2
  }
  
  baseplot=function(dat, clr, lbl){
    
    ggplot(data=dat)+
      geom_line(aes(x=year, y=L95, color=source))+
      geom_vline(aes(xintercept=bp[2]), linetype="dashed")+
      geom_rect(aes(xmin = min(AP$year), xmax = max(AP$year), ymin =mean(AP$Val)-sd(AP$Val), ymax = mean(AP$Val)+sd(AP$Val)),alpha = 0.005, fill = clr)+
      geom_segment(x=min(AP$year),xend=max(AP$year),y=mean(AP$Val), yend=mean(AP$Val), size=1, linetype = "dashed", color=clr)+
      geom_rect(aes(xmin = min(BRP$year), xmax = max(BRP$year), ymin =mean(BRP$Val)-sd(BRP$Val), ymax = mean(BRP$Val)+sd(BRP$Val)),alpha = 0.005, fill = "blue")+
      geom_segment(x=min(BRP$year),xend=max(BRP$year),y=mean(BRP$Val), yend=mean(BRP$Val), size=1, linetype = "dashed", color="blue")+
      theme_bw()+scale_x_continuous(breaks=seq(min(data$year),max(data$year),spaz))+
      theme_classic()+
      annotate("text", x=max(dat$year-spaz), y=max(dat$L95)-((max(dat$L95)-min(dat$L95))*0.02), label= lbl)
  }
  
  if(regr.p >= 0.05){
    
    baseplot(dataplot, "blue", "BPA result: AP = RP")
    
  }else if(regr.p < 0.05){
    
    if(meanAP > meanBRP){
      
      baseplot(dataplot, "green", "BPA result: AP > RP")
      
    }else{
      
      baseplot(dataplot, "red", "BPA result: AP < RP")
      
    }
    
  }
  
  ### Save ####
  
  writexl::write_xlsx(res,paste0(mydir_out, "/",  sel_ind, species,  "_results.xlsx") )
  saveRDS(list(lm0, my.seg), file=paste0(mydir_out, "/",  sel_ind, species,  "_models.RData"))
  ggsave(file=paste0(mydir_out, "/",  sel_ind, species, "_fitted_vs_observed.png"), width=12)
  
}
