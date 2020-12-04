library(tidyverse)
library(data.table)
library(writexl)
library(segmented)
library(MASS)
library(caret)
rm(list=ls()) 
set.seed(123)

###--- Set parameters
# setting the working directory and select the datasets
#### Settings
sspp= "RJC"

gsa="17" #assign gsa code to input file name (i.e. 9_11_ for two GSAs) 

sel="sizethreshold"## alternatives: "total"; "sizethreshold"

input_dir="~/CNR/MSFD/github/release"

variables="year" #c("year", "month")

indicator="Pmega" # alternatives= L95 ; Pmega

stand="N" # Y if you use data from gam prediction, "N" if you are using data from LFD.

stand_mon=6 # If line 23 is "Y", write here the month on which the data are predicted.
### get WD #####
dir_t=paste0(input_dir, "/","GSA",gsa,"/", sspp,"/", sep="")

dir.create(paste0(dir_t, "/seg_reg"))

#for(i in 1:length(species_list)){

  #species=species_list[i]
  
  ###  Run ####
species=sspp

pop=ifelse(sel=="total", "_whole_", "_mat_")
    
if(stand=="Y"){
  data <- read_csv(paste0(dir_t,"/gam/",indicator,pop, sspp,"_standardized_", stand_mon,".csv"))%>%
    dplyr::select("Indicator"= stand, variables)
}else{
  data <- read_csv(paste0(dir_t,"/",sspp,"_GSA",gsa,"_Indicators_", sel,".csv"))%>%
    dplyr::select("Indicator"= indicator, variables)
}

 
    
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
      
      best_model="year"
      lm0=glm(Indicator ~ year, data=data)
      
    }
      
      summary(lm0)
      png(filename=paste0(dir_t, "/seg_reg", "/", species, "_GSA", gsa,indicator,pop, "_summary_model.png"),height = 20, width = 30, units = "cm", res = 600)
      par(mfrow=c(2,2))
      plot(lm0)
      dev.off()
   
   

  
  # regression
  my.seg <- segmented(lm0, seg.Z = ~ year )
  
  png(filename=paste0(dir_t, "/seg_reg", "/", species, "_GSA", gsa,indicator, pop, "_segmented.png"),height = 20, width = 30, units = "cm", res = 600)
  par(mfrow=c(1,1))
  plot(my.seg)
  dev.off()
  
  
  summary(my.seg)
  
  ###### Test
  bp=my.seg$psi
  
  AP=data%>%
    dplyr::filter(year > round(as.numeric(bp[2])))
  BRP=data%>%
    dplyr::filter(year <= round(as.numeric(bp[2])))
  
  periods=tibble(Indicator=c(BRP$Indicator, AP$Indicator), dat=c(rep("BRP", nrow(BRP)), rep("AP", nrow(AP))))
  
  ## Identify larger variance for set order in var.test
  large_period=periods%>%
    dplyr::group_by(dat)%>%
    dplyr::summarise(var=var(Indicator))%>%
    arrange(desc(var))
  
  var1=periods%>%dplyr::filter(dat==large_period[[1,1]])
  var2=periods%>%dplyr::filter(dat!=large_period[[1,1]])
  
  
  variance_sign<-var.test(var1$Indicator, var2$Indicator, alternative="greater")[["p.value"]]
  
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
  
  ####   ----------- Trend analysis: linear regression on the last five years
  APTA<-data%>%dplyr::filter(year > (max(year)-6)) ## Last five years
  
  reg<-summary(lm(APTA$Indicator ~ APTA$year))
  
  if(reg[["coefficients"]][2,4] >= 0.05){
    print("No significant trend detected in TA")
    TA=0
  }else if(reg[["coefficients"]][2,4] <= 0.05 & reg[["coefficients"]][2,1]>0 ){
    print("Significant increasing trend detected in TA")
    TA=1
  }else if (reg[["coefficients"]][2,4] <= 0.05 & reg[["coefficients"]][2,1]<0 ){
    print("Significant decreasing trend detected in TA")
    TA=-1
  }
  
  
  #### Result table ####
  dev_expl=1-(summary(lm0)[["deviance"]]/summary(lm0)[["null.deviance"]])
  
  res=tibble(species=paste(species, gsa, sep="_"),formula= paste(best_model,collapse="+"), bp.year=bp[2], slope.1= slope(my.seg)$year[1,1], slope.2=slope(my.seg)$year[2,1], bp.fitted=dev_expl,mean.BRP= meanBRP, meanAP= meanAP, bp.test= regr.p , spearman.p= sprmn.p, rho.spearman=rho, trend_sign =reg[["coefficients"]][2,4], trend = reg[["coefficients"]][2,1])
  
  ### Final plot ####
  newdat=bind_cols(data, fit=my.seg$fitted.values)
  
  dataplot=bind_rows(newdat%>%dplyr::group_by(year)%>%dplyr::summarize(Indicator=mean(Indicator))%>%dplyr::mutate(source="observed"),
                     newdat%>%dplyr::group_by(year)%>%dplyr::summarize(Indicator=mean(fit))%>%dplyr::mutate(source="fitted"))
  
  
  if(nrow(dataplot%>%distinct(year)) > 20 ){
    spaz=5
  }else{
    spaz=2
  }
  
  baseplot=function(dat, clr, lbl){
    
    ggplot(data=dat %>%dplyr::filter(source=="observed"))+
      geom_point(aes(x=year, y=Indicator, color=source), size=3)+
      geom_vline(aes(xintercept=bp[2]), linetype="dashed")+
      geom_rect(aes(xmin = min(AP$year), xmax = max(AP$year), ymin =mean(AP$Val)-sd(AP$Val), ymax = mean(AP$Val)+sd(AP$Val)),alpha = 0.005, fill = clr)+ylab(indicator)+
      geom_segment(x=min(AP$year),xend=max(AP$year),y=mean(AP$Val), yend=mean(AP$Val), size=1, linetype = "dashed", color=clr)+
      geom_rect(aes(xmin = min(BRP$year), xmax = max(BRP$year), ymin =mean(BRP$Val)-sd(BRP$Val), ymax = mean(BRP$Val)+sd(BRP$Val)),alpha = 0.005, fill = "blue")+
      geom_segment(x=min(BRP$year),xend=max(BRP$year),y=mean(BRP$Val), yend=mean(BRP$Val), size=1, linetype = "dashed", color="blue")+
      theme_bw()+scale_x_continuous(breaks=seq(min(data$year),max(data$year),spaz))+
      theme_classic()+
      annotate("text", x=max(dat$year-spaz), y=max(dat$Indicator)-((max(dat$Indicator)-min(dat$Indicator))*0.02), label= lbl)+
      ggtitle(paste(sspp,"GSA",gsa,indicator,"segmented regression on", pop, "population"))+
      geom_smooth(data=dataplot%>%dplyr::filter(source=="fitted"),aes(year, Indicator, color=source), span = 0.5 ,stat="identity")
  }
  
  if(regr.p >= 0.05){
    
    pl=baseplot(dataplot, "blue", "BPA result: AP = RP")
    
  }else if(regr.p < 0.05){
    
    if(meanAP > meanBRP){
      
      pl=baseplot(dataplot, "green", "BPA result: AP > RP")
      
    }else{
      
      pl=baseplot(dataplot, "red", "BPA result: AP < RP")
      
    }
    
  }
  
  if(TA==0){
    pl+ 
      geom_smooth(method = "lm", se = TRUE, color="black", linetype="dashed", size=1.5, data=APTA,aes(year, Indicator))+
      annotate("text", x=max(dataplot$year-spaz), y=max(dataplot$Indicator)-((max(dataplot$Indicator)-min(dataplot$Indicator))*0.06), label= "No sign trend in recent years")
    
    }else if(TA==1){
        pl+
      geom_smooth(method = "lm", se = TRUE, color="green", linetype="dashed", size=1.5, data=APTA,aes(year, Indicator))+
        annotate("text", x=max(dataplot$year-spaz), y=max(dataplot$Indicator)-((max(dataplot$Indicator)-min(dataplot$Indicator))*0.06), label= "Significant increasing trend detected in recent years")
      
      }else if(TA==-1){
        pl+
        geom_smooth(method = "lm", se = TRUE, color="red", linetype="dashed", size=1.5, data=APTA,aes(year, Indicator))+
          annotate("text", x=max(dataplot$year-spaz), y=max(dataplot$Indicator)-((max(dataplot$Indicator)-min(dataplot$Indicator))*0.06), label= "Significant decreasing trend detected in recent years")
        
      }
    
  
  ### Save ####
  
  writexl::write_xlsx(res,paste0(dir_t, "/seg_reg", "/", species, "_GSA", gsa,indicator,pop, "_results.xlsx"))
  saveRDS(list(lm0, my.seg), file=paste0(dir_t, "/seg_reg", "/", species, "_GSA", gsa,indicator,pop, "_models.rds"))
  ggsave(file=paste0(dir_t, "/seg_reg", "/", species, "_GSA", gsa,indicator,pop, "_segmented_fit_vs_obs.png"), width=12)
  
#}
