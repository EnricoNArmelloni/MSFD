library(strucchange)
library(splitstackshape)
library(tidyverse)
library(readxl)
library(data.table)
rm(list=ls()) 


###--- Set parameters
sel_ind<-"BIOM"
species="EDT"
area="GSA18"
mydir_input="~/CNR/MSFD/github/MSFD/Data"
mydir_out="~/CNR/MSFD/github/MSFD/output" #### Directory for outputs
survey="MEDITS"


###-- Run
if(dir.exists(paste0(mydir_out,"/",species))==TRUE){
  NA
}else{
  dir.create(paste0(mydir_out,"/",species))
}

# Load data
if(survey=="MEDITS"){
  LFD<-read_csv(paste0("~/CNR/MSFD/github/MSFD/input/",area, "/", species,"/", species,"_", area, "__LFD.csv"))%>%dplyr::rename("Yr"="Year")
  L95haul= read_csv(paste0("~/CNR/MSFD/github/MSFD/input/",area, "/", species,"/", species,"_", area, "_L95_haul.csv"))%>%dplyr::rename("Year"="year")
  L95unc=read_csv(paste0("~/CNR/MSFD/github/MSFD/input/",area, "/", species,"/", species,"_", area,"_L95_plot.csv"))
  BIOMhaul= read_csv(paste0("~/CNR/MSFD/github/MSFD/input/",area, "/", species,"/", species,"_", area, "_BIOM_haul.csv"))#%>%dplyr::rename("Year"="year")
  BIOMunc=read_csv(paste0("~/CNR/MSFD/github/MSFD/input/",area, "/", species,"/", species,"_", area,"_BIOM_plot.csv"))%>%dplyr::rename("Year"="year", "Val"="total_biomass")
  
}else if (survey == "SOLEMON"){
  if(species =="MTS_simulated"){
    LFD<-read_delim("~/CNR/MSFD/Squilla case study/Numbers_at_lenght_MTS.csv",   ";", escape_double = FALSE, trim_ws = TRUE)%>%dplyr::filter(Seas==4)%>%dplyr::select(-Area, -Bio_Pattern,BirthSeas,-Settlement,-Platoon,-Morph,-Seas,-Era,-Sex,-BirthSeas, -Time)%>%dplyr::filter(`Beg/Mid`=="B")%>%dplyr::select(- `Beg/Mid`)%>%gather(value=Number, key=Lenght_class,- Yr)%>%group_by(Yr)%>%dplyr::mutate(Frequency=Number/sum(Number), Lenght=as.numeric(Lenght_class)*10)
  }else{
    LFD<-read_excel(paste0(mydir_input, "/LFDs/LFD_", species ,".xlsx"))%>%
      dplyr::mutate(Yr=as.numeric(substr(str_remove(Survey, "SOLEMON"), 1,4)))%>%
      dplyr::rename("Length"="LengthClass", "Frequency"="AbunIndex")
    #dplyr::mutate(Lenght_class=Lenght_class/10)
  }
}


pars<-read_excel(paste0(mydir_input, "/parameters.xlsx"))
Linf=pars$Linf[pars$species==species]
cutoff=1.1*((2/3)*Linf)

#if(is.na(pars$minyear[pars$species==species])==FALSE){
#  LFD<-LFD%>%dplyr::filter(Yr >= pars$minyear)
#}

#####---- Calculate Indicators

# L95
dati <- LFD %>%dplyr::mutate(Frequenza = round(Frequency*100))
dati <- expandRows(dati[complete.cases(dati),] , "Frequenza")

p<-c(0.95)
p_names <- map_chr(p, ~paste0( "perc", .x*100))
p_funs <- map(p, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% set_names(nm = p_names)

L95 <-dati%>%dplyr::group_by(Yr)%>%dplyr::summarize_at(vars(Length), funs(!!!p_funs))%>%
  dplyr::rename("Val"="perc95")


# Pmega
PMega<-LFD%>%dplyr::mutate(Frequency=ifelse(is.na(Frequency)==TRUE, 0, Frequency))%>%
  dplyr::mutate(ms=ifelse(Length >= cutoff, "Y", "N"))%>%
  group_by(Yr, ms)%>%dplyr::summarize(Parz=sum(Frequency*100))%>%
  group_by(Yr)%>%dplyr::mutate(Tot = sum(Parz), Val= Parz/Tot)%>%dplyr::filter(ms=="Y")%>%
  dplyr::select(Yr, Val)

# Biomass
Biom=BIOMunc%>%dplyr::select(Year,Val)%>%dplyr::rename("Yr"="Year")


##### Analysis
if(sel_ind=="L95"){
  Ind=L95
}else if(sel_ind=="PMega"){
  Ind=PMega
}else{
  Ind=Biom
}

if(is.na(pars$minyear[pars$species==species])==FALSE){
  Ind<-Ind%>%dplyr::filter(Yr >= pars$minyear[pars$species==species])
}

ggplot() + geom_line(data=Ind,aes(x=Yr,y=Val)) + geom_hline(data=Ind, aes(yintercept = mean(Ind$Val)), color="red")+ggtitle(paste0(sel_ind, " from observed population,", species," " ,area))+ylab("cm")+xlab("Year")
ggsave(paste0(mydir_out, "/", species, "/", sel_ind, species,  ".png"))

####   ----------- BP analysis
# Step 1: identify BP
dfpar<-Ind%>%dplyr::rename("Year"="Yr")
df<-ts(data=dfpar$Val)

bp.msfd <- breakpoints(df ~ 1)
if(exists("bp.msfd")==FALSE){
  print("BP analysis did not succeed")
  BPA=-2
}else{

png(filename=paste0(mydir_out, "/", species, "/", sel_ind, species,  "BP_selection.png"))
  plot(bp.msfd)
dev.off()


n_breaks<-length(breakpoints(bp.msfd)[["breakpoints"]])
bp<-breakpoints(bp.msfd)[[1]]

if((is.na(bp[1])==TRUE)==TRUE){
  print("BP analysis did not detect BPs")
  BPA=-2
  }else{
    
  

# Step 2:  test significance of difference between assessment period (Map) and reference period (RP). Under the "improve leading to recovery (IMP)" rationale, the the Map was reaching the assessment benchmark when it was significantly higher than the mean of the best period (BMRP)

brp<-c(1, bp, nrow(dfpar))
Lst<-list()

for(i in 1:n_breaks){
  RP<-tibble(Val=dfpar[brp[i]:brp[i+1],]$Val, Year=dfpar[brp[i]:brp[i+1],]$Year, RP=i)
  Lst[[i]]<-RP
}
BRP = rbindlist(Lst)%>%arrange(Year, RP)%>%dplyr::distinct(Year,.keep_all=TRUE)%>%group_by(RP)%>%summarize(mean=mean(Val))%>%dplyr::filter(mean==max(mean))
BRP = rbindlist(Lst)%>%arrange(Year, RP)%>%dplyr::distinct(Year,.keep_all=TRUE)%>%dplyr::filter(RP==BRP$RP)

AP=tibble(Val=dfpar[(bp[n_breaks]+1):nrow(dfpar),]$Val, Year=dfpar[(bp[n_breaks]+1):nrow(dfpar),]$Year) ## Assessment period is equal to last period of stability

variance_sign<-var.test(BRP$Val, AP$Val)[["p.value"]]

if(variance_sign >= 0.05){
  ttest<-t.test(BRP$Val , AP$Val,  var.equal = TRUE)
}else{
  ttest<-t.test(BRP$Val , AP$Val,  var.equal = FALSE)
}

meanRP<-ttest[["estimate"]][[1]]
meanAP<-ttest[["estimate"]][[2]]

if(ttest[["p.value"]] >= 0.05){
  print(paste("No significant differences detected in BP analysis, BPA assessment denote status as good as in the best period. RP =", round(meanRP,2),"; AP =", round(meanAP,2)))
  BPA= 0
}else if(ttest[["p.value"]] < 0.05 & meanAP > meanRP){
  print(paste("Significant differences detected in BP analysis, BPA assessment denote good status. RP =", round(meanRP,2),"; AP =", round(meanAP,2)))
  BPA= 1
}else if (ttest[["p.value"]] < 0.05 & meanAP < meanRP){
  print(paste("Significant differences detected in BP analysis, BPA assessment denote bad status. RP =", round(meanRP,2),"; AP =", round(meanAP,2)))
  BPA=-1
}
  } 
}
####   ----------- Trend analysis: linear regression on the last five years
APTA<-tibble(Val=dfpar[(nrow(dfpar)-4):nrow(dfpar),]$Val, Year=dfpar[(nrow(dfpar)-4):nrow(dfpar),]$Year) ## Last five years

reg<-summary(lm(APTA$Val ~ APTA$Year))

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

if(TA == 1 & BPA ==1 ){
  print("Final diagnosis: healthy stock")
  res="Final diagnosis: healthy stock"
} else if(TA == 1 & BPA !=1) {
  print("Final diagnosis: check management regimes")
  res="Final diagnosis: check management regimes"
}else if(TA == 0 & BPA ==0) {
  print("Final diagnosis: stock status stable")
  res="Final diagnosis: stock status stable"
  
}else if(TA != 1 & BPA ==-1) {
  print("Final diagnosis: decrease pressure")
  res="Final diagnosis: decrease pressure"
}else{
  print("Final diagnosis: case not considered")
  res="Final diagnosis: case not considered"
}

if(BPA == -2){
  print("No Plot to show")
}else{

### final plot

dfpar_mean=tibble(mean_95perc=mean(dfpar$Val), Percentile="Mean 0.95")

if(nrow(dfpar) > 20 ){
  spaz=5
}else{
  spaz=2
}

plotfin<-dfpar%>%ggplot(aes(x=Year, y=Val))+
  geom_rect(aes(xmin = min(BRP$Year), xmax = max(BRP$Year), ymin =mean(BRP$Val)-sd(BRP$Val), ymax = mean(BRP$Val)+sd(BRP$Val)),alpha = 0.01, fill = "darkgrey")+
  geom_segment(x=min(BRP$Year),xend=max(BRP$Year),y=mean(BRP$Val), yend=mean(BRP$Val), size=1, linetype = "dashed", color="blue")

if(BPA == 0){
  
  plotfin=plotfin+
    geom_rect(aes(xmin = min(AP$Year), xmax = max(AP$Year), ymin =mean(AP$Val)-sd(AP$Val), ymax = mean(AP$Val)+sd(AP$Val)),alpha = 0.01, fill = "blue")+
    geom_segment(x=min(AP$Year),xend=max(AP$Year),y=mean(AP$Val), yend=mean(AP$Val), size=1, linetype = "dashed", color="blue")+
    geom_line()+ 
    #geom_hline(data=dfpar_mean, aes(yintercept = mean_95perc, color=Percentile), size=0.02, linetype = "dashed")+ 
    geom_vline(xintercept = dfpar[bp,]$Year, linetype="dashed", size=0.01)+
    scale_x_continuous(breaks=seq(min(dfpar$Year),max(dfpar$Year),spaz))+
    theme_classic()+
    annotate("text", x=max(dfpar$Year-spaz), y=max(dfpar$Val)-((max(dfpar$Val)-min(dfpar$Val))*0.02), label= "BPA result: AP = RP") 
  
 }else if(BPA > 0){
  plotfin=plotfin+
    geom_rect(aes(xmin = min(AP$Year), xmax = max(AP$Year), ymin =mean(AP$Val)-sd(AP$Val), ymax = mean(AP$Val)+sd(AP$Val)),alpha = 0.01, fill = "green")+
    geom_segment(x=min(AP$Year),xend=max(AP$Year),y=mean(AP$Val), yend=mean(AP$Val), size=1, linetype = "dashed", color="green")+
    geom_line()+ 
    #geom_hline(data=dfpar_mean, aes(yintercept = mean_95perc, color=Percentile), size=0.02, linetype = "dashed")+ 
    geom_vline(xintercept = dfpar[bp,]$Year, linetype="dashed", size=0.01)+
    scale_x_continuous(breaks=seq(min(dfpar$Year),max(dfpar$Year),spaz))+
    theme_classic()+
    annotate("text", x=max(dfpar$Year-spaz),  y=max(dfpar$Val)-((max(dfpar$Val)-min(dfpar$Val))*0.02), label= "BPA result: AP > RP") 
  
 }else if(BPA < 0){
  plotfin=plotfin+
    geom_rect(aes(xmin = min(AP$Year), xmax = max(AP$Year), ymin =mean(AP$Val)-sd(AP$Val), ymax = mean(AP$Val)+sd(AP$Val)),alpha = 0.01, fill = "red")+
    geom_segment(x=min(AP$Year),xend=max(AP$Year),y=mean(AP$Val), yend=mean(AP$Val), size=1, linetype = "dashed", color="red")+
    geom_line()+ 
    #geom_hline(data=dfpar_mean, aes(yintercept = mean_95perc, color=Percentile), size=0.02, linetype = "dashed")+ 
    geom_vline(xintercept = dfpar[bp,]$Year, linetype="dashed", size=0.01)+
    scale_x_continuous(breaks=seq(min(dfpar$Year),max(dfpar$Year),spaz))+
    theme_classic()+
    annotate("text", x=max(dfpar$Year-spaz), y=max(dfpar$Val)-((max(dfpar$Val)-min(dfpar$Val))*0.02), label= "BPA result: AP < RP") 
 }


if(TA==0) {
 pf=plotfin+geom_smooth(method = "lm", se = TRUE, color="black", linetype="dashed", size=1.5, data=dfpar%>%dplyr::filter(Year >= max(Year)-4))+
   annotate("text", x=max(dfpar$Year-spaz), y=max(dfpar$Val)-((max(dfpar$Val)-min(dfpar$Val))*0.07), label= "TA result: no sign trend in recent years") 
}else if(TA == 1){
 pf=plotfin+geom_smooth(method = "lm", se = TRUE, color="green", size=1.5, data=dfpar%>%dplyr::filter(Year>= max(Year)-4))+
   annotate("text", x=max(dfpar$Year-spaz), y=max(dfpar$Val)-((max(dfpar$Val)-min(dfpar$Val))*0.07), label= "TA result: positive trend in recent years") 
}else if(TA == -1){
  pf=plotfin+geom_smooth(method = "lm", se = TRUE, color="red", size=1.5, data=dfpar%>%dplyr::filter(Year >= max(Year)-4))+
    annotate("text", x=max(dfpar$Year-spaz), y=max(dfpar$Val)-((max(dfpar$Val)-min(dfpar$Val))*0.07), label= "TA result: negative trend in recent years") 
}
pf+ggtitle(paste(species, area, sel_ind))+ annotate("text", x=max(dfpar$Year-spaz), y=max(dfpar$Val)-((max(dfpar$Val)-min(dfpar$Val))*0.12), label= res) 
ggsave(paste0(mydir_out, "/", species, "/", sel_ind, species,  "summary.png"), width=200, units="mm")  
}

#}


