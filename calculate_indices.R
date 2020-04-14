library(splitstackshape)
library(tidyverse)
library(readxl)

### set parameters
GSA<-"17"
species<- "Eledone cirrhosa"
sp_code<-"EOI"



#### Imput Data
Biomass <- read_csv("~/CNR/MSFD/data/MFSD_ANALISI MEDITS SOLEMON GSA 17/Output/GSA17/ITA_HRV_SVN/EOI/ELEDCIR_GSA17_ITA_HRV_SVN_BIOMASS.csv") ## Biomass trend as for MEDITS script
Catches <- read_excel("~/CNR/MSFD/data/Landings_CNR.xlsx", sheet = "Landings")%>%dplyr::filter(SPECIES == "EOI", AREA == "17")%>%dplyr::group_by(ANNO)%>%dplyr::summarize(Landing = sum(WEIGHT))
LFD<- read_csv("~/CNR/MSFD/data/MFSD_ANALISI MEDITS SOLEMON GSA 17/Output/GSA17/ITA_HRV_SVN/EOI/ELEDCIR_GSA_17_ITA_HRV_SVN_LFDTOT_10-800m.csv")
###

### Out creation
dir.create(file.path(paste0(sp_code, GSA)))


##### D3 - C1: catches/biomass
B<-Biomass %>%dplyr::filter(year %in% seq(min(Catches$ANNO), max(Catches$ANNO)))%>%dplyr::select(total_biomass, year)
C<-Catches%>%dplyr::rename("year"="ANNO")
df<-inner_join(B,C, by="year")%>%dplyr::select(year, Landing, total_biomass)%>%dplyr::mutate(Ratio= Landing/total_biomass)
qnt<-tibble(percent=quantile(df$Ratio, c(0.33, 0.66)), Percentile=c("0.33", "0.66"))
ggplot()+ geom_line(data=df, aes(x=year, y=Ratio))+ geom_hline(data=qnt, aes(yintercept=percent, color=Percentile))+ ylab("Commercial landings/ Biomass Index")+ ggtitle("D3C1", subtitle=paste(species, GSA, sep=" "))
ggsave(paste0(".","/",sp_code, GSA,"/", "D3C1.png"))

##### D3 - C2 biomass
qnt<-tibble(percent=quantile(Biomass$total_biomass, c(0.33, 0.66)), Percentile=c("0.33", "0.66"))
ggplot()+ geom_line(data=Biomass, aes(x=year, y= total_biomass))+ geom_hline(data=qnt, aes(yintercept=percent, color=Percentile))+ ylab("Biomass (kg/km2)")+ ggtitle("D3C2", paste(species, GSA, sep=" "))
ggsave(paste0(".","/",sp_code, GSA,"/", "D3C2.png"))

#### D3 - C3 95%? LFD from survey
datiraw <- LFD
LFD$Frequency <- (LFD$Frequency*10)
datiraw$Frequenza <- (round(LFD$Frequency))
dati_compl <- datiraw[complete.cases(datiraw),] 
dati <- expandRows(dati_compl, "Frequenza")

p<-c(0.05, 0.95)
p_names <- map_chr(p, ~paste0( "perc", .x*100))
p_funs <- map(p, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% set_names(nm = p_names)

percentiledf <-dati%>%dplyr::group_by(Year)%>%dplyr::summarize_at(vars(Length), funs(!!!p_funs))
percentiledf_mean <-na.omit(percentiledf)%>%dplyr::summarise(mean_95perc=mean(perc95))

ggplot() + geom_line(data=percentiledf,aes(x=Year,y=perc95)) + geom_hline(data=percentiledf_mean, aes(yintercept = mean_95perc), color="blue")+ ggtitle("D3C3", paste(species, GSA, sep=" "))+ ylab("0.95 percentile lenght (mm)")
ggsave(paste0(".","/",sp_code, GSA,"/", "D3C3.png"))

