#### The following code is built on

# ---
# title: "Computing MEDITS indexes and LFD"
# author: "Alessandro Mannini"
# date: "October 4th, 2018"
# ---
#    Copyright (C) <2014>  <Alessandro Mannini>
#    excluding code about read in stratification table of the survey: Copyright Tristan Rouyer
#    excluding code about bubble plot by year of density and biomass by haul: Copyright Isabella Bitetto and Maria Teresa Facchini
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/
# 

remove(list=ls())

library(matlab)
tic()
require(doBy)  # summaryBy
require(lattice)#cool graphics
require(gclus)#correlation output - Scatterplot Matrices
library(reshape)
require(data.table)
require(ggmap)
require(mapdata)
require(maps)
library(stringr)
library(tidyverse)
library(OpenStreetMap)
library(splitstackshape)

#### ----- Settings: set lines 42-64 basing on your data  ----- ####

input_dir="~/CNR/MSFD/github/release/data" ### folder where are stored TA TB TC asfis list

output_dir="~/CNR/MSFD/github/release" # folder for outputs

# need to exclude strata? If yes, activate line 50 and type the strata to include
selection=NA ## let this parameter unchanged for the present version
#selection=c("B", "C")

Lmat = 230 #Lmat in mm

Linf= 380

Lopt= (2/3)*Linf

survey="medits"

compile_AMSY="N" # parameter for future usage

sspp= "FLE" # three alpha code of your species

areacode=c("17") # use combine for many GSA

countrycode=c("ITA")# Use uppercase and use combine for pull together two or more countries

state="ITA" #assign country code to output file name (i.e ITA_HRV_ for 2 countries)

gsa="17" #assign gsa code to output file name (i.e. 9_11_ for two GSAs) 

lfstep=50 # set as 1 if you are working with crustaceans and 50 for fish and cephalopods

maxratiosampling=5 # Set value according to what do you think should be a reasonable max value sampling factor (e.g. maximum value of subsample apply to the total catch in a haul)

#### ---- Code start ------ ####

### Folder creation #####
setwd(input_dir)
dir_t=paste0(output_dir,"/","GSA",gsa,sep="")
dir.create(file.path(dir_t,sspp), recursive = T)

#### Read data ####
GSAtable=read.csv("gsa_coordinates.csv") # use to bubble plot legend position
alpha_code=sspp
mdts=fread("Sp_Medits.csv") ; mdts=mdts[,c(2,3,7)] ; colnames(mdts)[2]="Scientific_name"
asfis=fread("ASFIS_2017.csv"); asfis=asfis[,c(3,4)]
sp_list=merge(mdts,asfis,by="Scientific_name",all = T)
sp_list=subset(sp_list[!is.na(sp_list$code),])### NA values dismiss
sp_list=subset(sp_list[!is.na(sp_list$`3A_CODE`),])### NA values dismiss
sp_list$gen=str_sub(sp_list$code,start=1,end=4)
sp_list$spec=str_sub(sp_list$code,start=5,end=7)
gen=unique(sp_list$gen[which(sp_list$`3A_CODE`==alpha_code)])
spec=unique(sp_list$spec[which(sp_list$`3A_CODE`==alpha_code)])

# Extract data ####
TAn <- fread(paste0(survey, "_ta.csv"))
TAn <- as.data.frame(subset(TAn, (area %in% areacode) & (country %in% countrycode))) %>%dplyr::filter(validity == "V")
droplevels(TAn)

TBn <- fread(paste0(survey, "_tb.csv"))
TBn <- as.data.frame(subset(TBn, (area %in% areacode) & (country %in% countrycode) & (genus %in% gen) & (species %in% spec)))
droplevels(TBn)

TCn <- fread(paste0(survey, "_tc.csv"))
TCn <- as.data.frame(subset(TCn, (area %in% areacode) & (country %in% countrycode) & (genus %in% gen) & (species %in% spec)))
droplevels(TCn)



if(survey=="medits"){
  
  ### FIXING MISTAKES ###
  # Please activate code needed in your own stock #
  # There are mistakes in length class measure (e.g. HRV 2016 or ITA17 2017 for NEP)
  #TCn$length_class=ifelse(TCn$country=="HRV" & TCn$year==2016 & alpha_code=="NEP",TCn$length_class/10,TCn$length_class)
  #TCn$length_class=ifelse(TCn$country=="HRV" & alpha_code=="NEP" & TCn$length_class>100,TCn$length_class/10,TCn$length_class)
  #TCn$length_class=ifelse(TCn$country=="ITA" & alpha_code=="NEP" & TCn$length_class>100,TCn$length_class/10,TCn$length_class)
  #TCn$length_class=ifelse(TCn$country=="ESP" & alpha_code=="NEP" & TCn$length_class>100,TCn$length_class/10,TCn$length_class)
  TCn$length_class=ifelse(TCn$country=="ITA" & alpha_code=="HKE" & TCn$length_class>900,TCn$length_class/10,TCn$length_class)
  #TCn$length_class=ifelse(TCn$country=="ITA" & alpha_code=="DPS" & TCn$length_class>50,TCn$length_class/10,TCn$length_class)
  #TCn$length_class=ifelse(TCn$country=="ESP" & alpha_code=="ARA" & TCn$area==5 & TCn$length_class>80,TCn$length_class/10,TCn$length_class)
  #TCn$length_class=ifelse(TCn$country=="ESP" & alpha_code=="ARA" & TCn$area==6 & TCn$length_class>80,TCn$length_class/10,TCn$length_class)
  TCn$length_class=ifelse(TCn$country=="ITA" & alpha_code=="SYC" & TCn$length_class>800,TCn$length_class/10,TCn$length_class)
  # There are mistakes in length class value in MTS GSA 17 ITA Years 2012:2016 (converting TL in CL based on values find in paper)
  # TCn$length_class=ifelse(TCn$country=="ITA" & TCn$year==2012 & TCn$length_class>50,as.integer(exp((log(TCn$length_class)-log(7.3949))/0.8693),round(0)),TCn$length_class)
  # TCn$length_class=ifelse(TCn$country=="ITA" & TCn$year==2013 & TCn$length_class>50,as.integer(exp((log(TCn$length_class)-log(7.3949))/0.8693)),TCn$length_class)
  # TCn$length_class=ifelse(TCn$country=="ITA" & TCn$year==2014 & TCn$length_class>50,as.integer(exp((log(TCn$length_class)-log(7.3949))/0.8693)),TCn$length_class)
  # TCn$length_class=ifelse(TCn$country=="ITA" & TCn$year==2015 & TCn$length_class>50,as.integer(exp((log(TCn$length_class)-log(7.3949))/0.8693)),TCn$length_class)
  # TCn$length_class=ifelse(TCn$country=="ITA" & TCn$year==2016 & TCn$length_class>50,as.integer(exp((log(TCn$length_class)-log(7.3949))/0.8693)),TCn$length_class)
  
}




# read in stratification table of the survey Copyrigth Tristan Rouyer
stratification_scheme <- fread(paste0(survey,"_Strata.csv"))
stratum<-stratification_scheme[stratification_scheme$AREA %in% areacode,]#extraction based on gSA
stratum<-stratum[stratum$COUNTRY %in% countrycode,]# extraction based on country

## prepare the TA file for the next elaborations##
meandepth=(TAn$shooting_depth+TAn$hauling_depth)/2
sqkm=TAn$wing_opening/10000000*TAn$distance
id2=paste(TAn$country,TAn$area,TAn$year,TAn$haul_number,sep="")
TAn["strata"]=NA
TAn=cbind(TAn,meandepth,sqkm,id2)
TAn$strata[]=TAn$nstrate[]


if(survey=="medits"){
  
  ##Assign strata code based on the meandepth value and not on "codestrata"
  for (i in 1:length(TAn$strata))
    if(TAn$meandepth[i]>0 & TAn$meandepth[i] < 51){TAn$strata[i]="A"}else{
      if(TAn$meandepth[i]>=51 & TAn$meandepth[i] < 101){TAn$strata[i]="B"}else{
        if(TAn$meandepth[i]>=101 & TAn$meandepth[i] < 201){TAn$strata[i]="C"}else{
          if(TAn$meandepth[i]>=201 & TAn$meandepth[i] <501){TAn$strata[i]="D"}else{
            TAn$strata[i]="E"}}}}
  
} else if(survey=="solemon"){
  
  ##Assign strata code based on the meandepth value and not on "codestrata"
  for (i in 1:length(TAn$strata))
    if(TAn$meandepth[i]>0 & TAn$meandepth[i] < 31){TAn$strata[i]="A"}else{
      if(TAn$meandepth[i]>=31 & TAn$meandepth[i] < 51){TAn$strata[i]="B"}else{
        if(TAn$meandepth[i]>=51 & TAn$meandepth[i] < 121){TAn$strata[i]="C"}else{
          if(TAn$meandepth[i]>=121 & TAn$meandepth[i] <201){TAn$strata[i]="D"}else{
            TAn$strata[i]="E"}}}}
  
}


#prepare TB for next elaborations#
id2=paste(TBn$country,TBn$area,TBn$year,TBn$haul_number,sep="")
TBn=cbind(TBn,id2)
#create a new database merging TATB
TATBn=merge(TAn,TBn,by="id2",all=T)
TATBn$W_sqkm=TATBn$ptot/TATBn$sqkm/1000
TATBn$N_sqkm=TATBn$nbtot/TATBn$sqkm


# dir_t=getwd()
# dir.create(file.path(dir_t,"output"))
# dir.create(file.path(dir_t,paste0("output","/",state,"_",gsa,"_",gen,spec)))

####### Coordinates formatting
TA_lon=formatC(as.numeric(TATBn$shooting_longitude),width=7,format='f',digits=2,flag='0')


TA_gr_lon=substr(TA_lon, 1, 2)
TA_gr_lon=as.integer(TA_gr_lon)
which(is.na(TA_gr_lon))
TA_gr_lon[is.na(TA_gr_lon)] <- 0

TA_mi_lon=substr(TA_lon, 3, 4)
TA_mi_lon=as.integer(TA_mi_lon)
which(is.na(TA_mi_lon))
TA_mi_lon[is.na(TA_mi_lon)] <- 0

TA_se_lon=substr(TA_lon, 6, 7)
TA_se_lon=as.integer(TA_se_lon)
which(is.na(TA_se_lon))
TA_se_lon[is.na(TA_se_lon)] <- 0

TA_lon_fin=TA_gr_lon+(TA_mi_lon/60)+(TA_se_lon/6000)# Attention are in decimal of minutes. Need to find conversion factor maybe 6000
which(is.na(TA_lon_fin))
TATBn$TA_lon_fin=TA_lon_fin


TA_lat=formatC(as.numeric(TATBn$shooting_latitude),width=7,format='f',digits=2,flag='0')

TA_gr_lat=substr(TA_lat, 1, 2)
TA_gr_lat=as.integer(TA_gr_lat)
which(is.na(TA_gr_lat))
TA_gr_lat[is.na(TA_gr_lat)] <- 0

TA_mi_lat=substr(TA_lat, 3, 4)
TA_mi_lat=as.integer(TA_mi_lat)
which(is.na(TA_mi_lat))
TA_mi_lat[is.na(TA_mi_lat)] <- 0

TA_se_lat=substr(TA_lat, 6, 7)
TA_se_lat=as.integer(TA_se_lat)
which(is.na(TA_se_lat))
TA_se_lat[is.na(TA_se_lat)] <- 0

TA_lat_fin=TA_gr_lat+(TA_mi_lat/60)+(TA_se_lat/6000)# Attention are in decimal of minutes. Need to find conversion factor maybe 6000
which(is.na(TA_lat_fin))
TATBn$TA_lat_fin=TA_lat_fin
TATBn$TA_lon_fin=TA_lon_fin

TATBn$TA_lon_fin=ifelse(TATBn$shooting_quadrant=="7",-TATBn$TA_lon_fin,TATBn$TA_lon_fin)
TA_lon_fin=TATBn$TA_lon_fin
##############?
getwd()
setwd(file.path(dir_t,alpha_code))


# Exploring hauls time series ####
hauls_table=aggregate(TAn$haul_number,list(TAn$year,TAn$strata),length)
names(hauls_table)=c("year","stratum","hauls")
hauls_table2=hauls_table
hauls_table_tot=aggregate(hauls~year, hauls_table, sum)
hauls_table_tot$stratum=rep("total",nrow(hauls_table_tot))
hauls_table_tot=hauls_table_tot[,c(1,3,2)]
hauls_table=rbind(hauls_table,hauls_table_tot)

# Exploring vessels ####
vessels=TAn%>%dplyr::distinct(year, vessel,strata)

ggplot()+
  geom_point(data=vessels, aes(x=year, y=vessel))

ggsave(filename=paste(state,"_","GSA_",gsa,"exploring_vessels.png"),width = 10, height = 8, dpi = 150, units = "in")


# Exploring swept area ####
exploring_sweptarea=ggplot(TAn, aes(x=meandepth, y=sqkm,color=strata)) +
  geom_point(cex=0.9) + facet_wrap(~strata,scales ="free") +
  theme(legend.position="none") + geom_smooth()
exploring_sweptarea
ggsave(filename=paste(state,"_","GSA_",gsa,"exploring_sweptarea.png"),width = 10, height = 8, dpi = 150, units = "in", plot=exploring_sweptarea)

# Exploring survey period ####
TA_date=TAn%>%dplyr::select(year, month, day,strata)
TA_date$data=as.Date(paste0(TA_date$year,"-",TA_date$month,"-",TA_date$day))
invisible(TA_date[order(TA_date$data),])
TA_date$dayofyear=as.numeric(format(TA_date$data, "%j"))
TA_date=TA_date[!duplicated(TA_date$data), ]

survey_period=ggplot(TA_date, aes(x=month, y=dayofyear))+
  geom_rect(data=NULL,aes(xmin=min(month),xmax=max(month),ymin=0,ymax=90,
          fill="q1_winter"))+
  geom_rect(data=NULL,aes(xmin=min(month),xmax=max(month),ymin=91,ymax=180,
            fill="q2_spring"))+
  geom_rect(data=NULL,aes(xmin=min(month),xmax=max(month),ymin=181,ymax=270,
            fill="q3_summer"))+
  geom_rect(data=NULL,aes(xmin=min(month),xmax=max(month),ymin=271,ymax=365,
            fill="q4_fall"))+
  scale_fill_manual('Season_quarter',
                    values = c("lightblue","lightgreen","lightyellow","lightsalmon"))+
  geom_point(size=3,colour="red")+
  facet_wrap(~TA_date$year)+
  ggtitle(paste0("GSA",gsa,state)) +
  xlab("Month") + ylab("Day_of_year")

survey_period
ggsave(filename=paste(state,"_","GSA_",gsa,"survey_period.png"),width = 10, height = 8, dpi = 150, units = "in", plot=survey_period)

#### Create variables df
side_vars=hauls_table%>%
  dplyr::filter(stratum!="total")%>%
  dplyr::group_by(year)%>%
  dplyr::summarise(nhauls=sum(hauls))%>%
  left_join(TA_date%>%
              dplyr::group_by(year)%>%
              dplyr::summarise(mean_doy=round(mean(dayofyear))), by="year")

if(is.na(selection)==F){
  
  side_vars_sel=hauls_table%>%
    dplyr::filter(stratum %in% selection)%>%
    dplyr::group_by(year)%>%
    dplyr::summarise(nhauls=sum(hauls))%>%
    left_join(TA_date%>%
                dplyr::filter(strata %in% selection)%>%
                dplyr::group_by(year)%>%
                dplyr::summarise(mean_doy=round(mean(dayofyear))), by="year")
}

# Species occurence ########
tothaulyear=tapply(TATBn$haul_number.x,TATBn$year.x,length)
TATBn['pos']=1
TATBn['neg']=0
TATBn['occurence']=NA
TATBn$nbtot[is.na(TATBn$nbtot)]=0

# Abundance Indexes by strata####
strata=aggregate(TAn$sqkm,list(TAn$year,TAn$strata),sum)
names(strata)=c("year","stratum","sweptarea")
weigthbystrata=aggregate(TATBn$ptot,list(TATBn$year.y,TATBn$strata),sum)
names(weigthbystrata)=c("year","stratum","weigth")
numberbystrata=aggregate(TATBn$nbtot,list(TATBn$year.y,TATBn$strata),sum)
names(numberbystrata)=c("year","stratum","number")
biomass=merge(strata,weigthbystrata,by=c("year","stratum"))
biomass$kg_km2=((biomass$weigth/biomass$sweptarea)/1000)
density=merge(strata,numberbystrata,by=c("year","stratum"))
density$n_km2=((density$number/density$sweptarea))


ggplot(data=biomass)+
  geom_line(aes(x=year, y=kg_km2))+facet_wrap(~stratum)

ggsave(filename=paste(state,"_","GSA_",gsa,"biom_by_strata.png"),width = 10, height = 8, dpi = 150, units = "in")

# Area and weights of each stratum####
area.str=stratum %>%
  dplyr::group_by(CODESTRATA)%>%
  dplyr::summarise(area=sum(AREASTRATA))%>%
  dplyr::mutate(W=area/sum(area))%>%
  dplyr::rename("stratum"="CODESTRATA")

if(is.na(selection)==F){
  area.sel=area.str%>%
    dplyr::filter(stratum %in% selection)%>%
    dplyr::mutate(W=area/sum(area))
  
}
# Indexes by macrostrata (shelf, slope and total)####
# BIOMASS ####
#Total#
Biomass_tot=biomass%>%
  left_join(area.str, by="stratum")%>%
  dplyr::mutate(biomStrata=kg_km2*W )%>%
  dplyr::group_by(year)%>%
  dplyr::summarise(total_biomass=sum(biomStrata))

# selected strata
if(is.na(selection)==F){
  
  Biomass_sel=biomass%>%
    dplyr::filter(stratum %in% selection)%>%
    left_join(area.sel, by="stratum")%>%
    dplyr::mutate(biomStrata=kg_km2*W )%>%
    dplyr::group_by(year)%>%
    dplyr::summarise(sel_biomass=sum(biomStrata))
}


# Correction factor f ###
SW=strata%>%
  left_join(area.str, by="stratum")%>%
  dplyr::mutate(f=sweptarea/as.integer(area))

if(is.na(selection)==F){
  
  SW_sel=strata%>%
    dplyr::filter(stratum%in%selection)%>%
    left_join(area.sel, by="stratum")%>%
    dplyr::mutate(f=sweptarea/as.integer(area))
}



# Elaboration ###
TATBn$id3=paste0(TATBn$year.x,TATBn$strata,sep="")
biomass$id3=paste0(biomass$year,biomass$stratum,sep="")
hauls_table2$id3=paste0(hauls_table2$year,hauls_table2$stratum,sep="")
TATBBI=merge(TATBn,hauls_table2%>%dplyr::select(-year, -stratum),by="id3",all=TRUE)
TATBBI=merge(TATBBI,biomass,by="id3",all.y=TRUE)



TATBBI$var1=TATBBI$sqkm*(TATBBI$W_sqkm-TATBBI$kg_km2)^2
var1=aggregate(TATBBI$var1,list(TATBBI$year,TATBBI$strata),sum,na.rm=TRUE)
var1$id3=paste0(var1$Group.1,var1$Group.2,sep="")

var2a=merge(var1,hauls_table2,by="id3",all=T)
var2a$variance=(1/(var2a$hauls-1))*var2a$x
var2a$devst=sqrt(var2a$variance)


biom_var=var2a%>%
  left_join(SW, by=c("year", "stratum"))%>%
  group_by(stratum)%>%
  dplyr::mutate(Vt=(variance*W^2)/(sweptarea*(1-f)))%>%
  ungroup()%>%
  dplyr::group_by(year)%>%
  dplyr::summarize(var_tot=sum(Vt))%>%
  dplyr::mutate(stdev_tot=sqrt(var_tot))

BIOMASS=Biomass_tot%>%left_join(biom_var, by="year")
BIOMASS[is.na(BIOMASS)]=0


plo_biom=ggplot(BIOMASS, aes(x=year, y=total_biomass)) + 
  geom_line(linetype = "solid",size=1.25,color="tomato") +
  geom_point(color="tomato")+geom_hline(aes(yintercept = mean(total_biomass)), color="tomato")+
  geom_errorbar(aes(ymin=total_biomass-stdev_tot, ymax=total_biomass+stdev_tot), width=.2,col="2",
                position=position_dodge(0.05),linetype = "dashed",size=0.75)+
  ylab("kg/km2")+
  ggtitle(paste0(gen,spec,"_","GSA",gsa,"_",state,"Total_biomass"))

plo_biom

if(is.na(selection)==F){
  selection_txt=paste0(selection, collapse=" , ")
  
  biom_var_sel=var2a%>%
    dplyr::filter(stratum %in% selection)%>%
    left_join(SW_sel, by=c("year", "stratum"))%>%
    dplyr::group_by(stratum)%>%
    dplyr::mutate(Vt=(variance*W^2)/(sweptarea*(1-f)))%>%
    ungroup()%>%
    dplyr::group_by(year)%>%
    dplyr::summarize(var_sel=sum(Vt))%>%
    dplyr::mutate(stdev_sel=sqrt(var_sel))
  
  BIOMASS=BIOMASS %>%
    dplyr::left_join(Biomass_sel, by="year")%>%
    left_join(biom_var_sel, by="year")
  
  dataplot=BIOMASS%>%
    gather(key="set", value="biomass", total_biomass, sel_biomass)%>%
    gather(key="set_2", value="st_dev", stdev_tot, stdev_sel)%>%
    dplyr::filter(set=="total_biomass" & set_2=="stdev_tot"|set=="sel_biomass" & set_2=="stdev_sel")
  
  plo_biom=ggplot(dataplot, aes(x=year, y=biomass, color=set)) + 
    geom_line(linetype = "solid",size=1.25) +
    geom_point()+
    geom_hline(data=dataplot%>%dplyr::filter(set=="total_biomass"),aes(yintercept = mean(biomass), color=set))+
    geom_hline(data=dataplot%>%dplyr::filter(set=="sel_biomass"),aes(yintercept = mean(biomass), color=set))+
    geom_errorbar(aes(ymin=biomass-st_dev, ymax=biomass+st_dev), width=.2,
                  position=position_dodge(0.05),linetype = "dashed",size=0.75)+
    ylab("kg/km2")+
    ggtitle(paste0(gen,spec,"_","GSA",gsa,"_",state," Biomass total and selected strata: ", selection_txt))
  
  plo_biom

}

write.csv(BIOMASS,file=paste0(alpha_code,"_GSA",gsa,"_","BIOM.csv"))

ggsave(filename=paste(state,"_","GSA_",gsa,"Total_biomass.png"),width = 10, height = 8, dpi = 150, units = "in")

# Compile AMSY info
if(compile_AMSY=="Y"){
  
  dir.create(paste0(dir_t,"/", alpha_code, "/AMSY"))
  
  format_AMSY=BIOMASS %>%dplyr::select(Year=year, CPUE=total_biomass)%>%
    dplyr::mutate(Stock=alpha_code, Catch=NA)%>%
    dplyr::select(Stock, Year, Catch, CPUE)
  
  newf=data.frame(Year=seq(min(format_AMSY$Year), max(format_AMSY$Year)))
  
  format_AMSY=newf%>%left_join(format_AMSY, by="Year")%>%
    dplyr::mutate(Stock=alpha_code)
  
  for(i in 1:nrow(format_AMSY)){
    if(is.na(format_AMSY[i,4])==TRUE){
      
      if(is.na(format_AMSY[i+1,4])==FALSE){
        
        format_AMSY[i,4]=(format_AMSY[i+1,4] + format_AMSY[i-1,4])/2
        
      }else if(is.na(format_AMSY[i+1,4])==FALSE) {
        print("check gaps in biomass")
        
      }
    }
  }
  
  
  
 write.csv(format_AMSY%>%
             dplyr::select(Stock, Year, Catch, CPUE), file=paste0(dir_t,"/", alpha_code, "/AMSY/MSFD_Stocks_CPUE_medits",alpha_code, ".csv"), row.names = F)
  
}


# Standardized LFDs by km2 ####
raise=TCn$pfrac/TCn$pechan
nblonraise=TCn$nblon*raise
TCn=cbind(TCn,raise,nblonraise)
# TCtoCheck <- TCn[which(TCn$raise > maxratiosampling),]
# write.table(TCtoCheck,file=paste("TCtoCheck_",state,gsa,gen,spec,".csv"),sep=";",row.names=F)
id2=paste(TCn$country,TCn$area,TCn$year,TCn$haul_number,sep="")
TCn=cbind(TCn,id2)
TATCn=merge(TAn,TCn,by=c("id2","year","haul_number", "area"),all=T)

# Checking TC raising factor and compare TC and TB total weight and number ###
TCtoCheck <- TCn[which(TCn$raise > maxratiosampling),]
#write.table(TCtoCheck,file=paste("TCtoCheck_",state,gsa,gen,spec,".csv"),sep=";",row.names=F)

tempwgB <- TBn%>% 
  dplyr::group_by(country,year,area, haul_number)%>% 
  dplyr::summarize(totwgB=sum(ptot))

tempnbB <- TBn%>% 
  dplyr::group_by(country,year,area, haul_number)%>% 
  dplyr::summarize(totnbB=sum(nbtot))

tempTB <- merge(tempwgB,tempnbB,by=c("country","year","area", "haul_number"),all=T)

tempwgC <- TCn%>% dplyr::group_by(country,area,year,haul_number)%>% dplyr::summarize(totwgC=mean(pfrac))
tempnbC <- TCn%>% dplyr::group_by(country,area,year,haul_number)%>% dplyr::summarize(totnbC=sum(nblonraise))
tempTC <- merge(tempwgC,tempnbC,by=c("country","area","year","haul_number"),all=T)

TBTCcheck <- merge(tempTB,tempTC,by=c("country", "area", "year", "haul_number"))
TBTCcheck$wgratio <- TBTCcheck$totwgB/TBTCcheck$totwgC
TBTCcheck$nbratio <- TBTCcheck$totnbB/TBTCcheck$totnbC

TBTCcheck[which(TBTCcheck$wgratio > 1 |TBTCcheck$nbratio > 1),]
#write.csv(TBTCcheck,file=paste0(gen,spec,"_GSA_",gsa,"_",state,"_","TBTCcheck.csv",sep=""),row.names = F)
#write.csv(TBTCcheck[which(TBTCcheck$wgratio > 1 |TBTCcheck$nbratio > 1),],file=paste0("TBTCtoCheck_",state,gsa,gen,spec,".csv",sep=""),row.names = F)

# TOTAL ###
TATCn1=TATCn[,c("year","area", "haul_number","strata","sqkm","length_class","nblonraise")]
st=TATCn1$nblonraise/TATCn1$sqkm  #### frequency/swkm x cala
TATC1=cbind(TATCn1,st)

lf=data.frame(aggregate(TATC1$nblonraise,list(TATC1$length_class,TATC1$year,TATC1$strata),sum))
names(lf)=c("LC","year","stratum","Value")
#lf$id3=paste0(lf$Year,lf$Stratum,sep="")
#SW

#### Weighting by strata
LFD_fun=function(lfdata, type){
  LFD=lfdata%>%dplyr::mutate(lfd=W*st)%>%
    dplyr::group_by(year, LC)%>%
    dplyr::summarise(Frequency=sum(lfd))%>%
    dplyr::rename("Length"="LC")
  
  LFD=LFD[complete.cases(LFD),]
  
  lclasses <- as.data.frame(rep(seq(min(LFD$Length),max(LFD$Length),lfstep),times=max(LFD$year)-min(LFD$year)+1))
  names(lclasses) <- "Length"
  lclasses$year <- rep(min(LFD$year):max(LFD$year),each=length(unique(lclasses$Length)))
  
  tempLFD <- merge(lclasses,LFD,by=c("Length","year"),all=T)
  
  tempLFD[is.na(tempLFD)] <- 0
  
  tempLFD <- tempLFD[order(tempLFD$year,tempLFD$Length),]
  
  # LFDJ=subset(LFD,LFD$Length>160) if need to plot subset
  ggplot(tempLFD, aes(y=Frequency, x=Length,col=year))+ geom_bar(stat= "identity")+
    facet_wrap(~year)+
    ggtitle(paste(gen,spec,"LFDs_10-800m_GSA",gsa,state, type))+xlab("Length")+ylab("n/km2")+
    theme(legend.position="none")
  
  if( type=="Total" & is.na(Lmat)==FALSE){
    
    Lref=data.frame(Ind=c("Lmat", "Lopt"), val=c(Lmat, Lopt))
    
    ggplot(tempLFD, aes(y=Frequency, x=Length))+ geom_bar(stat= "identity", fill="deepskyblue3")+
      geom_vline(data=Lref, aes(xintercept = val, linetype=Ind), color="black", size=0.8)+
      #geom_vline(xintercept = Lopt,  color="black", size=1,linetype=1)+
      facet_wrap(~year)+
      ggtitle(paste(gen,spec,"LFDs_10-800m_GSA",gsa,state, type))+xlab("Length")+ylab("n/km2")+
      theme(legend.position="bottom")
    
  }
    
    
  
  ggsave(filename=paste0(state,"_","GSA_",gsa,"LFD_10-800m", type,".png"),width = 10, height = 8, dpi = 150, units = "in")
  
  write.csv(tempLFD,file=paste0(alpha_code,"_GSA",gsa,"_",type, "_LFD.csv",sep=""),row.names = F)
  
  return(tempLFD)
  
}

## Total
lfst=merge(lf,SW,by=c("year","stratum"),all=T)%>%
  dplyr::mutate(st=Value/sweptarea)

LFD=LFD_fun(lfst, "Total")

## Lmat
if(is.na(Lmat)==F){
  
  lfst_threshold=lfst%>%dplyr::filter(LC >= Lmat)
  
  LFD_threshold=LFD_fun(lfst_threshold, "Lmat")
  
}

## Selection
#if(is.na(selection)==F){
#  
#  lfst_sel=lf%>%dplyr::filter(stratum %in%selection)%>%
#    dplyr::left_join(SW_sel, by=c("year", "stratum"))%>%
#    dplyr::mutate(st=Value/sweptarea)
#  
#  LFD_sel=LFD_fun(lfst_sel, "Selection")
#  
#}






##### L95 ######
L95_fun=function(lfdata, swdata, type, side_vrs){
  p<-c(0.95)
  p_names <- map_chr(p, ~paste0( "perc", .x*100))
  p_funs <- purrr::map(p, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% set_names(nm = p_names)
  
  # Tot
  dati <- lfdata %>%dplyr::mutate(Frequenza = round(Frequency*100))
  dati <- expandRows(dati[complete.cases(dati),] , "Frequenza")
  L95 <-dati%>%dplyr::group_by(year)%>%dplyr::summarize_at(vars(Length), funs(!!!p_funs))%>%
    dplyr::rename("Val"="perc95")
  
  Pmega=dati%>%
    dplyr::mutate(mega=ifelse(Length > Lopt, "Y","N"))%>%
    group_by(year, mega)%>%
    dplyr::summarise(int= sum(Frequency))%>%
    spread(key=mega, value=int)%>%
    replace(is.na(.), 0)%>%
    dplyr::mutate(Pmega=Y/(Y+N))
  
  
  ### calculate variability
  #haul
  dati <- TATC1 %>%
    dplyr::mutate(Frequenza = round(st))
  dati <- expandRows(dati[complete.cases(dati),] , "Frequenza")
  L95_haul <-dati%>%dplyr::group_by(year, haul_number,area, strata)%>%dplyr::summarize_at(vars(length_class), funs(!!!p_funs))%>%
    dplyr::rename("Indicator"="perc95", "stratum"="strata")
  
  #strata 
  L95strata<-TATC1 %>%dplyr::group_by(length_class, year,strata )%>%
    dplyr::summarize(value=sum(nblonraise))%>%
    dplyr::rename("stratum"="strata")%>%
    dplyr::right_join(swdata, by=c("year", "stratum"))%>%
    dplyr::mutate(st_stratum= value/sweptarea)
  
  dati <- L95strata %>%dplyr::mutate(Frequenza = round(st_stratum*100))
  dati <- expandRows(dati[complete.cases(dati),] , "Frequenza")
  L95_strat <-dati%>%dplyr::group_by(year, stratum)%>%dplyr::summarize_at(vars(length_class), funs(!!!p_funs))%>%
    dplyr::rename("Val_strat"="perc95")
  
  L95_unc=left_join(L95_haul, L95_strat, by=c("year", "stratum"))%>%
    left_join(., hauls_table2, by=c("year", "stratum"))%>%
    dplyr::mutate(var=(Indicator-Val_strat)^2)%>%
    dplyr::group_by(year,stratum)%>%
    dplyr::mutate(var=(sum(var)/(hauls-1)))%>%
    dplyr::mutate(devst=(sqrt(var)))%>%
    dplyr::distinct(year, stratum,.keep_all=TRUE)%>%
    dplyr::select(year, stratum, Val_strat, devst, var)
  
  # uncertainty based on method for other variables
  uncert_indicator=L95_unc%>%
    right_join(swdata, by=c("year", "stratum"))%>%
    dplyr::mutate(Vt=var*W^2)%>%
    dplyr::filter(!is.na(Vt))%>%
    dplyr::group_by(year)%>%
    dplyr::summarise(var_tot=sum(Vt))%>%
    dplyr::mutate(stdev_tot=sqrt(var_tot))
  
  
  
  
  
  L95_plot=left_join(L95, uncert_indicator, by=c("year"))%>%
    left_join(side_vrs, by=c("year"))%>%
    left_join(Pmega%>%dplyr::select(year,Pmega), by="year")
  
  ggplot(L95_plot, aes(x=year, y=Val)) + 
    geom_line(linetype = "solid",size=1,col="deepskyblue3") +
    geom_point()+
    geom_errorbar(aes(ymin=Val-stdev_tot, ymax=Val+stdev_tot), width=.2,col="red",
                  position=position_dodge(0.05),linetype = "dashed",size=0.75)+
    ylab("mm")+
    ggtitle(paste0(gen,spec,"_","GSA",gsa,"_",state,"_L95", type))
  
 return(L95_plot)
  
 
  
}

pl1=L95_fun(LFD, SW, "Total", side_vars)


write.csv(pl1%>%dplyr::rename("L95"="Val")%>%
            dplyr::mutate(month=lubridate::month(as.Date(mean_doy, origin=paste0(as.numeric(year), "-01-01")))),file=paste0(alpha_code,"_GSA",gsa,"_","Indicators_total.csv",sep=""),row.names = F)

#if(is.na(selection)==F){
#  
#  pl2=L95_fun(LFD_sel, SW_sel, "Selection", side_vars_sel) %>%
#    dplyr::mutate(month=lubridate::month(as.Date(mean_doy, origin=paste0(as.numeric(year), "-01-01"))))
#  
#  write.csv(pl2%>%dplyr::rename("Indicator"="Val"),file=paste0(alpha_code,"_GSA",gsa,"_","L95_selection.csv",sep=""),row.names = F)
#  
#  pl3=pl1%>%
#    dplyr::select(year, Val,stdev_tot)%>%
#    dplyr::mutate(type="tot")%>%
#    bind_rows(pl2%>%
#                dplyr::select(year,Val,stdev_tot)%>%
#                dplyr::mutate(type="selection"))
#  
#  ggplot(pl3, aes(x=year, y=Val, color=type)) + 
#    geom_line(linetype = "solid",size=1.25) +
#    geom_point()+
#    geom_hline(data=pl3%>%dplyr::filter(type=="tot"),aes(yintercept = mean(Val), color=type))+
#    geom_hline(data=pl3%>%dplyr::filter(type=="selection"),aes(yintercept = mean(Val), color=type))+
#    geom_errorbar(aes(ymin=Val-stdev_tot, ymax=Val+stdev_tot), width=.2,
#                  position=position_dodge(0.05),linetype = "dashed",size=0.75)+
#    ylab("mm")+
#    ggtitle(paste0(gen,spec,"_","GSA",gsa,"_",state," L95 total and selected strata:" , selection_txt))
#  
#  
#
#}

if(is.na(Lmat)==F){
  
  pl2=L95_fun(LFD_threshold, SW, "Lmat", side_vars) %>%
    dplyr::mutate(month=lubridate::month(as.Date(mean_doy, origin=paste0(as.numeric(year), "-01-01"))))
  
  write.csv(pl2%>%dplyr::rename("L95"="Val"),file=paste0(alpha_code,"_GSA",gsa,"_","Indicators_sizethreshold.csv",sep=""),row.names = F)
  
  pl3=pl1%>%
    dplyr::select(year, Val,Pmega,stdev_tot)%>%
    dplyr::mutate(type="tot")%>%
    bind_rows(pl2%>%
                dplyr::select(year,Val,Pmega, stdev_tot)%>%
                dplyr::mutate(type="Lmat"))
  
  pl95=ggplot(pl3, aes(x=year, y=Val, color=factor(type))) +
    geom_line(linetype = "solid", size=1.25) +
    scale_colour_manual(values=c("deepskyblue3", "tomato"))+
    geom_point(data=pl3, aes(x=year, y=Val, color=factor(type)))+
    geom_hline(data=pl3%>%dplyr::filter(type=="tot"),aes(yintercept = mean(Val)), color="tomato")+
    geom_hline(data=pl3%>%dplyr::filter(type=="Lmat"),aes(yintercept = mean(Val)), color="deepskyblue3")+
    geom_errorbar(data=pl3, aes(ymin=Val-stdev_tot, ymax=Val+stdev_tot), width=.2,
                  position=position_dodge(0.05),linetype = "dashed",size=0.75)+
    ylab("mm")+
    ggtitle(paste0(gen,spec,"_","GSA",gsa,"_",state," L95 total and on Lmat:" , Lmat, " mm"))+
    guides(color = guide_legend(title = "Data"))
  
  ppmega=ggplot(pl3, aes(x=year, y=Pmega, color=factor(type))) +
    geom_line(linetype = "solid", size=1.25) +
    scale_colour_manual(values=c("deepskyblue3", "tomato"))+
    geom_point(data=pl3, aes(x=year, y=Pmega, color=factor(type)))+
    geom_hline(data=pl3%>%dplyr::filter(type=="tot"),aes(yintercept = mean(Pmega)), color="tomato")+
    geom_hline(data=pl3%>%dplyr::filter(type=="Lmat"),aes(yintercept = mean(Pmega)), color="deepskyblue3")+
    geom_hline(yintercept=0.3, color="black")+
    annotate(geom="text", x=1997, y=0.35, label="Ref.point = 0.3", lynetype=2)+
    ylab("% individuals > Lopt")+
    ggtitle(paste0(gen,spec,"_","GSA",gsa,"_",state," Pmega total and on Lmat:" , Lmat, " mm"))+
    guides(color = guide_legend(title = "Data"))
  
  
  
}

ggsave(plot =pl95, filename=paste0(state,"_","GSA_",gsa,"L95uncertainty.png"),width = 10, height = 8, dpi = 150, units = "in")

ggsave(plot =ppmega, filename=paste0(state,"_","GSA_",gsa,"Pmega.png"),width = 10, height = 8, dpi = 150, units = "in")


toc()
# ######### END OF SCRIPT ###############

