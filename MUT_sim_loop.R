# load the function modified
#library(devtools)
#install_github("duplisea/dublogistic")
#install_github("marchtaylor/fishdynr")
require(fishdynr)
require(tidyverse)
library(ggpubr)
require(magrittr)


remove(list=ls())
source("VirtualPop_new3.R")

# FMSY ----
# test the simulation at FMSY = 0.45
harvest=seq(0,2,0.05)

input_dir="~/CNR/MSFD/github/release"

sspp= "MUT" ## set species 3 alpha code

gsa="17" # Assign gsa code to output file name (i.e. 9_11_ for two GSAs) 

sel="sizethreshold" ## set kind of input file. Alternatives: "selection" ; "total; sizethreshold"

indicator="L95" # alternatives= L95 ; Pmega

lmat=12.2
### set WD #####
dir_t=paste0(input_dir,"/","GSA",gsa,"/", sspp,"/", sep="")
dir.create(paste0(dir_t, "/pop_sym"))
dir.create(paste0(dir_t, "/pop_sim"))

for(i in 1:length(harvest)){
  
  cat(paste(i," of " ,length(harvest)))
  
  MUT_FMSY <- virtualPop_new3(
    tincr = 1/12,
    K.mu = 0.247, K.cv = 0.05,
    Linf.mu = 29.185, Linf.cv = 0.05,
    ts = 0, C = 0,
    LWa = 0.000012, LWb = 3.015,
    Lmat = 10.6, wmat = 1,
    rmax = 500000, beta = 1,
    harvest_rate = harvest[i], 
    L50 = 13.73,
    l1 = 8,
    bin.size = 1,
    wqs = 2.19,
    timemin = 0, timemax = 35, timemin.date = as.Date("1985-01-01"),
    N0 = 500000,
    fished_t = seq(15, 35, 1/12),
    lfqFrac = 1,
    Lrec = 4,
    Lopt = 20,
    age = c(0, 1, 2, 3, 4), mortality = c(1.41, 0.71, 0.52, 0.42, 0.37),
    repro_wt = c(0.0628707, 0.159090909, 0.194775435, 0.579536, 0.898477, 0.798611, 
                 0.654121864, 0.137440758, 0.015881, 0.05703, 0.020263213, 0.039121),
    selectivity = "logistic",
    growth_type = "Fast_growth_K",
    progressBar = TRUE)
  
  saveRDS(MUT_FMSY, file = paste0(dir_t, "pop_sym/","MUT_sim_F", harvest[i], ".rds"))
  
}


temp=list.files(path=paste0(dir_t, "pop_sym/"))
setwd(paste0(dir_t, "pop_sym/"))
myf=lapply(temp, readRDS )



res=tibble(time=as.Date(character()), month=as.numeric(), F_ref=as.character(), trend=as.numeric() , seasonal=as.numeric(), effect=as.numeric(), error=as.numeric())

biom=tibble(length=as.numeric(), n=as.numeric(), biomass=as.numeric(), F_ref=as.character())

sea_pop=tibble(biomass=as.numeric(), n=as.integer(), length=as.numeric(), month=as.numeric(), F_ref=as.character())

catch=tibble(size=as.integer(), n=as.numeric(), n_rel=as.numeric(), F_ref=as.character())


for(i in 1:length(temp)){
  
  MUT_FMSY=myf[[i]]
  
  ### sea pop
  
  sea_pop_F03 <- MUT_FMSY$sea_pop
  
  for (j in 1:length(sea_pop_F03)){
    
    sea_pop_F03[[j]]$month <- j +1
    
  }
  
  sea_F03_df <- sea_pop_F03 %>% map_df(as_tibble)
  sea_F03_df$month <- ifelse(sea_F03_df$month == 13, 1, sea_F03_df$month)

  sea_F03_df=sea_F03_df%>%
    select(biomass, n, length, month)%>%
    dplyr::mutate(F_ref=str_remove(str_remove(temp[i], "MUT_sim_"), ".rds"))
  
  sea_pop=bind_rows(sea_pop, sea_F03_df)
  
  ##### Biomass
  sea_pop_F00 <- MUT_FMSY$sea_pop
  
  sea_pop_F00_df <- sea_pop_F00 %>% map_df(as_tibble)
  sea_pop_F00_df %<>% 
    mutate(tot_biomass = sum(biomass),
           tot_n = sum(n)) %>% 
    group_by(length) %>% 
    summarise(n = sum(n)/max(tot_n),
              biomass = sum(biomass)/ max(tot_biomass))%>%
    dplyr::mutate(F_ref=str_remove(str_remove(temp[i], "MUT_sim_"), ".rds"))
  
  biom=bind_rows(biom, sea_pop_F00_df)
  
  ##### Catches
  
  catch_F03 <- MUT_FMSY$inds[230:241]
  if(is.null(catch_F03)==TRUE){
    
    print("no catches")
    
  }else{
    
    catch_F03 <- bind_rows(catch_F03)# %>% map_df(as_tibble)
    catch_F03 %<>% group_by(length) %>% 
      summarise(n = sum(n)) %>% 
      ungroup()
    
    catch_F03= catch_F03 %>% dplyr::rename("size" = "length") %>% 
      mutate(n_rel = n/sum(n), F_ref=str_remove(str_remove(temp[i], "MUT_sim_"), ".rds"))
    
    catch =  bind_rows(catch, catch_F03)
    
    
  }
  

  ##### Indicator
  
  
  
  datt=bind_rows(MUT_FMSY[["pop"]])
  
  myts <- ts(data = datt$L95, frequency = 12, start = datt$dates[1])
  
  ## weights for moving avg 
  fltr <- c(1/2, rep(1, times = 11), 1/2)/12
  
  myts.trend <- stats::filter(myts, filter = fltr, method = "convo", sides = 2) ## get the trend 
  
  dat_plot=tibble(time=datt$dates, month=c(rep(seq(1:12), 35),1), F_ref=str_remove(str_remove(temp[i], "MUT_sim_"), ".rds") , trend=as.numeric(myts.trend))
  
 
  ## seasonal effect over time 
  myts.1T <- myts - myts.trend
  
  
  dat_plot$seasonal=as.numeric(myts.1T)
  
  ll <- length(myts.1T) ## frequency (ie, 12) 
  ff <- frequency(myts.1T) ## number of periods (years); %/% is integer division 
  periods <- ll%/%ff ## index of cumulative month 
  index <- seq(1, ll, by = ff) - 1 ## get mean by month 
  mm <- numeric(ff)
  
  for (k in 1:ff) { mm[k] <- mean(myts.1T[index + k], na.rm = TRUE)
  }
  ## subtract mean to make overall mean=0 
  mm <- mm - mean(mm)
  
  dat_plot_2=tibble(month=seq(1:12), effect=as.numeric(mm))
  
  dat_plot=dat_plot%>%left_join(dat_plot_2, by="month")
  
  myts.seas <- ts(rep(mm, periods + 1)[seq(ll)], start = start(myts.1T), frequency = ff)
  
  myts.err <- myts - myts.trend - myts.seas
  
  dat_plot$error=as.numeric(myts.err)
  
  res=bind_rows(res,dat_plot) 
  
  
}


Fstep=seq(0,1.9,0.1)
F01=0.41

ggplot(data=res%>%dplyr::filter(as.numeric(str_remove(F_ref, "F"))%in%Fstep))+
  geom_line(aes(x=time, y=trend, color=F_ref)) +  ylab("mm")+
  ggtitle(paste0(sspp,"_","GSA",gsa,"_L95 trend,simulated population"))+
  scale_y_continuous(breaks=seq(5,20,1))+ labs(color = "F")+
  scale_x_date(breaks = seq(min(res$time), max(res$time),366*2), labels=lubridate::year(seq(min(res$time), max(res$time),366*2)))

ggsave(paste0(dir_t, "pop_sim/L95_trend.png"), width = 30, height = 15, units = "cm")

L95Ref=res%>%
  dplyr::filter(lubridate::year(time) > 2015)%>%
  dplyr::filter(!is.na(trend))%>%
  dplyr::group_by(F_ref)%>%
  mutate(val=mean(trend))%>%
  distinct(F_ref, val)%>%
  dplyr::filter(as.numeric(str_remove(F_ref, "F"))< F01)

write.csv(L95Ref, paste0(dir_t, "pop_sim/L95RP.csv"), row.names = F)


ggplot(data=res%>%dplyr::filter(as.numeric(str_remove(F_ref, "F"))%in%Fstep))+geom_line(aes(x=month, y=effect, color=F_ref))+
  ylab("effect (mm)")+
  ggtitle(paste0(sspp,"_","GSA",gsa,"_L95 seasonal effect,simulated population"))+
  scale_x_continuous(breaks=seq(1,12,1))+ labs(color = "F")

ggsave(paste0(dir_t, "pop_sim/L95_seasonal.png"), width = 30, height = 15, units = "cm")

biom %>% dplyr::filter(as.numeric(str_remove(F_ref, "F"))%in%Fstep)%>%
  ggplot(aes(x = length, y  = biomass, color=F_ref)) +
  geom_point() +
  geom_smooth(method = "gam", se=FALSE) +
  geom_vline(xintercept = 20.84, color = "firebrick1", linetype = "twodash") +
  annotate(geom = "text", x = 20.84, y = 0, label = "L[opt]", output = 'character', parse = TRUE) +
  geom_vline(xintercept = lmat, color = "firebrick1", linetype = "twodash") +
  annotate(geom = "text", x = lmat, y = 0, label = "L[mat]", output = 'character', parse = TRUE) +
  labs(x = "Length", y = "Relative Biomass")+ggtitle(paste0(sspp,"_","GSA",gsa,"cumulated LFD, simulated population"))

ggsave(paste0(dir_t, "pop_sim/cumulated_lfd.png"), width = 30, height = 15, units = "cm")

#### Catches comparison
library(readxl)
MUT_comm <- read_excel("~/CNR/Stock Assessment/GFCM/2019/MUT/Data/Input_ss3/Input_ENAmese1.xlsx", 
                       sheet = "LFD")%>%
  dplyr::filter(fleet==1)%>% 
  dplyr::select(-month, -fleet,-sex,-part,-Nsamp)%>%
  gather(key = "len", value = "n", -"#_yr")%>%
  dplyr::mutate(len=as.numeric(len), n=as.numeric(n))%>%
  dplyr::rename("year"="#_yr" )


MUT_comm <- na.omit(MUT_comm)

MUT_comm2 <- MUT_comm %>% filter(year >= 2016) %>% 
  mutate(data = "CAMPBIOL")

MUT_comm2 %<>%
  group_by(year) %>% 
  mutate(n_rel = n/sum(n)) %>% 
  ungroup()


#### MEDITS comparison
MUT_MEDITS <- read_csv("~/CNR/MSFD/github/release/GSA17/MUT/MUT_GSA17_Total_LFD.csv")
mon <- read_csv("~/CNR/MSFD/github/release/GSA17/MUT/MUT_GSA17_Indicators_total.csv")%>%select(year, month)

MUT_MEDITS2=MUT_MEDITS%>% 
  left_join(mon, by="year")%>%
  #filter(Length > 122) %>% 
  group_by(year, month) %>% 
  mutate(n_tot = sum(Frequency)) %>% 
  ungroup() %>% 
  group_by(Length, year, month) %>% 
  mutate(n_rel = Frequency/max(n_tot))%>% 
  mutate(Survey = "Medits")

sea_pop2=sea_pop%>%
  filter(length >=4)%>%
  dplyr::filter(as.numeric(str_remove(F_ref, "F"))%in%Fstep)%>%
  right_join(mon, by="month")%>%
  group_by(year, F_ref)%>%
  dplyr::mutate(n_rel=n/sum(n), length=length*10)

  

catch2=catch%>%dplyr::filter(as.numeric(str_remove(F_ref, "F"))%in%Fstep)

p_c=ggplot(catch2)+geom_line(aes(x=size, y=n_rel, color=F_ref), size=1)+
  geom_bar(data = MUT_comm2, mapping = aes(x = len , y = n_rel ), stat = "identity",
           color = "grey30", alpha = 0.5)+
  facet_wrap(~year)+ggtitle(paste0(sspp,"_","GSA",gsa," catches LFD, observed (CAMPBIOL: grey bars) vs simulated"))+
  ylab("%")+guides(color = guide_legend(nrow = 2))


p_s=ggplot(sea_pop2)+
  geom_bar(data = MUT_MEDITS2, mapping = aes(x = Length , y = n_rel, fill=Survey), stat = "identity",
           color = "grey30", alpha = 0.5)+
  geom_line(aes(x=length, y=n_rel, color=F_ref), size=1)+
  facet_wrap(~year)+
  ggtitle(paste0(sspp,"_","GSA",gsa," population LFD, observed (MEDITS: pink bars) vs simulated"))+
  ylab("%")+xlab("Length (mm)")+guides(color = guide_legend(nrow = 2))


ggarrange(p_c, p_s, nrow = 2, common.legend = TRUE, legend = "bottom", heights = c(0.5, 1.5))

ggsave(paste0(dir_t, "pop_sim/lfd_comparison.png"), width = 30, height = 30, units = "cm")





################### ------- end 
