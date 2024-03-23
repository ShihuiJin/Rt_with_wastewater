#data importation
#covid19 cases (local)
local=as.data.frame(read.csv('local.csv'))
#covid19 cases (imported)
import=as.data.frame(read.csv('imported.csv'))
#covid19 cases (hospitalized, new, starting from May, 2020)
hosp.new=as.data.frame(read.csv('new.hosp.csv'))
t.h=which(local[,1]==hosp.new[1,1])
#covid19 cases (hospitalized, active, starting from May, 2020)
hosp.act=as.data.frame(read.csv('active.hosp.csv'))
#covid19 cases (deaths)
death=as.data.frame(read.csv('death.csv'))

library(lubridate)
library(dplyr)
library(rstan)
#turn into day of a week
weekday=wday(as.Date(local[,1], tryFormats = '%d/%m/%y'),week_start = getOption("lubridate.week.start", 1))

t10=t1=which(local[,1]=='1/1/22'); t20=t2=which(local[,1]=='31/12/22')
N=t20-t10+1 #total number of days included in the dataset

#distributions
n.age=1
source('dist.r')
#merge age groups 
{
  age.group=c(2:11)
  local.age=rowSums(local[,age.group])
  import.age=rowSums(import[,age.group])
  hosp.age=rowSums(hosp.new[,age.group])
  death.age=rowSums(death[,age.group])
}


#wastewater data
ww.c1=as.data.frame(read.csv('Community_forSharing.csv')[,2:9])
time0=strptime(paste0('23/1/2020',':','0:00'), format = '%d/%m/%Y:%H:%M')
ww.c=ww.c1%>%
  mutate(s.time=strptime(paste0(Date,':',SAMPLING.START.TIME), format = '%d/%m/%Y:%H:%M'), #start time
         e.time=strptime(paste0(Date,':',SAMPLING.END.TIME), format = '%d/%m/%Y:%H:%M'))%>% #end time
  mutate(d.time=as.numeric(difftime(e.time, s.time, units = 'mins')), #duration of sampling time
         rs.time=as.numeric(difftime(s.time, time0, units=c('days')))) #relative start time
location=unique(ww.c$SITE)
#correct the record?
idx.c=which(ww.c$d.time<(-500)); ww.c$d.time[idx.c]=ww.c$d.time[idx.c]+1440
idx.c=which(ww.c$d.time<0); ww.c$d.time[idx.c]=-ww.c$d.time[idx.c]
idx.c=which(ww.c$d.time==0); ww.c$d.time[idx.c]=1
idx.c=which(is.na(ww.c$d.time)); ww.c$d.time[idx.c]=1
idx.c=which(is.na(ww.c$rs.time)); ww.c$d.time[idx.c]=1;ww.c$s.time[idx.c]=strptime(paste0(ww.c$Date[idx.c],':12:00'), format = '%d/%m/%Y:%H:%M')
ww.c$rs.time[idx.c]=with(ww.c,as.numeric(difftime(s.time[idx.c], time0, units=c('days'))))
ww.c$Covid.result=as.factor(ww.c$Covid.result)
#missing data imputation (gam)
ww.c$rm.time=ww.c$rs.time+ww.c$d.time/1440
ww.c$Covid.copy.L=as.numeric(ww.c$Covid.copy.L)
ww.c$Covid.CT=as.numeric(ww.c$Covid.CT)

#impute missing data
fit=lm(Covid.copy.L~exp(-Covid.CT)-1, data=na.omit(ww.c))
ww.c$pred.L=rep(0.001, nrow(ww.c))
idx=with(ww.c, which(is.na(Covid.copy.L)==0))
ww.c$pred.L[idx]=ww.c$Covid.copy.L[idx]
idx=with(ww.c, which(is.na(Covid.CT)==0&is.na(Covid.copy.L)))
ww.c$pred.L[idx]=predict(fit, newdata=ww.c[idx, ])

#aggregate by day
#viral load (L) for a single site
L_day_site=function(site){
  load=matrix(NA,1227,1440) #from Jan 23, 2020 to June 2, 2023, a matrix for each minute
  temp=ww.c%>%filter(SITE==site)%>%filter(is.na(Covid.copy.L)==0)
  if(nrow(temp)==0) return(rep(NA, nrow(load)))
  for(t in 1:nrow(temp)){
    r=ceiling(temp$rs.time[t])
    v=temp$pred.L[t]
    #v=temp$Covid.copy.L[t]
    idx1=1+round((temp$rs.time[t]-r+1)*1440)
    if(idx1>1440) idx1=idx1-1440
    idx2=idx1-1+temp$d.time[t]
    if(idx2<=1440){
      load[r,idx1:idx2]=v
    }else{
      load[r,idx1:1440]=load[r+1,1:(idx2-1440)]=v
    }
  }
  as.vector(do.call('rbind',lapply(1:nrow(load), function(i){
    if(mean(is.na(load[i,]))==1) return(NA)
    mean(load[i,],na.rm=T)
  })))
}
#status: 1 for negative, 2 for equivocal, 3 for positive
L_day_site_sta=function(site){
  load=matrix(NA,1227,1440) #from Jan 23, 2020 to June 2, 2023, a matrix for each minute
  temp=ww.c%>%filter(SITE==site)%>%
    mutate(sta.numeric=as.numeric(as.factor(Covid.result)))
  for(t in 1:nrow(temp)){
    r=ceiling(temp$rs.time[t])
    v=temp$sta.numeric[t];
    if(v<3) v=3-v
    idx1=1+round((temp$rs.time[t]-r+1)*1440)
    if(idx1>1440) idx1=idx1-1440
    idx2=idx1-1+temp$d.time[t]
    if(idx2<=1440){
      load[r,idx1:idx2]=v
    }else{
      load[r,idx1:1440]=load[r+1,1:(idx2-1440)]=v
    }
  }
  as.vector(do.call('rbind',lapply(1:nrow(load), function(i){
    if(mean(is.na(load[i,]))==1) return(NA)
    round(mean(load[i,],na.rm=T))
  })))
}
ww.c.location=as.matrix(do.call('cbind',lapply(location, function(c) L_day_site(c))))
#ww.c.location.sta=as.matrix(do.call('cbind',lapply(location, function(c) L_day_site_sta(c))))
ww.ave.c=rowMeans(ww.c.location,na.rm=T) 
ww.c.location.na=rowSums(is.na(ww.c.location)==0)
ww.c.location.1=as.matrix(do.call('rbind',lapply(1:nrow(ww.c.location), function(i){
  c(ww.c.location[i,is.na(ww.c.location[i,])==0],
    rep(1e-10,sum(is.na(ww.c.location[i,]))))
})))
#M=ncol(ww.c.location)
M=max(ww.c.location.na)

#subsampling
{
  # s.max=50
  # ww.c.location.ns=ww.c.location.na #numbers of subsamples
  # ww.c.location.ns[which(ww.c.location.ns>s.max)]=s.max
  # ww.c.location.s=lapply(1:100, function(n){
  #   set.seed(i)
  #   as.matrix(do.call('rbind',lapply(1:nrow(ww.c.location), function(i){
  #     if(ww.c.location.ns[i]==0) return(rep(1e-10,M))
  #     idx.s=sample(which(is.na(ww.c.location[i,])==0), ww.c.location.ns[i], replace=F)
  #     c(ww.c.location[i,idx.s],rep(1e-10,M-ww.c.location.ns[i]))
  #   })))
  # })
  # save(ww.c.location.ns, ww.c.location.s, file=paste0('output.0821/subsample/subsample.',s.max))
}



