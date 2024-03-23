#output analysis
t10=1044;t20=1225 #Dec 2022 to May 2023 
t10=710; t20=1074 #2022

#comparison functions
{
  #comparison (cases/hosp/etc)
  comparison=function(estim, obs){
    #pred has six columns: mean, median, 95%CrI, IQR
    #point estimate is posterior mean
    #obs is observation
    obs.0=which(abs(obs)>=1e-6)
    mae=mean(abs(estim[,1]-obs))
    mape=mean(abs(estim[obs.0,1]-obs[obs.0])/obs[obs.0])
    cover.iq=mean(estim[,5]<=obs&estim[,6]>=obs) #IQR coverage
    cover.95=mean(estim[,3]<=obs&estim[,4]>=obs) #95% CrI coverage
    c(mae, mape*100, cover.iq*100, cover.95*100)
  }
  #viral load comparison
  comparison.v=function(estim.raw, estim.log, t0=10){
    #pred has six columns: mean, median, 95%CrI, IQR
    #point estimate is posterior median
    #obs is observation
    m.viral=as.vector(do.call('cbind',lapply((t1+t0):t2, function(i){
      median(ww.c.location[i,],na.rm=T)
    })))
    lm.viral=as.vector(do.call('cbind',lapply((t1+t0):t2, function(i){
      median(log(ww.c.location[i,]),na.rm=T)
    })))
    obs.0=which(abs(lm.viral)>=1e-6)
    mae=mean(abs(estim.log[,1]-lm.viral))
    mape=mean(abs(estim.log[obs.0,1]-lm.viral[obs.0])/lm.viral[obs.0])
    #IQR coverage   
    cover.iq=do.call('mean',lapply((t1+t0):t2, function(t){
      mean(estim.raw[t-t0-t1+1,5]<=ww.c.location[t,]&estim.raw[t-t0-t1+1,6]>=ww.c.location[t,], na.rm=T)
    }))
    #95% CrI coverage   
    cover.95=do.call('mean',lapply((t1+t0):t2, function(t){
      mean(estim.raw[t-t0-t1+1,3]<=ww.c.location[t,]&estim.raw[t-t0-t1+1,4]>=ww.c.location[t,], na.rm=T)
    }))
    c(mae, mape*100, cover.iq, cover.95)
  }
  
}


n.sim=5e3; t0=10
draws=extract(fit)

#Rt
Rt=as.matrix(do.call('rbind',lapply(1:ncol(draws$Rt), function(i){
  temp=draws$Rt[,i]
  c(mean(temp), quantile(temp,c(0.5,0.025,0.975,0.25,0.75)))
})))

#viral load
Vt.raw=as.matrix(do.call('rbind',lapply(1:(t2-t1+1), function(i){
  temp=rlnorm(n.sim*10, log(draws$E_viral[,i]), draws$nu[,i])
  c(mean(temp), quantile(temp,c(0.5,0.025,0.975,0.25,0.75)))
})))
Vt.l=as.matrix(do.call('rbind',lapply(1:(t2-t1+1), function(i){
  temp=log(draws$E_viral[,i])
  c(mean(temp), quantile(temp,c(0.5,0.025,0.975,0.25,0.75)))
})))
#comparison with observations
comparison.v(Vt.raw[-c(1:t0),], Vt.l[-c(1:t0),])

#case
It1.raw=as.matrix(do.call('rbind',lapply(1:(t2-t1+1), function(i){
  temp1=draws$E_cases[,i];temp2=draws$phi[,1]
  temp=rnbinom(length(temp1), size=temp2, mu=temp1)
  c(mean(temp), quantile(temp,c(0.5,0.025,0.975,0.25,0.75)))
})))
comparison(It1.raw[-c(1:t0),], local$Total[(t1+t0):t2])

#hospitalization
Ht.raw=as.matrix(do.call('rbind',lapply(1:(t2-t1+1), function(i){
  temp1=draws$E_hosp[,i];temp2=draws$phi[,4]
  temp=rnbinom(length(temp1), size=temp2, mu=temp1)
  c(mean(temp), quantile(temp,c(0.5,0.025,0.975,0.25,0.75)))
})))
comparison(Ht.raw[-c(1:t0),], hosp.new$Total[(t1+t0):t2-t.h+1])

#day-of-the-week effects for reported case counts
r=as.matrix(do.call('rbind',lapply(1:7, function(i){
  temp=draws$r[,i]
  c(mean(temp), quantile(temp,c(0.5,0.025,0.975,0.25,0.75,0,1)))
})))
r1=round(r,2)

#time-varying case-hospitalization rates
{
  tw=7;N1=ceiling((t2-t1+1)/tw)
  temp=as.matrix(do.call('rbind',lapply(1:N1, function(i){
    temp1=draws$rh[,i]
    c(mean(temp1), quantile(temp1,c(0.5,0.025,0.975,0.25,0.75)))
  })))
  #average case-hospitalization rate throughout the inference window
  temp1=(rowSums(draws$rh[1:(N1-1),])*7+((t2-t1)%%tw+1)*draws$rh[N1,])/(t2-t1+1)
  c(mean(temp1), quantile(temp1,c(0.5,0.025,0.975,0.25,0.75)))
  #fatality
  temp=as.matrix(do.call('rbind',lapply(1:N1, function(i){
    temp1=draws$rd[,i]
    c(mean(temp1), quantile(temp1,c(0.5,0.025,0.975,0.25,0.75)))
  })))
  #average case-fatality rate throughout the inference window
  temp1=(rowSums(draws$rd[1:(N1-1),])*7+((t2-t1)%%tw+1)*draws$rd[N1,])/(t2-t1+1)
  c(mean(temp1), quantile(temp1,c(0.5,0.025,0.975,0.25,0.75)))
}

#psi
Psi=as.matrix(do.call('rbind',lapply(1:ceiling((t2-t1+1)/tw), function(i){
  temp=draws$psi[,i]
  temp1=c(mean(temp), quantile(temp,c(0.5,0.025,0.975,0.25,0.75)))
  tm=tw
  if(i*tw>(t2-t1+1)) tm=t2-t1+1-(i-1)*tw
  as.matrix(do.call('rbind',lapply(1:tm, function(j){
    temp1
  })))
})))

