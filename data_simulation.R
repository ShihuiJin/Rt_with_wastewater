#simulation data
nday=180;tw=7
R=1
for(i in 2:nday)
{
  R=c(R,exp(log(R[i-1])+rnorm(1,0,0.02)))
}
plot.ts(R)
library(mgcv)
t=1:nday
R1=as.numeric(predict(gam(R~s(t,k=30))))
plot.ts(R1)

R=exp(sin(seq(0,15,15/179))/6*exp(rnorm(180,0,0.5)))


I=round(rnorm(5,2000,50))
for(i in 1:nday)
{
  n1=min(i+4, n.w)
  I=c(I,R1[i]*round(sum(I[(i+5-n1):(i+4)]*rev(w[1:n1])) ))
}

I1=I[-c(1:5)]
c1=c() #reported case counts
r1=c(0.6,1.25,1.15,1.1,1.1,1,0.8)
for(i in 1:nday)
{
  n1=min(i+4, n.g)
  m.c=psi[ceiling(i/tw)]*r1[(i-1)%%tw+1]*sum(I[(i+5-n1):(i+4)]*rev(g[1:n1]))
  c1=c(c1,rpois(1,m.c))
}

#relative reporting rate
#v1: slight change
psi=exp(sin(0.5+seq(0,15,15/179))/15*exp(rnorm(180,0,0.5)))
#v2: median change
psi=exp(sin(0.5+seq(0,5,5/179))/1.5*exp(rnorm(180,0,0.5)))
#v3: significant change
psi=exp(sin(0.5+pi+seq(0, 30, length.out=nday))/2*exp(rnorm(180,0,0.5)))

psi0=as.numeric(predict(gam(psi~s(t,k=30))))
psi=psi0[seq(4,nday,7)];psi=psi/mean(psi)


#shedding median
n2=50 #sites sampled each day
viral=matrix(1e-6,nday,n2)
for(i in 1:nday)
{
  n1=min(i+4, n.v)
  m.v=sum(I[(i+5-n1):(i+4)]*rev(v[1:n1]))
  viral[i,]=exp(rnorm(n2,log(m.v),1))
}

#data list for the model
N1=ceiling(nday/tw) 
data_list_w <- list(
  N0 = 3, N = nday, EpidemicStart=10,
  tw= tw,
  N1 = N1, biweeks=as.vector(do.call('cbind',lapply(1:N1, function(i) rep(i,tw))))[1:nday],
  weekday=rep(1:7,N1)[1:nday],
  nsd=nday-1, sampleday=2:nday,  #index of sampled days
  f = f, n_f = n.f,  # infection to death distribution
  g = g, n_g = n.g,  # infection to report distribution
  h = h, n_h = n.h,  # infection to hospitalization distribution
  SI = w, n_SI = n.w, # serial interval
  v = v, n_v = n.v, # viral load distribution
  viral = log(viral),  # average community viral load (log scale)
  viral_obs=rep(n2,nday), 
  cases = c1,  # total local cases
  M=n2
)



