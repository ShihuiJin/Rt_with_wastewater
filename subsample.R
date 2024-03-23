#wastewater sub-sampling

#alternative sampling days
{
  stan_code_w <- '
data {
  int <lower=1> N0; // number of days for which to impute infections
  int<lower=1> N; // days of observed data 
  int n_SI;
  real SI[n_SI];
  int n_v;
  real v[n_v]; // distribution of viral load against time
  int<lower=1> M; // number of total communities included in the samples
  real viral[N, M]; // log viral load in community wastewater
  int <lower=0> viral_obs[N]; //idx for test result (numeric)
  int EpidemicStart; // start time for inference
  int nsd; //number of sampled days
  int sampleday[nsd]; //sampled days
 }

parameters {
  /*vector<lower=0>[N1] psi; // scaling parameter for cases
  real<lower=0> psi_sigma;*/ // scaling parameter for cases
  vector<lower=0>[N] nu; // sd for LN distribution of viral load
  real<lower=0> nu0; //mean of sd for LN distribution of viral load
  real<lower=0> nu_sigma; // sd of sd for LN distribution of viral load
  real<lower=0> kappa; //sd for inital case counts
  real<lower=0> y; //initial case count
  real<lower=0> tau; //expected initial case count
  vector<lower=0>[N] Rt;
  real<lower=0> R0; //intial Rt
  real<lower=0> sigma; //Rt variation
}

transformed parameters {
    real convolution;
    vector[N] prediction = rep_vector(0,N);
    vector<lower=0>[N] E_viral = rep_vector(0,N);
    // expected number of infections on day t
    prediction[1:N0] = rep_vector(y,N0); // learn the number of cases in the first N0 days
    for (i in (N0+1):N) {
        convolution=0;
        for(j in 1:min(n_SI, i-1)) {
          convolution += prediction[i-j]*SI[j]; 
        }
        prediction[i] = Rt[i] * convolution;
    }
    
    E_viral[1]= 1e-9;
    for(i in 2:N){
      for(j in 1:min(n_v, i-1)){
          E_viral[i] += prediction[i-j]*v[j];
      }
    }
 }
model {
  tau ~ exponential(0.03);
    y ~ exponential(1.0/tau);
  kappa ~ normal(0,1);
  R0 ~ normal(2, 1);
  Rt[1] ~ normal(R0, kappa); 
  sigma~gamma(0.01,0.01);
  nu0~gamma(0.01,0.01);
  nu_sigma~gamma(0.01,0.01);
  nu~normal(nu0,nu_sigma);
  /*psi[1]~gamma(0.01,0.01);
  psi_sigma~gamma(0.01,0.01);*/
  for(t in 2:N){
    Rt[t] ~ lognormal(log(Rt[t-1]),sigma);
  }
  for(t in sampleday){
    viral[t,1:viral_obs[t]] ~ normal(log(E_viral[t]), nu[t]);
  }
  /*for(t in 2:N1){
      psi[t] ~ lognormal(log(psi[t-1]),psi_sigma);
  }*/
  
}
'
}

  tw=7; N1=ceiling((t2-t1+1)/tw)
  sampleday=which(weekday[t1:t2]%in%c(1:5)) #second subsampling strategy
  sampleday=which(weekday[t1:t2]%in%c(1,4)) #third subsampling strategy
  sampleday=which(weekday[t1:t2]%in%c(1)) #fourth subsampling strategy
  data_list_w <- list(
    N0 = 3, N = t2-t1+1, EpidemicStart=10,
    # tw= tw,
    # N1 = N1, biweeks=as.vector(do.call('cbind',lapply(1:N1, function(i) rep(i,tw))))[1:(t2-t1+1)],
    # weekday=weekday[t1:t2],
    nsd=length(sampleday), sampleday=sampleday,  #index of sampled days
    f = f, n_f = n.f,  # infection to death distribution
    g = g, n_g = n.g,  # infection to report distribution
    h = h, n_h = n.h,  # infection to hospitalization distribution
    SI = w, n_SI = n.w, # serial interval
    v = v, n_v = n.v, # viral load distribution
    viral = log(ww.c.location.1[t1:t2,1:M]+1e-10),  # average community viral load (log scale)
    viral_obs=ww.c.location.na[t1:t2], 
    cases = local.age[t1:t2],  # total local cases
    import = import.age[t1:t2],  # total imported cases
    hosp = hosp.age[t1:t2-t.h+1],  # total imported cases
    deaths = death.age[t1:t2],  # total deaths
    M=M
  )
  
  
  model <- stan_model(model_code = stan_code_w)
  fit <- sampling(model, data = data_list_w, chains = 1, iter = 1e4)

#sampling subsets
  #subsample data
  s.max=50; load(paste0('output.1018.subsample/subsample_list.',s.max))
  t1=741; sampleday=which(ww.c.location.ns[[1]][t1:t2]>0)[-1]
  # tw=7; N1=ceiling((t2-t1+1)/tw)
  data_list_w=list(
    N0 = 3, N = t2-t1+1, EpidemicStart=10,
    # tw= tw,
    # N1 = N1, biweeks=as.vector(do.call('cbind',lapply(1:N1, function(i) rep(i,tw))))[1:(t2-t1+1)],
    # weekday=weekday[t1:t2],
    nsd=length(sampleday), sampleday=sampleday,  #index of sampled days
    f = f, n_f = n.f,  # infection to death distribution
    g = g, n_g = n.g,  # infection to report distribution
    h = h, n_h = n.h,  # infection to hospitalization distribution
    SI = w, n_SI = n.w, # serial interval
    v = v, n_v = n.v, # viral load distribution
    viral = log(ww.c.location.s[[1]][t1:t2,]+1e-10),  # average community viral load (log scale)
    viral_obs=ww.c.location.ns[[1]][t1:t2], 
    cases = local.age[t1:t2],  # total local cases
    import = import.age[t1:t2],  # total imported cases
    hosp = hosp.age[t1:t2-t.h+1],  # total imported cases
    deaths = death.age[t1:t2],  # total deaths
    M=s.max
  )
  
  for(i in problem){
    if(i%%5==0) print(paste0(i,': ',Sys.time()))
    data_list_w$viral=log(ww.c.location.s[[i]][t1:t2,]+1e-10)
    temp1=ww.c.location.ns[[i]][t1:t2]; sampleday=which(temp1>0)[-1]
    data_list_w$viral_obs=temp1
    data_list_w$sampleday=sampleday; data_list_w$nsd=length(sampleday)
    fit <- sampling(model, data = data_list_w, chains = 1, iter = 1e4)
    comp=paste0('9ws-',s.max,'-',i)
    save(fit, file=paste0('output.1018.subsample/rstan.',paste0(comp,collapse = '')))
  }
  
  