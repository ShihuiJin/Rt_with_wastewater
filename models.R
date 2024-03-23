#wastewater model
{
  stan_code<- '
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

#case model
{
stan_code <- '
data {
  int <lower=1> N0; // number of days for which to impute infections
  int<lower=1> N; // days of observed data 
  int cases[N]; // reported local case counts
  //int import[N]; // reported imported case counts
  int n_g;
  real g[n_g]; // distribution of time from infection to report 
  int n_SI;
  real SI[n_SI]; // fixed pre-calculated SI 
  /*int n_h;
  real h[n_h]; // distribution of time from infection to hospitalization 
  int hosp[N]; // reported new hospitalizations
  int n_f;
  real f[n_f]; // distribution of time from infection to death, by age group 
  int deaths[N];*/ // reported deaths
  int weekday[N]; // which day of a week 
  int EpidemicStart; // start time for inference
 }

parameters {
  real<lower=0> kappa; //sd for inital case counts
  real<lower=0> y; //initial case count
  vector<lower=0>[4] phi; //negative binomial parameter
  real<lower=0> tau; //expected initial case count
  vector<lower=0>[N] Rt;
  real<lower=0> R0; //intial Rt
  real<lower=0> sigma; //Rt variation
  //real<lower=0>[N] I_import; //latent importation
  simplex[7] r11; //scaled relative reporting rate for each week, local cases
  //simplex[7] r12; //scaled relative reporting rate for each week, imported cases
  /*simplex[7] r21; //scaled relative reporting rate for each week, hospitalised cases
  real<lower=0> rh; // case-hospitalisation rates by age
  real<lower=0> rd;*/ // case-fatality rates by age
}

transformed parameters {
    real convolution;
    vector[N] prediction = rep_vector(0,N);
    vector[N] E_cases  = rep_vector(0,N);
    /*vector[N] E_import  = rep_vector(0,N);
    vector[N] E_hosp  = rep_vector(0,N);
    vector[N] E_deaths  = rep_vector(0,N);*/
    vector[7] r = r11*7;
    /*vector[7] r1 = r12*7;
    vector[7] r2  = r21*7;*/
    // expected case counts on day t
    prediction[1:N0] = rep_vector(y,N0); // learn the number of cases in the first N0 days
    for (i in (N0+1):N) {
        convolution=0;
        for(j in 1:min(n_SI, i-1)) {
          //convolution += (prediction[i-j]+I_import[i-j])*SI[j]; 
          convolution += prediction[i-j]*SI[j]; 
        }
        prediction[i] = Rt[i] * convolution;
    }
    {
      E_cases[1]=1e-9;
      /*E_import[1]=1e-9;
      E_hosp[1]=1e-9;
      E_deaths[1]= 1e-9;*/
      for (i in 2:N){
        E_cases[i]=0;
        //E_import[i]=0;
        for(j in 1:min(n_g, i-1)){
          E_cases[i] += prediction[i-j]*g[j];
          //E_import[i] += I_import[i-j]*g[j];
        }
        E_cases[i]*=r[weekday[i]];
        //E_import[i]*=r1[weekday[i]];
        /*E_hosp[i]=0;
        for(j in 1:min(n_h, i-1)){
          E_hosp[i] += prediction[i-j]*h[j];
        }
        E_hosp[i]*=r2[weekday[i]]*rh; //different hospitalization rates
        E_deaths[i]= 0;
        for(j in 1:min(n_f, i-1)){
          E_deaths[i] += prediction[i-j]*f[j]; //different distribution from onset to death for different age groups
        }
        E_deaths[i]*=rd;*/
      }
    }
 }
model {
  /*rd ~ normal(0,0.1);
  rh ~ normal(0,0.1);*/
  tau ~ exponential(0.03);
  r11 ~ dirichlet(rep_vector(1, 7));
  //r12 ~ dirichlet(rep_vector(1, 7));
  //r21 ~ dirichlet(rep_vector(1, 7));
  y ~ exponential(1.0/tau);
  kappa ~ normal(0,1);
  R0 ~ normal(2, 1);
  Rt[1] ~ normal(R0, kappa); 
  sigma~gamma(0.01,0.01);
 
  for(t in 2:N){
    Rt[t] ~ lognormal(log(Rt[t-1]),sigma);
  }
  {
    phi ~ normal(0,5);
    /*I_import[1:2]~normal(150,75);
    for(t in 3:N){
      I_import[t]~normal(import[t-2], max(import[t-2]/2,1));
    }*/
    for(i in EpidemicStart:N){
      cases[i] ~ neg_binomial_2(E_cases[i],phi[1]);
      //import[i] ~ neg_binomial_2(E_import[i],phi[2]);
      /*deaths[i] ~ neg_binomial_2(E_deaths[i],phi[3]); 
      hosp[i] ~ neg_binomial_2(E_hosp[i],phi[4]);*/
    }
  }
}
'
}

#case and wastewater model
{
stan_code <- '
data {
  int <lower=1> N0; // number of days for which to impute infections
  int<lower=1> N; // days of observed data 
  int cases[N]; // reported local case counts
  //int import[N]; // reported imported case counts
  int n_g;
  real g[n_g]; // distribution of time from infection to report 
  int n_SI;
  real SI[n_SI]; // fixed pre-calculated SI 
  /*int n_h;
  real h[n_h]; // distribution of time from infection to hospitalization 
  int hosp[N]; // reported new hospitalizations
  int n_f;
  real f[n_f]; // distribution of time from infection to death, by age group 
  int deaths[N];*/ // reported deaths
  int n_v;
  real v[n_v]; // distribution of viral load against time
  int<lower=1> M; // number of total communities included in the samples
  real viral[N, M]; // log viral load in community wastewater
  int <lower=0> viral_obs[N]; //idx for test result (numeric)
  int weekday[N]; // which day of a week 
  int EpidemicStart; // start time for inference
  int tw; //number of days in a short window when reporting rate is constant
  int N1; // number of biweeks
  int biweeks[N]; //biweek index
  int nsd; //number of sampled days
  int sampleday[nsd]; //sampled days
 }

parameters {
  vector<lower=0>[N1] psi; // scaling parameter for cases
  real<lower=0> psi_sigma; // scaling parameter for cases
  vector<lower=0>[N] nu; // sd for LN distribution of viral load
  real<lower=0> nu0; //mean of sd for LN distribution of viral load
  real<lower=0> nu_sigma; // sd of sd for LN distribution of viral load
  real<lower=0> kappa; //sd for inital case counts
  real<lower=0> y; //initial case count
  vector<lower=0>[4] phi; //negative binomial parameter
  real<lower=0> tau; //expected initial case count
  vector<lower=0>[N] Rt;
  real<lower=0> R0; //intial Rt
  real<lower=0> sigma; //Rt variation
  //vector<lower=0>[N] I_import; //latent importation
  simplex[7] r11; //scaled relative reporting rate for each week, local cases
  /*simplex[7] r12; //scaled relative reporting rate for each week, imported cases
  simplex[7] r21; //scaled relative reporting rate for each week, hospitalised cases
  real<lower=0> rh; // case-hospitalisation rates by age
  real<lower=0> rd;*/ // case-fatality rates by age
}

transformed parameters {
    real convolution;
    vector[N] prediction = rep_vector(0,N);
    vector[N] E_cases  = rep_vector(0,N);
    /*vector[N] E_import  = rep_vector(0,N);
    vector[N] E_hosp  = rep_vector(0,N);
    vector[N] E_deaths  = rep_vector(0,N);*/
    vector<lower=0>[N] E_viral = rep_vector(0,N);
    vector[7] r = r11*7;
    /*vector[7] r1 = r12*7;
    vector[7] r2  = r21*7;*/
    // expected case counts on day t
    prediction[1:N0] = rep_vector(y,N0); // learn the number of cases in the first N0 days
    for (i in (N0+1):N) {
        convolution=0;
        for(j in 1:min(n_SI, i-1)) {
          //convolution += (prediction[i-j]+I_import[i-j]/psi[biweeks[i-j]])*SI[j]; 
          convolution += prediction[i-j]*SI[j]; 
        }
        prediction[i] = Rt[i] * convolution;
    }
    {
      E_cases[1]=1e-9;
      /*E_import[1]=1e-9;
      E_hosp[1]=1e-9;
      E_deaths[1]= 1e-9;*/
      for (i in 2:N){
        E_cases[i]=0;
        //E_import[i]=0;
        for(j in 1:min(n_g, i-1)){
          E_cases[i] += prediction[i-j]*g[j];
          //E_import[i] += I_import[i-j]*g[j];
        }
        E_cases[i]*=r[weekday[i]]*psi[biweeks[i]];
        /*E_import[i]*=r1[weekday[i]];
        E_hosp[i]=0;
        for(j in 1:min(n_h, i-1)){
          E_hosp[i] += prediction[i-j]*h[j];
        }
        E_hosp[i]*=r2[weekday[i]]*psi[biweeks[i]]; //different hospitalization rates
        E_hosp[i]*=r2[weekday[i]]*rh*psi[biweeks[i]]; //different hospitalization rates
        E_deaths[i]= 0;
        for(j in 1:min(n_f, i-1)){
          E_deaths[i] += prediction[i-j]*f[j]; //different distribution from onset to death for different age groups
        }
        E_deaths[i]*=rd*psi[biweeks[i]];*/
      }
    }
    E_viral[1]= 1e-9;
    for(i in 2:N){
      for(j in 1:min(n_v, i-1)){
          E_viral[i] += prediction[i-j]*v[j];
      }
      //E_viral[i]*=psi;
    }
 }
model {
  /*rd ~ normal(0,0.1);
  rh ~ normal(0,0.1);*/
  tau ~ exponential(0.03);
  r11 ~ dirichlet(rep_vector(1, 7));
  /*r12 ~ dirichlet(rep_vector(1, 7));
  r21 ~ dirichlet(rep_vector(1, 7));*/
  y ~ exponential(1.0/tau);
  kappa ~ normal(0,1);
  R0 ~ normal(2, 1);
  Rt[1] ~ normal(R0, kappa); 
  sigma~gamma(0.01,0.01);
  nu0~gamma(0.01,0.01);
  nu_sigma~gamma(0.01,0.01);
  nu~normal(nu0,nu_sigma);
  psi[1]~gamma(0.01,0.01);
  psi_sigma~gamma(0.01,0.01);
  for(t in 2:N){
    Rt[t] ~ lognormal(log(Rt[t-1]),sigma);
  }
  for(t in sampleday){
    viral[t,1:viral_obs[t]] ~ normal(log(E_viral[t]), nu[t]);
  }
  for(t in 2:N1){
      psi[t] ~ lognormal(log(psi[t-1]),psi_sigma);
  }
  {
    phi ~ normal(0,5);
    /*I_import[1:2]~normal(150,75);
    for(t in 3:N){
      I_import[t]~normal(import[t-2], max(import[t-2]/2,1));
    }*/
    for(i in EpidemicStart:N){
      cases[i] ~ neg_binomial_2(E_cases[i],phi[1]);
      /*import[i] ~ neg_binomial_2(E_import[i],phi[2]);
      deaths[i] ~ neg_binomial_2(E_deaths[i],phi[3]);
      hosp[i] ~ neg_binomial_2(E_hosp[i],phi[4]);*/
    }
  }
}
'
} 
  
#hospital and wastewater model
{
stan_code <- '
data {
  int <lower=1> N0; // number of days for which to impute infections
  int<lower=1> N; // days of observed data 
  /*int cases[N]; // reported local case counts
  int import[N]; // reported imported case counts
  int n_g;
  real g[n_g];*/ // distribution of time from infection to report 
  int n_SI;
  real SI[n_SI]; // fixed pre-calculated SI 
  int n_h;
  real h[n_h]; // distribution of time from infection to hospitalization 
  int hosp[N]; // reported new hospitalizations
  /*int n_f;
  real f[n_f]; // distribution of time from infection to death, by age group 
  int deaths[N];*/ // reported deaths
  int n_v;
  real v[n_v]; // distribution of viral load against time
  int<lower=1> M; // number of total communities included in the samples
  real viral[N, M]; // log viral load in community wastewater
  int <lower=0> viral_obs[N]; //idx for test result (numeric)
  int weekday[N]; // which day of a week 
  int EpidemicStart; // start time for inference
  int tw; //number of days in a short window when reporting rate is constant
  int N1; // number of biweeks
  int biweeks[N]; //biweek index
  int nsd; //number of sampled days
  int sampleday[nsd]; //sampled days
 }

parameters {
  vector<lower=0>[N1] psi; // scaling parameter for cases
  real<lower=0> psi_sigma; // scaling parameter for cases
  vector<lower=0>[N] nu; // sd for LN distribution of viral load
  real<lower=0> nu0; //mean of sd for LN distribution of viral load
  real<lower=0> nu_sigma; // sd of sd for LN distribution of viral load
  real<lower=0> kappa; //sd for inital case counts
  real<lower=0> y; //initial case count
  vector<lower=0>[4] phi; //negative binomial parameter
  real<lower=0> tau; //expected initial case count
  vector<lower=0>[N] Rt;
  real<lower=0> R0; //intial Rt
  real<lower=0> sigma; //Rt variation
  /*vector<lower=0>[N] I_import; //latent importation
  simplex[7] r11; //scaled relative reporting rate for each week, local cases
  simplex[7] r12;*/ //scaled relative reporting rate for each week, imported cases
  simplex[7] r21; //scaled relative reporting rate for each week, hospitalised cases
  /*real<lower=0> rh; // case-hospitalisation rates by age
  real<lower=0> rd;*/ // case-fatality rates by age
}

transformed parameters {
    real convolution;
    vector[N] prediction = rep_vector(0,N);
    /*vector[N] E_cases  = rep_vector(0,N);
    vector[N] E_import  = rep_vector(0,N);*/
    vector[N] E_hosp  = rep_vector(0,N);
    //vector[N] E_deaths  = rep_vector(0,N);
    vector<lower=0>[N] E_viral = rep_vector(0,N);
    /*vector[7] r = r11*7;
    vector[7] r1 = r12*7;*/
    vector[7] r2  = r21*7;
    // expected case counts on day t
    prediction[1:N0] = rep_vector(y,N0); // learn the number of cases in the first N0 days
    for (i in (N0+1):N) {
        convolution=0;
        for(j in 1:min(n_SI, i-1)) {
          //convolution += (prediction[i-j]+I_import[i-j]/psi[biweeks[i-j]])*SI[j]; 
          convolution += prediction[i-j]*SI[j]; 
        }
        prediction[i] = Rt[i] * convolution;
    }
    {
     /*E_cases[1]=1e-9;
      E_import[1]=1e-9;*/
      E_hosp[1]=1e-9;
      //E_deaths[1]= 1e-9;
      for (i in 2:N){
        /*E_cases[i]=0;
        E_import[i]=0;
        for(j in 1:min(n_g, i-1)){
          E_cases[i] += prediction[i-j]*g[j];
          E_import[i] += I_import[i-j]*g[j];
        }
        E_cases[i]*=r[weekday[i]]*psi[biweeks[i]];
        E_import[i]*=r1[weekday[i]];*/
        E_hosp[i]=0;
        for(j in 1:min(n_h, i-1)){
          E_hosp[i] += prediction[i-j]*h[j];
        }
        E_hosp[i]*=r2[weekday[i]]*psi[biweeks[i]]; //different hospitalization rates
        /*E_hosp[i]*=r2[weekday[i]]*rh*psi[biweeks[i]]; //different hospitalization rates
        E_deaths[i]= 0;
        for(j in 1:min(n_f, i-1)){
          E_deaths[i] += prediction[i-j]*f[j]; //different distribution from onset to death for different age groups
        }
        E_deaths[i]*=rd*psi[biweeks[i]];*/
      }
    }
    E_viral[1]= 1e-9;
    for(i in 2:N){
      for(j in 1:min(n_v, i-1)){
          E_viral[i] += prediction[i-j]*v[j];
      }
      //E_viral[i]*=psi;
    }
 }
model {
  /*rd ~ normal(0,0.1);
  rh ~ normal(0,0.1);*/
  tau ~ exponential(0.03);
  /*r11 ~ dirichlet(rep_vector(1, 7));
  r12 ~ dirichlet(rep_vector(1, 7));*/
  r21 ~ dirichlet(rep_vector(1, 7));
  y ~ exponential(1.0/tau);
  kappa ~ normal(0,1);
  R0 ~ normal(2, 1);
  Rt[1] ~ normal(R0, kappa); 
  sigma~gamma(0.01,0.01);
  nu0~gamma(0.01,0.01);
  nu_sigma~gamma(0.01,0.01);
  nu~normal(nu0,nu_sigma);
  psi[1]~gamma(0.01,0.01);
  psi_sigma~gamma(0.01,0.01);
  for(t in 2:N){
    Rt[t] ~ lognormal(log(Rt[t-1]),sigma);
  }
  for(t in sampleday){
    viral[t,1:viral_obs[t]] ~ normal(log(E_viral[t]), nu[t]);
  }
  for(t in 2:N1){
      psi[t] ~ lognormal(log(psi[t-1]),psi_sigma);
  }
  {
    phi ~ normal(0,5);
    /*I_import[1:2]~normal(150,75);
    for(t in 3:N){
      I_import[t]~normal(import[t-2], max(import[t-2]/2,1));
    }*/
    for(i in EpidemicStart:N){
      /*cases[i] ~ neg_binomial_2(E_cases[i],phi[1]);
      import[i] ~ neg_binomial_2(E_import[i],phi[2]);
      deaths[i] ~ neg_binomial_2(E_deaths[i],phi[3]);*/
      hosp[i] ~ neg_binomial_2(E_hosp[i],phi[4]);
    }
  }
}
'
}

#clinical data model 
{
stan_code <- '
data {
  int <lower=1> N0; // number of days for which to impute infections
  int<lower=1> N; // days of observed data 
  int<lower=1> Nd; // days of observed deaths 
  int cases[N]; // reported local case counts
  //int import[N]; // reported imported case counts
  int n_g;
  real g[n_g]; // distribution of time from infection to report 
  int n_SI;
  real SI[n_SI]; // fixed pre-calculated SI 
  int n_h;
  real h[n_h]; // distribution of time from infection to hospitalization 
  int hosp[N]; // reported new hospitalizations
  int n_f;
  real f[n_f]; // distribution of time from infection to death, by age group 
  int deaths[Nd]; // reported deaths, no data for May 2023
  int weekday[N]; // which day of a week 
  int EpidemicStart; // start time for inference
  int tw; //number of days in a short window when reporting rate is constant
  int N1; // number of biweeks
  int biweeks[N]; //biweek index
 }

parameters {
  real<lower=0> kappa; //sd for inital case counts
  real<lower=0> y; //initial case count
  vector<lower=0>[4] phi; //negative binomial parameter
  real<lower=0> tau; //expected initial case count
  vector<lower=0>[N] Rt;
  real<lower=0> R0; //intial Rt
  real<lower=0> sigma; //Rt variation
  //vector<lower=0>[N] I_import; //latent importation
  simplex[7] r11; //scaled relative reporting rate for each week, local cases
  //simplex[7] r12; //scaled relative reporting rate for each week, imported cases
  simplex[7] r21; //scaled relative reporting rate for each week, hospitalised cases
  vector<lower=0>[N1] rh; // case-hospitalisation rates
  real<lower=0> rh_sigma; // sd for variation of rh
  vector<lower=0>[N1] rd; // case-fatality rates
  real<lower=0> rd_sigma; // sd for variation of rd
}

transformed parameters {
    real convolution;
    vector[N] prediction = rep_vector(0,N);
    vector[N] E_cases  = rep_vector(0,N);
    //vector[N] E_import  = rep_vector(0,N);
    vector[N] E_hosp  = rep_vector(0,N);
    vector[N] E_deaths  = rep_vector(0,N);
    //vector<lower=0>[N] E_viral = rep_vector(0,N);
    vector[7] r = r11*7;
    //vector[7] r1 = r12*7;
    vector[7] r2  = r21*7;
    // expected case counts on day t
    prediction[1:N0] = rep_vector(y,N0); // learn the number of cases in the first N0 days
    for (i in (N0+1):N) {
        convolution=0;
        for(j in 1:min(n_SI, i-1)) {
          //convolution += (prediction[i-j]+I_import[i-j])*SI[j]; 
          convolution += prediction[i-j]*SI[j]; 
        }
        prediction[i] = Rt[i] * convolution;
    }
    {
      E_cases[1]=1e-9;
      //E_import[1]=1e-9;
      E_hosp[1]=1e-9;
      E_deaths[1]= 1e-9;
      for (i in 2:N){
        E_cases[i]=0;
        //E_import[i]=0;
        for(j in 1:min(n_g, i-1)){
          E_cases[i] += prediction[i-j]*g[j];
          //E_import[i] += I_import[i-j]*g[j];
        }
        E_cases[i]*=r[weekday[i]];
        //E_import[i]*=r1[weekday[i]];
        E_hosp[i]=0;
        for(j in 1:min(n_h, i-1)){
          E_hosp[i] += prediction[i-j]*h[j];
        }
        E_hosp[i]*=r2[weekday[i]]*rh[biweeks[i]]; //different hospitalization rates
        E_deaths[i]= 0;
        for(j in 1:min(n_f, i-1)){
          E_deaths[i] += prediction[i-j]*f[j]; //different distribution from onset to death for different age groups
        }
        E_deaths[i]*=rd[biweeks[i]];
      }
    }
 }
model {
  rd[1] ~ normal(0,0.1);
  rh[1] ~ normal(0,0.1);
  tau ~ exponential(0.03);
  r11 ~ dirichlet(rep_vector(1, 7));
  //r12 ~ dirichlet(rep_vector(1, 7));
  r21 ~ dirichlet(rep_vector(1, 7));
  y ~ exponential(1.0/tau);
  kappa ~ normal(0,1);
  R0 ~ normal(2, 1);
  Rt[1] ~ normal(R0, kappa); 
  sigma~gamma(0.01,0.01);
  rd_sigma~gamma(0.01,0.01);
  rh_sigma~gamma(0.01,0.01);
  for(t in 2:N){
    Rt[t] ~ lognormal(log(Rt[t-1]),sigma);
  }
  for(t in 2:N1){
    rh[t] ~ lognormal(log(rh[t-1]),rh_sigma);
    rd[t] ~ lognormal(log(rd[t-1]),rd_sigma);
  }
  {
    phi ~ normal(0,5);
    /*I_import[1:2]~normal(150,75);
    for(t in 3:N){
      I_import[t]~normal(import[t-2], max(import[t-2]/2,1));
    }*/
    for(i in EpidemicStart:N){
      cases[i] ~ neg_binomial_2(E_cases[i],phi[1]);
      //import[i] ~ neg_binomial_2(E_import[i],phi[2]);
      hosp[i] ~ neg_binomial_2(E_hosp[i],phi[4]);
    }
    for(i in EpidemicStart:Nd){
      deaths[i] ~ neg_binomial_2(E_deaths[i],phi[3]); 
    }
  }
}
'
}

#full model
{
stan_code <- '
data {
  int <lower=1> N0; // number of days for which to impute infections
  int<lower=1> N; // days of observed data 
  int<lower=1> Nd; // days of observed deaths 
  int cases[N]; // reported local case counts
  //int import[N]; // reported imported case counts
  int n_g;
  real g[n_g]; // distribution of time from infection to report 
  int n_SI;
  real SI[n_SI]; // fixed pre-calculated SI 
  int n_h;
  real h[n_h]; // distribution of time from infection to hospitalization 
  int hosp[N]; // reported new hospitalizations
  int n_f;
  real f[n_f]; // distribution of time from infection to death, by age group 
  int deaths[Nd]; // reported deaths
  int n_v;
  real v[n_v]; // distribution of viral load against time
  int<lower=1> M; // number of total communities included in the samples
  real viral[N, M]; // log viral load in community wastewater
  int <lower=0> viral_obs[N]; //idx for test result (numeric)
  int weekday[N]; // which day of a week 
  int EpidemicStart; // start time for inference
  int tw; //number of days in a short window when reporting rate is constant
  int N1; // number of biweeks
  int biweeks[N]; //biweek index
  int nsd; //number of sampled days
  int sampleday[nsd]; //sampled days
 }

parameters {
  vector<lower=0>[N1] psi; // scaling parameter for cases
  real<lower=0> psi_sigma; // scaling parameter for cases
  vector<lower=0>[N] nu; // sd for LN distribution of viral load
  real<lower=0> nu0; //mean of sd for LN distribution of viral load
  real<lower=0> nu_sigma; // sd of sd for LN distribution of viral load
  real<lower=0> kappa; //sd for inital case counts
  real<lower=0> y; //initial case count
  vector<lower=0>[4] phi; //negative binomial parameter
  real<lower=0> tau; //expected initial case count
  vector<lower=0>[N] Rt;
  real<lower=0> R0; //intial Rt
  real<lower=0> sigma; //Rt variation
  //vector<lower=0>[N] I_import; //latent importation
  simplex[7] r11; //scaled relative reporting rate for each week, local cases
  //simplex[7] r12; //scaled relative reporting rate for each week, imported cases
  simplex[7] r21; //scaled relative reporting rate for each week, hospitalised cases
  vector<lower=0>[N1] rh; // case-hospitalisation rates
  real<lower=0> rh_sigma; // sd for variation of rh
  vector<lower=0>[N1] rd; // case-fatality rates
  real<lower=0> rd_sigma; // sd for variation of rd
}

transformed parameters {
    real convolution;
    vector[N] prediction = rep_vector(0,N);
    vector[N] E_cases  = rep_vector(0,N);
    //vector[N] E_import  = rep_vector(0,N);
    vector[N] E_hosp  = rep_vector(0,N);
    vector[N] E_deaths  = rep_vector(0,N);
    vector<lower=0>[N] E_viral = rep_vector(0,N);
    vector[7] r = r11*7;
    //vector[7] r1 = r12*7;
    vector[7] r2  = r21*7;
    // expected case counts on day t
    prediction[1:N0] = rep_vector(y,N0); // learn the number of cases in the first N0 days
    for (i in (N0+1):N) {
        convolution=0;
        for(j in 1:min(n_SI, i-1)) {
          //convolution += (prediction[i-j]+I_import[i-j]/psi[biweeks[i-j]])*SI[j]; 
          convolution += prediction[i-j]*SI[j]; 
        }
        prediction[i] = Rt[i] * convolution;
    }
    {
      E_cases[1]=1e-9;
      //E_import[1]=1e-9;
      E_hosp[1]=1e-9;
      E_deaths[1]= 1e-9;
      for (i in 2:N){
        E_cases[i]=0;
        //E_import[i]=0;
        for(j in 1:min(n_g, i-1)){
          E_cases[i] += prediction[i-j]*g[j];
          //E_import[i] += I_import[i-j]*g[j];
        }
        E_cases[i]*=r[weekday[i]]*psi[biweeks[i]];
        //E_import[i]*=r1[weekday[i]];
        E_hosp[i]=0;
        for(j in 1:min(n_h, i-1)){
          E_hosp[i] += prediction[i-j]*h[j];
        }
        //E_hosp[i]*=r2[weekday[i]]*psi[biweeks[i]]; //different hospitalization rates
        E_hosp[i]*=r2[weekday[i]]*rh[biweeks[i]]*psi[biweeks[i]]; //different hospitalization rates
        E_deaths[i]= 0;
        for(j in 1:min(n_f, i-1)){
          E_deaths[i] += prediction[i-j]*f[j]; //different distribution from onset to death for different age groups
        }
        E_deaths[i]*=rd[biweeks[i]]*psi[biweeks[i]];
      }
    }
    E_viral[1]= 1e-9;
    for(i in 2:N){
      for(j in 1:min(n_v, i-1)){
          E_viral[i] += prediction[i-j]*v[j];
      }
      //E_viral[i]*=psi;
    }
 }
model {
  rd[1] ~ normal(0,0.1);
  rh[1] ~ normal(0,0.1);
  tau ~ exponential(0.03);
  r11 ~ dirichlet(rep_vector(1, 7));
  /*r12 ~ dirichlet(rep_vector(1, 7));
  r21 ~ dirichlet(rep_vector(1, 7));*/
  y ~ exponential(1.0/tau);
  kappa ~ normal(0,1);
  R0 ~ normal(2, 1);
  Rt[1] ~ normal(R0, kappa); 
  sigma~gamma(0.01,0.01);
  nu0~gamma(0.01,0.01);
  nu_sigma~gamma(0.01,0.01);
  nu~normal(nu0,nu_sigma);
  psi[1]~gamma(0.01,0.01);
  psi_sigma~gamma(0.01,0.01);
  rh_sigma~gamma(0.01,0.01);
  for(t in 2:N){
    Rt[t] ~ lognormal(log(Rt[t-1]),sigma);
  }
  for(t in sampleday){
    viral[t,1:viral_obs[t]] ~ normal(log(E_viral[t]), nu[t]);
  }
  for(t in 2:N1){
    psi[t] ~ lognormal(log(psi[t-1]),psi_sigma);
    rh[t] ~ lognormal(log(rh[t-1]),rh_sigma);
    rd[t] ~ lognormal(log(rd[t-1]),rd_sigma);
  }
  {
    phi ~ normal(0,5);
    /*I_import[1:2]~normal(150,75);
    for(t in 3:N){
      I_import[t]~normal(import[t-2], max(import[t-2]/2,1));
    }*/
    for(i in EpidemicStart:N){
      cases[i] ~ neg_binomial_2(E_cases[i],phi[1]);
      //import[i] ~ neg_binomial_2(E_import[i],phi[2]);
      hosp[i] ~ neg_binomial_2(E_hosp[i],phi[4]);
    }
    for(i in EpidemicStart:Nd){
      deaths[i] ~ neg_binomial_2(E_deaths[i],phi[3]);
    }
  }
}
'
} 

  
#data
  tw=7; N1=ceiling((t2-t1+1)/tw)
  sampleday=which(ww.c.location.na[(t1+1):t2]>0)+1
  data_list_w <- list(
    N0 = 3, N = t2-t1+1, EpidemicStart=10,
    tw= tw, 
    N1 = N1, biweeks=as.vector(do.call('cbind',lapply(1:N1, function(i) rep(i,tw))))[1:(t2-t1+1)],
    weekday=weekday[t1:t2],
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
    Nd = t2-t1+1, #number of days with death counts available
    M=M
  )
 
  
library(rstan)
model <- stan_model(model_code = stan_code)
fit <- sampling(model, data = data_list_w, chains = 1, iter = 1e4)
  