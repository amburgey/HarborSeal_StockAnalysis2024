//
// This Stan program defines a simple state space model of
// density dependent population growth, for a population with
// multiple stocks (adapted Jeff Laake model 2022)
//
// Input data
data {
  int<lower=0> Nstk;              // Number of stocks (sub-populations) 
  int<lower=0> Nyrs;              // Number of years population dynamics
  int<lower=0> Nsrv;              // Number of population surveys (in total, all stocks)  
  array[Nsrv] int<lower=0> yr ;   // Year of each survey
  array[Nsrv] int<lower=0> st ;   // Stock ID of each survey  
  array[Nsrv] int<lower=0> Obs ;  // Observed population count (adjusted) for each survey  
  vector<lower=1>[Nstk] N_min ;   // Minimum count for each stock (used to "ground" log N0 values )  
  vector<lower=1>[Nstk] N_max ;   // Maximum count for each stock (used to "ground" log K values )
  vector<lower=0>[Nsrv] a ;       // Alpha parameter for beta distribution of correction factor informed by Huber et al. and London et al.
  vector<lower=0>[Nsrv] b ;       // Beta parameter for beta distribution of correction factor informed by Huber et al. and London et al.
 }
//
transformed data {
  vector[Nstk] log_N_min = log(N_min) ;  // minimum modeled abundance to help set prior
  vector[Nstk] log_N_max = log(N_max) ;  // maximum modeled abundance to help set prior
}
// The parameters accepted by the model. 
parameters {
  vector<lower=0,upper=1>[Nstk] rmax;
  real<lower=0,upper=10> z;           // growth rate inflection (theta of theta-logistic)
  vector[Nstk] log_K ;                // log K for each stock
  vector<lower=1>[Nstk] log_N0 ;      // initial log abundance for each stock at year 1
  real<lower=-1,upper=5> sigma_K ;     // variation in log K density among stocks (was lower=0,upper=3)
  real<lower=-0.3> sigma_r ;          // stochastic variation in instantaneos growth rate
  matrix[Nstk,Nyrs-1] eps ;           // stochastic effects by stock and year 
  real<lower=0, upper=500> phi ;        // inverse scale param for negative binom: observer error (did take up to 200; lower=0,upper=150)
  vector<lower=0,upper=1>[Nsrv] cf ;  // correction factor for animals not on shore
}
// Transformed and derived parameters
transformed parameters {
  matrix[Nstk,Nyrs] N ;         // estimated abundance
  vector<lower=0>[Nstk] K ;     // abundance at K by stock
  K = exp(log_K);
  for(i in 1:Nstk){
    N[i,1] = exp(log_N0[i]) ;
    for(t in 2:Nyrs){
        real lambda ;
        lambda = exp(rmax[i] * (1 - pow( N[i,t-1] / K[i] , z) ) + sigma_r * eps[i,t-1] ) ;
        N[i,t] = N[i,t-1] * lambda ;
    }
  }
}
// The model to be estimated. 
model {
  // Observed data:  
  // (neg binom distribution for discreet counts: alternatively, for continuous
  // non-integer values, could use log-normal or gamma distribtion )
  for(j in 1:Nsrv){
    cf[j] ~ beta(a[j], b[j]) ;
    Obs[j] ~ neg_binomial_2( N[st[j] , yr[j]] * cf[j] , phi) ;   // observations are the real abundance modified by the correction factor for each observation for each stock in each year
  }
  // Priors:
  rmax ~ beta(0.25*10,10) ;           // weakly informed prior for rmax, biological feasibility
  z ~ gamma(2.5,1.5) ;                 // based on Laake et al. and previous literature
  sigma_K ~ cauchy(0,0.1) ;            // vague shrinkage half-cauchy prior
  sigma_r ~ cauchy(0,0.01) ;           // vague shrinkage half-cauchy prior
  log_N0 ~ normal(log_N_min, 1) ;     // weakly informed prior for N0
  phi ~ cauchy(0,1) ;                // vague shrinkage half-cauchy prior (cauchy(0,2.5))
  // Random effects:
  log_K ~ normal(log_N_max, sigma_K) ;
  for(i in 1:Nstk){
    eps[i] ~ normal(0,1) ;
  }
}
// Post-fitting derived values 
generated quantities{
  vector[Nstk] MNPL ;
  vector[Nstk] OSP ;
  // NOTE: very close approximation: MNPL = pow(K,0.9945) * pow((z + 1),(-1/z)) ;
  for(i in 1:Nstk){
    real GR1 = -1 ;
    real GR2 = 0 ;
    real n = round(0.25 * K[i]) ;
    MNPL[i] = n ;
    while (GR2 > GR1){
      n = n + 1 ;
      GR1 = GR2 ;
      GR2 = n * exp(rmax[i] * (1 - pow( n / K[i] , z) ) ) - n ;
      MNPL[i] = n ;
    }
  }
  for(i in 1:Nstk){
    OSP[i] = N[i,Nyrs]/MNPL[i] ;    //calculate OSP as ratio of abundance in most recent year to MNPL
  }
}
