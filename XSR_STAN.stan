// *************************************************************
//  Example Stan code for a state-space model proposed by 
//    Su (2023) “Management performance evaluation of state-space models for 
//                 Pacific pink salmon stock-recruitment analysis” 
//                 submitted for peer review and publication
// STAN: XSR -- Extended State-Space Stock-Recruitment model
// Dr. Zhenming Su
//  IFR, ANN ARBOR, MI
//  Development started on 12/20/2013
//  Added more priors on 5/28/2022
//  Last revised 12/20/2013: 
// *************************************************************
data {
  int<lower=1>  Y;             // number of years (generations)
  real<lower=1E-10> Sobs[Y];   // escapement or observed spawners
  int<lower=1>  K;             // maturity or return age
  real<lower=1E-10> Cobs[Y-K]; // observed catch  
}
parameters{
  real a_t[Y-K+1]; // Ricker stock-recruitment model parameter: a = log(alpha)
  real<lower=0> sd_a;
  real beta;
  real<lower=0, upper=20> sd_S;  // observation error sd of log(Sobs)
  real<lower=0, upper=1> sd_C;    // observation error sd of log(Cobs)

  real<lower=0> sd_R;  // process error sd of log(R)
  real<lower=0> sd_lamb;   // process error sd of lambda
  real lambda[Y-K+1];      // state variable: logit of harvest rate
  real<lower=0> R[Y-K];    // state variable: true unknown recruitment
  real<lower=0> S0[K];
}
transformed parameters{
  real logR[Y-K]; 
  real HR[Y-K+1]; 
  real C[Y-K]; 
  real S[Y]; 

  real r_pred;

  // first lag years
  for (t in 1:K)
  {
	 S[t] = S0[t];
  }
  
  HR[1] = inv_logit(lambda[1]);  //harvest rate at t = 1

  for (t in 1:(Y-K)) { 
	  logR[t] = a_t[t] + log(S[t]) - beta * S[t];
    C[t] = R[t] * HR[t]; 
	  S[t+K] = R[t] - C[t];

    HR[t+1] = inv_logit(lambda[t+1]);
  }
}
model 
{
  real sd_r;
  // first lag years
  for (t in 1:K)
  {
     //S0[t] ~ lognormal(0.0, 100);
	 S0[t] ~ lognormal(0.0, 100);
	 Sobs[t] ~ lognormal(log(S[t]), sd_S);
  }
	
  sd_R ~ cauchy(0, 10);
  sd_a ~ cauchy(0, 10);
  sd_lamb ~ cauchy(0, 10);

  // lambda[1]
  lambda[1] ~ normal(0.0, 2); //STAN uses N(mu, sd)
  a_t[1] ~ normal(0, 100);
  for (t in 1:(Y-K)) { 
    R[t] ~ lognormal(logR[t], sd_R);  
	  a_t[t+1] ~ normal(a_t[t], sd_a);  
    lambda[t+1] ~ normal(lambda[t], sd_lamb);
    Sobs[t+K] ~ lognormal(log(S[t+K]), sd_S);      
    Cobs[t] ~ lognormal(log(C[t]), sd_C);      
  } 
}
generated quantities { 
 // real alpha;  // Ricker alpha: alpha = exp(a):
  real sig2_E; // observation error variance for log(Sobs)
  real sig2_C; // observation variance for log(C)
  real tau2_R; // process error variance for log(R)
  real tau2_a; // process error variance for log(R)
  real tau2_lambda; // process error variance for lambda = logit(HR)
  real SMSY; real hMSY; real MSY;
   
  //  alpha = exp(a);
  sig2_E = sd_S * sd_S;
  sig2_C = sd_C * sd_C;
  tau2_R = sd_R * sd_R;   
  tau2_a = sd_a * sd_a;
  tau2_lambda = sd_lamb * sd_lamb;

 // Management-oriented quantities
 // SMSY = a/beta*(0.5-0.07*a);
 // hMSY = 0.5*a-0.07*a*a;
 // MSY  = alpha*SMSY*exp(-beta*SMSY)-SMSY;
 // print("  alpha=", alpha, "  beta=", beta, "  SMSY=", SMSY, "  hMSY=", hMSY);
}
