unit StateSpace_StockRecruitment_Model_XSR;
//*****************************************************
//  Code used in the MCMC sampling algorithm for a state-space model proposed by 
//    Su (2023) “Management performance evaluation of state-space models for 
//                 Pacific pink salmon stock-recruitment analysis” 
//                 submitted for peer review and publication
//  (c) Zhenming Su, 2023
//*****************************************************

implementation
//*********************************************************
function logSR(S :double; a :double; beta :double) :double;
begin
  //Ricker Model: logR = log(alpha*S*exp(-beta*S)) = log(alpha)+log(S) - beta*S
   result := a + ln(S) - beta * S;
end;

function sigma2E_IG_update(ny :integer; const Et, St :TVector):double;
var
  ssq, c_shape, c_inv_scale :double;
  t :integer;
begin
  // Posterior for the escapement obser. error variance
  // based on a IG prior for sig_E^2:
  //   Gamma density function
  //   f(x) = (rate^shape)/Gamma(shape) * x^(shape-1) exp(-rate * x)
  // IG(alpha, beta) = beta/Gamma(alpha)

  // Rate: inv_scale, beta
  c_inv_scale := 0.001;
  // Shape: alpha
  c_shape := 0.001;
    
  ssq := 0;
  for t := 0 to ny-1 do
  begin
    ssq := ssq + sqr(ln(Et[t]/St[t]));
  end;
  
  result := (c_inv_scale + 0.5*ssq)/Random_Gamma(c_shape + 0.5*ny, MRNGRandom);
end;

function sigma2E_unif_update(ny :integer; const Et, St :TVector) :double;
var
  ssq, rv :double;
  t :integer;
begin
  // Posterior for the escapement obser. error standard deviation (sig_E)
  //   based on a uniform prior for sig_E.
  ssq := 0;
  for t :=  0 to ny-1 do
  begin
    ssq := ssq + sqr(ln(Et[t]/St[t]));
  end;
  
  rv := (0.5*ssq)/Random_Gamma(0.5*(ny-1), MRNGRandom);
  result := rv;
end;

//*********************************************************
function sigma2C_IG_update(ny, return_age :integer; const Cobs, Ct :TVector):double;
var
  ssq, c_inv_scale, c_shape :double;// m, v
  t :integer;
begin
  // Posterior for the catch obser. error variance (sig_C^2)
  //   based on a IG prior for sig_C^2
  // Gamma density function
  // f(x) = (rate^shape)/Gamma(shape) * x^(shape-1) exp(-rate * x)
  // IG(alpha, beta) = beta/Gamma(alpha)
  // rate: inv_scale, beta
  c_inv_scale := 0.04;
  // shape: alpha
  c_shape := 4;
  
  ssq := 0;
  for t := return_age to ny-1 do
  begin
    ssq := ssq + sqr(ln(Cobs[t]/Ct[t]));
  end;
  result := (c_inv_scale + 0.5*ssq)/Random_Gamma(c_shape + 0.5*(ny-return_age), MRNGRandom);
end;

function sigma2C_unif_update(ny, return_age :integer; const Cobs, Ct :TVector) :double;
var
  ssq, rv :double;
  t :integer;
begin
  // uniform prior for sig_C
  ssq := 0;
  for t :=  return_age to ny-1 do
  begin
    ssq := ssq + sqr(ln(Cobs[t]/Ct[t]));
  end;
  
  rv := (0.5*ssq)/Random_Gamma(0.5*(ny-(return_age+1)), MRNGRandom);
  result := rv
end;

function tau2R_IG_update(ny, return_age :integer;
                const St, Rt :TVector; const a_t :TVector; b :double) :double;
var
  ssq, c_inv_scale, c_shape :double;
  t :integer;
begin
  // IG prior for the recruitment process variance tau^2_R
  // IG prior for tau^2_R
  // Gamma density function
  // f(x) = (rate^shape)/Gamma(shape) * x^(shape-1) exp(-rate * x)
  // IG(alpha, beta) = beta/Gamma(alpha)

  // informative prior
  // rate: inv_scale, beta
  c_inv_scale := 0.390625;  
  // shape: alpha
  c_shape := 1.5625;
  
  ssq := 0;
  //years return_age..ny
  for t := return_age to ny-1 do
  begin
    ssq := ssq + sqr(ln(Rt[t])-logSR(St[t-return_age], a_t[t-return_age], b));
  end;
  result := (c_inv_scale + 0.5*ssq)/Random_Gamma(c_shape + 0.5*(ny-return_age), MRNGRandom);
end;

function tau2R_unif_update(ny, return_age :integer; const St, Rt, a_t :TVector;
                 b :double) :double;
var
  ssq :double;
  t :integer;
begin
  //Uniform prior for the sd of recruitment error: tau_R

  ssq := 0;
  //years return_age..ny
  for t := return_age to ny-1 do
  begin
    ssq := ssq + sqr(ln(Rt[t])-logSR(St[t-return_age], a_t[t-return_age], b));
  end;
  
  result := (0.5*ssq)/Random_Gamma((0.5*(ny-(return_age+1))), MRNGRandom);
end;

function tau2lamb_IG_update(ny, return_age :integer;
                  const lambda_t :TVector) :double;
var
  ssq :double;
  t :integer;
begin
  // IG prior for inv-logit(harvest rate) process variance tau^2
  ssq := 0;
  //years: [return_age+2, ny]
  for t := return_age+1 to ny-1 do
  begin
    ssq := ssq + sqr(lambda_t[t] - lambda_t[t-1]);
  end;
  result := (c_IG_Prior + 0.5*ssq)/Random_Gamma(c_IG_Prior + 0.5*(ny-return_age-1), MRNGRandom);
end;

function tau2lamb_unif_update(ny, return_age :integer; const lambda_t :TVector) :double;
var
  ssq :double;
  t :integer;
begin
  //uniform prior for inv-logit(harvest rate) process sd tau
  ssq := 0;
  //years: [return_age+2, ny]
  for t := return_age+1 to ny-1 do
  begin
    ssq := ssq + sqr(lambda_t[t] - lambda_t[t-1]);
  end;
  result :=  (0.5*ssq)/Random_Gamma((0.5*(ny-(return_age+2))), MRNGRandom);
end;

function tau2lamb_IG_HM_update(ny, return_age :integer; mu_lamb :double;
                      const lambda_t :TVector) :double;
var
  ssq :double;
  t :integer;
begin
  // IG prior for inv-logit(harvest rate) process variance tau^2
  ssq := 0;
  //years: [return_age+2, ny]
  for t := return_age to ny-1 do
  begin
    ssq := ssq + sqr(lambda_t[t] - mu_lamb);
  end;
  result := (c_IG_Prior + 0.5*ssq)/Random_Gamma(c_IG_Prior + 0.5*(ny-return_age), MRNGRandom);
end;

function tau2lamb_unif_HM_update(ny, return_age :integer;
                 mu_lamb :double; const lambda_t :TVector) :double;
var
  ssq :double;
  t :integer;
begin
  // Uniform prior for inv-logit(harvest rate) process error sd tau_lambda
  ssq := 0;
  //years: [return_age+2, ny]
  for t := return_age to ny-1 do
  begin
    ssq := ssq + sqr(lambda_t[t] - mu_lamb);
  end;
  result :=  (0.5*ssq)/Random_Gamma((0.5*(ny-(return_age+1))), MRNGRandom);
end;

function tau2a_IG_update(ny, return_age :integer;
                           const a_t :TVector) :double;
var
  ssq :double;
  t :integer;
begin
  // IG prior for RW Ricker 'a' process variance tau^2
  //years: [2, ny-return_age]
  ssq := 0;
  for t := 1 to ny-return_age-1 do
  begin
    ssq := ssq + sqr(a_t[t] - a_t[t-1]);
  end;
  result := (c_IG_Prior+0.5*ssq)/Random_Gamma(c_IG_Prior+0.5*(ny-return_age-1), MRNGRandom);
end;

function tau2a_unif_update(ny, return_age :integer; const a_t :TVector) :double;
var
  ssq :double;
  t :integer;
begin
  //uniform prior for RW Ricker 'a' process process sd
  //years: [2, ny-return_age]
  ssq := 0;
  for t := 1 to ny-return_age-1 do
  begin
    ssq := ssq + sqr(a_t[t] - a_t[t-1]);
  end;
  result := (0.5*ssq)/Random_Gamma((0.5*(ny-return_age-2)), MRNGRandom);
end;

// Random walk Ricker log productivity para 'a_t'
procedure at_update(return_age, ny :integer;
              beta, tau2_r, tau2_a :double; var a_t :TVector;
              const Rt, St :TVector);
var
  mu, V, mu_a, V_a, InvVn :double;
  t :integer;
begin
  for t := 0 to ny-return_age-1 do
  begin
    //First brood year
    if (t=0) then
    begin
      mu_a := 1.5;
      V_a  := 4.0;
      
      InvVn := 1/(tau2_a*tau2_r + V_a*tau2_r + V_a*tau2_a);
      mu := (mu_a*tau2_a*tau2_r + a_t[1]*V_a*tau2_r
              + (ln(Rt[2])-ln(St[0])+beta*St[0])*V_a*tau2_a
            )*InvVn;
      V :=  V_a*tau2_a*tau2_r * InvVn;
      a_t[t] := Random_Normal(mu, sqrt(V), MRNGRandom);
    end;

    //years (return_age+1) to ny-return_age-1
    if ((t >= 1) and (t < ny-return_age-1)) then
    begin
      InvVn := 1/(tau2_a+2*tau2_r);
      mu :=
          ((ln(Rt[t+return_age])-ln(St[t])+beta*St[t])*tau2_a
              + (a_t[t+1]+a_t[t-1])*tau2_r
           ) * InvVn;
      V :=  tau2_r*tau2_a * InvVn;
      a_t[t] := Random_Normal(mu, sqrt(V), MRNGRandom);
    end;

    //year = ny-return_age
    if (t = ny-return_age-1) then
    begin
      InvVn := 1/(tau2_a+tau2_r);
      mu :=
          ((ln(Rt[t+return_age])-ln(St[t])+beta*St[t])*tau2_a
              + a_t[t-1]*tau2_r
           ) * InvVn;
      V :=  tau2_r*tau2_a * InvVn;
      a_t[t] := Random_Normal(mu, sqrt(V), MRNGRandom);
    end;
  end;
end;

function a_update(ny, return_age :integer;
              a, beta, tau2 :double;
              const Rt, St :TVector):double;
var
  mu2, V2, mu, V, ssq, mu_a, V_a :double;
  t :integer;
begin
  // const 'a'
  //first return_age years
  ssq := 0;
  //years (return_age+1) to ny
  for t := return_age to ny-1 do
  begin
    ssq := ssq+(ln(Rt[t])-(logSR(St[t-return_age], a, beta) - a));
  end;

  // prior
  mu_a := 0;
  V_a := 1000;
  // likelihood
  mu2 := ssq/(ny-return_age);
  V2  := tau2/(ny-return_age);
  // posterior
  mu := (mu_a * V2 + mu2 * V_a)/(V_a + V2);
  V :=  V_a * V2/(V_a + V2);
  result := Random_Normal(mu, sqrt(V), MRNGRandom);
end;

function b_Ricker_update(ny, return_age :integer; a_t :TVector;
          beta, tau2_R :double; const St, Rt :TVector) :double;
var
  sum0, sum1 :double;
  mu_b, V_b, mu2, V2, mu, V, rv_b :double;
  t :integer;
begin
  sum0 := 0;
  sum1 := 0;
  //other years
  for t := return_age to ny-1 do
  begin
    sum0 := sum0 + St[t-return_age]
                 * (a_t[t-return_age]+ln(St[t-return_age])-ln(Rt[t]));
    sum1 := sum1 + St[t-return_age]*St[t-return_age];
  end;

  //normal prior parameter
  mu_b := 0;
  V_b  := 10;
  
  // posterior
  mu2 := sum0/sum1;
  V2 := tau2_R/sum1;
  mu := (mu_b*V2+mu2*V_b)/(V_b+V2);
  V := V_b*V2/(V_b+V2);
  rv_b := Random_Normal(mu, sqrt(V), MRNGRandom);
  result := rv_b;//Random_Normal(mu, sqrt(V), MRNGRandom);
end;

function logpost_S1(a, b, tau2, sig2, E1, R2, S1 :double) :double;
var
  s0, V0 :double;
begin
  s0 := 0;
  V0 := 10;
  logpost_S1 := - 0.5*sqr(ln(S1)-s0)/V0    // prior
                - 0.5*sqr(ln(E1/S1))/sig2  // likelihood
                - 0.5*sqr(ln(R2)-logSR(S1, a, b))/tau2; //hierarchical prior
end;

function S1_update(iter :integer; S1, R2,
          a, beta, sig2, tau2 :double; const E1 :double;
            var sd_jump, p_jump, SumRate, meanRA :double;
                var tune_done, Adapting :boolean) :double;
var
  log_post_old, log_post_new, ratio,
  S1star :double;
begin
  //draw log(S1): proposal distribution is for log(S1)
  S1star := exp(Random_Normal(ln(S1), sd_jump, MRNGRandom));

  log_post_old := logpost_S1(a, beta,tau2, sig2, E1, R2, S1);
  log_post_new := logpost_S1(a, beta, tau2, sig2, E1, R2, S1star);

  if not ((IsNan(log_post_new))
       or IsInfinite(log_post_new)) then
  begin
    ratio := exp(log_post_new - log_post_old);
    if MRNGRandom < ratio then
      S1 := S1star;
    p_jump := min(ratio,1);
  end
  else
  begin
    p_jump := 0;
  end;

  result := S1;

  if Adapting then
     Adaptive_Tuning_1D(iter, p_jump, SumRate, sd_jump, meanRA, tune_done);
     //AMwG(iter, p_jump, SumRate, sd_jump, meanRA, tune_done);
end;


function logpost_Rt_ind(t :integer; ny, return_age :integer;
                    att, b, sig2E, sig2C, tau2, Rtt, Rt, Ct,
                    Et, Cobs :double) :double;
var
  logpost, St :double;
begin
  St := Rt - Ct;
  if (t >= ny-return_age) then
    //likehood
    logpost := - 0.5 * sqr(ln(Et/St))/(sig2E)
               - 0.5 * sqr(ln(Cobs/Ct))/(sig2E)
  else
    // likehood
    logpost := - 0.5 * sqr(ln(Et/St))/(sig2E)
               - 0.5 * sqr(ln(Cobs/Ct))/(sig2E)
               - 0.5 * sqr(ln(Rtt) - logSR(St, att, b))/(tau2);  // hierarchcial prior p(r[t+k])

  result := logpost;
end;

function logpost_Rt_Metrop(t :integer; ny, return_age :integer;
                      const a_t :TVector; b, sig2E, sig2C, tau2R :double;
                      const St_k, Rtt :double; Et, Cobs :double;
                      Rt, ht :double):double;
var
  logpost, St, Ct :double;
begin
  Ct := Rt * ht;
  St := Rt - Ct;

  if (t >= ny-return_age) then
  begin
     logpost := - 0.5*sqr(ln(Et/St))/(sig2E)
                - 0.5*sqr(ln(Cobs/Ct))/(sig2C)
                - 0.5*sqr(ln(Rt)-logSR(St_k,a_t[t-return_age],b))/(tau2R);
  end
  else
  begin
     logpost := - 0.5*sqr(ln(Et/St))/(sig2E)
                - 0.5*sqr(ln(Cobs/Ct))/(sig2C)
                - 0.5*sqr(ln(Rt)-logSR(St_k, a_t[t-return_age], b))/(tau2R)
                - 0.5*sqr(ln(Rtt)-logSR(St, a_t[t], b))/(tau2R);
  end;
  result := logpost;
end;

procedure StateUpdate_R_Metrop(iter, ny, return_age :integer; const a_t :TVector;
  b, sig2E, sig2C, tau2R, tau2lamb :double; const Cobs, Et :TVector;
    var Rt, St, Ct, ht :TVector; var sd_jump :TVector; var p_jump :TVector;
       var SumRate, meanRA :TVector;
         var tune_done :array of boolean;
           var Adapting :boolean) ;
var
  log_post_old, log_post_new, ratio,
  Rold, R_star, Rtt :double;
  t :integer;
begin
 for t := return_age to ny-1 do
 begin
  Rold := Rt[t];
  R_star := exp(Random_Normal(ln(Rold), sd_jump[t], MRNGRandom));

  if (t < ny-return_age) then
       Rtt := Rt[t+return_age]
  else
       Rtt := 0;
       
  log_post_old := logpost_Rt_Metrop(t, ny, return_age, a_t, b,
                    sig2E, sig2C, tau2R, St[t-return_age], Rtt,
                    Et[t], Cobs[t], Rold, ht[t]);
  log_post_new := logpost_Rt_Metrop(t, ny, return_age, a_t, b,
                    sig2E, sig2C, tau2R, St[t-return_age], Rtt,
                    Et[t], Cobs[t], R_star, ht[t]);

  if not ((IsNan(log_post_new))
       or IsInfinite(log_post_new)) then
  begin
    ratio := exp(log_post_new - log_post_old);
    if MRNGRandom < ratio then
    begin
      Rt[t] := R_star;
      Ct[t] := Rt[t]*ht[t];
      St[t] := Rt[t]-Ct[t];
    end;
    p_jump[t] := min(ratio,1);
  end
  else
  begin
    p_jump[t] := 0;
  end;

  if Adapting then
       Adaptive_Tuning_1D(iter, p_jump[t], SumRate[t],
          sd_jump[t], meanRA[t], tune_done[t]);
      // AMwG(iter, p_jump[t], SumRate[t], sd_jump[t], meanRA[t], tune_done[t]);
 end; // t
end;

function logpost_Lamb_t_Metrop(t, return_age :integer; ny :integer;
                           sig2E, sig2C, tau2Lamb :double;
                           const Et, Cobs :double; Rt :double;
                           Lamb_m1, Lamb_p1, Lamb_star :double) :double;
var
  logpost, St, ht, Ct, V1 :double;
begin
  ht := exp(Lamb_star)/(1 + exp(Lamb_star));
  Ct := Rt * ht;
  St := Rt - Ct;

  if (t = ny-1) then
  begin
     logpost := - 0.5*sqr(ln(Et/St))/(sig2E) //- ln(Rt[t])
                - 0.5*sqr(ln(Cobs/Ct))/(sig2C)
                - 0.5*sqr(Lamb_star-Lamb_m1)/(tau2Lamb);
  end
  else
  begin
    if (t = return_age) then
    begin
     // Initial lambda
     // lambda[1] ~ N(0, V1)
     V1 := 2;
     logpost := - 0.5*sqr(Lamb_star)/V1 // prior of lambda[1]
                - 0.5*sqr(ln(Et/St))/(sig2E) //- ln(Rt[t])
                - 0.5*sqr(ln(Cobs/Ct))/(sig2C)
                - 0.5*sqr(Lamb_p1-Lamb_star)/(tau2Lamb);

    end
    else
    begin
     logpost := - 0.5*sqr(ln(Et/St))/(sig2E) //- ln(Rt[t])
                - 0.5*sqr(ln(Cobs/Ct))/(sig2C)
                - 0.5*sqr(Lamb_star-Lamb_m1)/(tau2Lamb)
                - 0.5*sqr(Lamb_p1-Lamb_star)/(tau2Lamb);
    end;
  end;
  result := logpost;
end;

function logpost_Lamb_t_Metrop_FixedEff(t, return_age :integer; ny :integer;
                           sig2E, sig2C :double;
                           const Et, Cobs :double; Rt :double;
                           Lamb_star :double) :double;
var
  logpost, St, ht, Ct, mu, V :double;
begin
  ht := exp(Lamb_star)/(1 + exp(Lamb_star));
  Ct := Rt * ht;
  St := Rt - Ct;

  mu := 0;
  V  := 2;
  logpost := - 0.5*sqr(ln(Et/St))/(sig2E) //- ln(Rt[t])
             - 0.5*sqr(ln(Cobs/Ct))/(sig2C)
             - 0.5*sqr(Lamb_star-mu)/V;

  result := logpost;
end;

function logpost_Lamb_t_Metrop_HM(t, return_age :integer; ny :integer;
                           sig2E, sig2C, tau2Lamb :double;
                           const Et, Cobs :double; Rt :double;
                           Lamb_star, Lamb_mu :double) :double;
var
  logpost, St, ht, Ct :double;
begin
  ht := exp(Lamb_star)/(1 + exp(Lamb_star));
  Ct := Rt * ht;
  St := Rt - Ct;

  logpost := - 0.5*sqr(ln(Et/St))/(sig2E) //- ln(Rt[t])
             - 0.5*sqr(ln(Cobs/Ct))/(sig2C)
             - 0.5*sqr(Lamb_star-Lamb_mu)/tau2Lamb;

  result := logpost;
end;

procedure StateUpdate_Lamb_t_Metrop(iter, ny, return_age :integer;
       sig2E, sig2C, tau2R, tau2lamb, mu_lambda :double;
       const Cobs, Et :TVector;
       var Lamb_t, Rt, St, Ct, ht :TVector;
       var sd_jump :TVector; var p_jump :TVector;
       var SumRate, meanRA :TVector;
       var tune_done :array of boolean;
       var Adapting :boolean);
var
  log_post_old, log_post_new, ratio,
  lamb_old, lamb_star, Lamb_p1 :double;
  t :integer;
begin
 for t := return_age to ny-1 do
 begin
  lamb_old := Lamb_t[t];
  lamb_star := Random_Normal(lamb_old, sd_jump[t], MRNGRandom);

  if (lamb_star <= -20) or (lamb_star > 20) then
     p_jump[t] := 0
  else
  begin
    if t < ny-1 then
      Lamb_p1 := Lamb_t[t+1]
    else
      Lamb_p1 := 0;

     if FormSRsimu.RadioButtonLambdaFixedEff.Checked then
    begin
       log_post_old := logpost_Lamb_t_Metrop_FixedEff(t, return_age, ny,
                         sig2E, sig2C, Et[t], Cobs[t], Rt[t],
                         lamb_old);

       log_post_new := logpost_Lamb_t_Metrop_FixedEff(t, return_age, ny,
                         sig2E, sig2C, Et[t], Cobs[t], Rt[t],
                         lamb_star);
    end
    else if FormSRsimu.RadioButtonLambdaPriorHM.Checked then
    begin
       log_post_old := logpost_Lamb_t_Metrop_HM(t, return_age, ny,
                         sig2E, sig2C, tau2lamb,
                         Et[t], Cobs[t], Rt[t],
                         lamb_old, mu_lambda);

       log_post_new := logpost_Lamb_t_Metrop_HM(t, return_age, ny,
                         sig2E, sig2C, tau2lamb,
                         Et[t], Cobs[t], Rt[t],
                         lamb_star, mu_lambda);
    end
    else
    begin  // random walk
       log_post_old := logpost_Lamb_t_Metrop(t, return_age, ny,
                         sig2E, sig2C, tau2lamb,
                         Et[t], Cobs[t], Rt[t],
                         Lamb_t[t-1], Lamb_p1,
                         lamb_old);

       log_post_new := logpost_Lamb_t_Metrop(t, return_age, ny,
                         sig2E, sig2C, tau2lamb,
                         Et[t], Cobs[t], Rt[t],
                         Lamb_t[t-1], Lamb_p1,
                         lamb_star);
    end;

    if not ((IsNan(log_post_new))
         or IsInfinite(log_post_new)) then
    begin
      ratio := exp(log_post_new - log_post_old);
      if MRNGRandom < ratio then
      begin
         Lamb_t[t] := lamb_star;
         ht[t] := exp(Lamb_t[t])/(1 + exp(Lamb_t[t]));
         Ct[t] := Rt[t]*ht[t];
         St[t] := Rt[t]-Ct[t];
      end;
      p_jump[t] := min(ratio, 1);
    end
    else
    begin
      Lamb_t[t] := lamb_old;
      p_jump[t] := 0;
    end;
  end;

  if Adapting then
       Adaptive_Tuning_1D(iter, p_jump[t], SumRate[t],
          sd_jump[t], meanRA[t], tune_done[t]);
      // AMwG(iter, p_jump[t], SumRate[t], sd_jump[t], meanRA[t], tune_done[t]);
 end; // t
end;

function mu_lamb_update(ny, return_age :integer; const lambda :TVector;
             tau2_lamb :double) :double;
var
  sum0 :double;
  mu_1, V_1, mu2, V2, mu, V, rv_b :double;
  t :integer;
begin
  //first year
  sum0 := 0;
  //other years
  for t := return_age to ny-1 do
  begin
    sum0 := sum0 + lambda[t];
  end;

  //normal prior parameter
  mu_1 := 0;
  V_1  := 2;

  // likelihood
  //mu2 := sum0/(T-return_age);
  //V2 := tau2_lamb/(T-return_age);
  // Changed T to ny   March 17, 2019
  mu2 := sum0/(ny-return_age);
  V2 := tau2_lamb/(ny-return_age);
  // posterior
  mu := (mu_1*V2+mu2*V_1)/(V_1+V2);
  V := V_1*V2/(V_1+V2);

  // make a draw of new b
  //repeat
  rv_b := Random_Normal(mu, sqrt(V), MRNGRandom);
  //until rv_b>0;

  result := rv_b;//Random_Normal(mu, sqrt(V), MRNGRandom);
end;

procedure GibbsOneIter(imcmc, ny,
           return_age :integer; RickerAVarying :Boolean;
           var a0, beta, sig2E, sig2C, tau2_a, tau2R,
               tau2lamb, mu_lambda :double;
           const Et, Cobs :TVector;
           var a_t :TVector;
           var Ct, St, Rt, lambda_t :TVector;
           var sd_jumpR, p_jumpR :TVector;
           var sdS1, P_JumpS1 :TVector;
           var SumRateS1, meanRAS1 :TVector;
           var tune_doneS1 :array of boolean;
           var Adapting :boolean;
           var SumRateR, meanRAR :TVector;
           var tune_done_R :array of boolean;
           var sd_jump_lamb :TVector;
           var p_jump_lamb :TVector;
           var SumRate_lamb, meanRA_lamb :TVector;
           var tune_done_lamb :array of boolean;
           var ifail :boolean);
var
  theta, ExPara :TVector;
  t, k :Integer;
  logRt, ht :TVector;
  data :TMatrix;
  rho, a, Sstar :double;
begin
 try
  setlength(ht, ny);
  begin
    for t := return_age to ny-1 do
    begin
      ht[t] := exp(lambda_t[t])/(1 + exp(lambda_t[t]));
      Ct[t] := Rt[t]*ht[t];
      St[t] := Rt[t]-Ct[t];
    end;
  end;

  if FormSRsimu.RadioButtonUnifPriorTau.checked then
     tau2R := tau2R_unif_update(ny, return_age, St, Rt, a_t, beta)
  else
     tau2R := tau2R_IG_update(ny, return_age, St, Rt, a_t, beta);

  if FormSRsimu.RadioButtonLambdaFixedEff.Checked then
     tau2lamb := 0
  else if FormSRsimu.RadioButtonLambdaPriorHM.Checked then
     tau2lamb := tau2lamb_unif_HM_update(ny, return_age, mu_lambda, lambda_t)
  else
     tau2lamb := tau2lamb_unif_update(ny, return_age, lambda_t);

  if RickerAVarying then
      tau2_a := tau2a_unif_update(ny, return_age, a_t)
  else
      tau2_a := 0;

  if not FormSRsimu.CheckBoxFixedTau2E.Checked then
    if FormSRsimu.RadioButtonIGPriorSig2.checked then
       sig2E := sigma2E_IG_update(ny, Et, St)
    else
       sig2E := sigma2E_unif_update(ny, Et, St)
  else
    sig2E := sqr(StrToFloat(FormSRsimu.EditFixedSigma_E_XSR.text));

  if FormSRsimu.CheckBoxFixedTau2C.checked then
     sig2C := sqr(StrToFloat(FormSRsimu.EditFixed_SigmaC_XSR.text)) 
  else
  begin
    if FormSRsimu.RadioButtonUnifPriorSig.checked then
      sig2C := sigma2C_unif_update(ny, return_age, Cobs, Ct)
    else 
      sig2C := sigma2C_IG_update(ny, return_age, Cobs, Ct);
  end;

  if RickerAVarying then  // Time-varying Ricker 'a'
    at_update(return_age, ny, beta, tau2R, tau2_a, a_t, Rt, St)
  else
  begin
    // constant Ricker 'a'
    a := a_update(ny, return_age, a_t[0], beta, tau2R, Rt, St);
    for t := 0 to ny-return_age-1 do
       a_t[t] := a;
  end;

  beta := b_Ricker_update(ny, return_age, a_t, beta, tau2R, St, Rt);

  for k := 0 to return_age-1 do
  begin
    // initial S1 to Sk
    St[k] := S1_update(imcmc, St[k], Rt[k+return_age],
                a_t[k], beta, sig2E, tau2R, Et[k],
                sdS1[k], P_JumpS1[k], SumRateS1[k], meanRAS1[k],
                tune_doneS1[k], Adapting);
  end;

  if FormSRsimu.RadioButtonLambdaPriorHM.Checked then
  begin
     mu_lambda := mu_lamb_update(ny, return_age, lambda_t, tau2lamb)
  end;
  
  StateUpdate_Lamb_t_Metrop(imcmc, ny, return_age,
       sig2E, sig2C, tau2R, tau2lamb, mu_lambda,
       Cobs, Et, lambda_t, Rt, St, Ct, ht,
       sd_jump_lamb, p_jump_lamb,
       SumRate_lamb, meanRA_lamb, tune_done_lamb,
       Adapting);

  for t := return_age to ny-1 do
  begin
      ht[t] := exp(lambda_t[t])/(1 + exp(lambda_t[t]));
      Ct[t] := Rt[t]*ht[t];
      St[t] := Rt[t]-Ct[t];
  end;

  StateUpdate_R_Metrop(imcmc, ny, return_age, a_t, beta,
            sig2E, sig2C, tau2R, tau2lamb,
            Cobs, Et, Rt, St, Ct, ht, sd_jumpR, p_jumpR,
            SumRateR, meanRAR, tune_done_R, Adapting);

  if FormSRsimu.CheckBoxRealTimeES.checked then
  begin
     FormSRsimu.Series1.Clear;
     FormSRsimu.Series2.Clear;
     FormSRsimu.Series3.Clear;
     for t := 0 to ny-return_age-1 do
     begin
       FormSRsimu.Series1.addxy(t+1, Et[t]);
       FormSRsimu.Series2.addxy(t+1, St[t]);
       FormSRsimu.Series3.addxy(t+1, a_t[t]);
     end;
     application.ProcessMessages;
  end;
  if FormSRsimu.CheckBoxRealTimeR.checked then
  begin
     FormSRsimu.Series4.Clear;
     FormSRsimu.Series5.Clear;
     FormSRsimu.Series6.Clear;
     FormSRsimu.Series7.Clear;
     for t := return_age to ny -1 do
     begin
       FormSRsimu.Series4.addxy(t+1, ht[t]);
       FormSRsimu.Series5.addxy(t+1, Cobs[t]);
       FormSRsimu.Series6.addxy(t+1, Ct[t]);
       FormSRsimu.Series7.addxy(t+1, Rt[t]);
     end;
     application.ProcessMessages;
  end;
 finally
   raise Exception.Create('Error occured') at @GibbsOneIter;
 end;
end;

procedure Gibbs_sampling(niter, NBurnIn, n_thin, return_age :integer;
            const para :TVector; const Et :TVector; const Cobs, RTrue :TVector;
            var theta_matrix, a_last_matrix :TMatrix;
            var Recruit_Pred :TMatrix; var ifail :boolean);
var
   ny, i, t, niters, k :Integer;
   RickerAVarying :boolean;
   a, beta, sig2E, sig2C, tau2a, tau2R, tau2lamb, mu_lambda,
    Eint, MSE, hr :double;
   at_pred, logR, RPred :double;
   Rt, St, Ct, a_t, lambda_t :TVector;
   ExPara :TVector;
   dataEmpty :TMatrix;
   SdS1, P_JumpS1, SumRateS1, meanRAS1 :TVector; // initial states
   tune_doneS1:array of boolean;
   Adapting :boolean;
   All_doneRt, All_doneS1, All_done_lamb :integer;
   sd_jumpR :TVector;
   SumRateR, meanRAR :TVector;
   tune_done_R :array of boolean;
   p_jumpR :TVector;

   sd_jump_lamb :TVector;
   p_jump_lamb :TVector;
   SumRate_lamb, meanRA_lamb :TVector;
   tune_done_lamb :array of boolean;

   n_tuning_interval :integer;

   tuningscale :double;
   V_pd :TVector;
   sumx :double;
   str :string;
begin
//a = log(alpha), beta, sig2E, sig2C, tau2a, tau2R, tau2Lamb
 a    := para[0];
 beta := para[1];
 sig2E := para[2];
 sig2C := para[3];
 tau2a := para[4];
 tau2R := para[5];
 tau2lamb:= para[6];
 if FormSRsimu.RadioButtonLambdaPriorHM.Checked then
   mu_lambda := para[7];
 
 //initialize state variables
 ny := length(Et);
 setlength(Rt, ny);
 setlength(St, ny);
 setlength(Ct, ny);
 setlength(lambda_t, ny);
 setlength(a_t, ny-return_age);

 for t := 0 to ny-1 do
 begin
   Eint := Et[t];
   if (Cobs[t]<=0) and (t>return_age) then
      Cobs[t] := 0.0001;
   if Et[t]<=0 then
      Et[t] := 0.0001;

   if Eint <= 0 then
        Eint := Et[t];
   Rt[t] := Eint + Cobs[t];
   St[t] := Eint;
   Ct[t] := Cobs[t];
 end;
 for k := 0 to return_age-1 do
 begin
    Rt[k] := 0;
    Ct[k] := 0;
 end;

 for k := return_age to ny-1 do
 begin
    hr := Ct[k]/Rt[k];
    if hr >= 1.0 then
      hr := 0.9
    else if hr <= 0.0 then
      hr := 0.01;
      
    lambda_t[k] := ln(hr/(1-hr))+Random_Normal(0, sqrt(0.01), MRNGRandom); //RandG(0, 0.01);;
 end;

 if FormSRsimu.CheckBoxRicker_a_const.checked then
 begin
   RickerAVarying := FALSE;
   tau2a := 0;
   // Constant Ricker 'a'
   for k := 0 to ny-return_age-1 do
   begin
      a_t[k] := a;
   end;
 end
 else
 begin
   // time varying Ricker 'a'
   RickerAVarying := TRUE;
   for k := 0 to ny-return_age-1 do
   begin
      a_t[k] := Random_Normal(a, 0.1, MRNGRandom);;
   end;
 end;
 //****************************************
 // Init adaptive tuning parameters for
 //  Metropolis samplers
 //****************************************
 n_tuning_interval := StrToInt(FormSRsimu.EditAdaptInterval.text);

 i := 0;  Adapting := True;
 FormSRsimu.Label_Iter.Caption := 'Adaptive';
 FormSRsimu.LabelCurrEst.Caption := 'Mean Accept Rates';
 
 // S1
 setlength(SdS1, return_age);
 setlength(SumRateS1, return_age);
 setlength(tune_doneS1, return_age);
 setlength(P_JumpS1, return_age);
 setlength(meanRAS1, return_age);
 for k := 0 to return_age-1 do
 begin
     SdS1[k] := abs(0.4 * ln(St[k]));
     SumRateS1[k] := 0;
     meanRAS1[k] := 0;
     P_JumpS1[k] := 0;
     tune_doneS1[k] := false;
 end;

 setlength(SumRate_lamb, ny);
 setlength(p_jump_lamb, ny);
 setlength(tune_done_lamb, ny);
 setlength(meanRA_lamb, ny);
 setlength(sd_jump_lamb, ny);
 for t := 0 to ny-1 do
 begin
    SumRate_lamb[t] := 0;
    tune_done_lamb[t] := false;
    sd_jump_lamb[t] := 0.1;
    SumRate_lamb[t] := 0;
    meanRA_lamb[t] := 0;
    P_Jump_lamb[t] := 0;
    tune_done_lamb[t] := false;
 end;

 // --------------------------------------
 // use normal proposal density to update Rt
 // --------------------------------------
 // only Rt and S1-Sk need to be updated by
 //  Metropolis step and need to be tuned
  setlength(SumRateR, ny);
  setlength(p_jumpR, ny);
  setlength(tune_done_R, ny);
  setlength(meanRAR, ny);
  for t := 0 to ny-1 do
  begin
    SumRateR[t] := 0;
    tune_done_R[t] := false;
  end;
  sd_jumpR := scale(0.3, Rt);

  All_doneRt := 0;
  All_doneS1 := 0;
  All_done_lamb := 0;
  while ((All_doneRt = 0) or (All_done_lamb = 0) or (All_doneS1 = 0)) do
  begin
      GibbsOneIter(i, ny, return_age, RickerAVarying,
        a, beta, sig2E, sig2C, tau2a, tau2R, tau2lamb, mu_lambda,
        Et, Cobs, a_t, Ct, St, Rt, lambda_t, sd_jumpR, p_jumpR,
        SdS1, P_JumpS1, SumRateS1, meanRAS1, tune_doneS1, Adapting,
        SumRateR, meanRAR, tune_done_R,
        sd_jump_lamb, p_jump_lamb, SumRate_lamb,
        meanRA_lamb, tune_done_lamb, ifail);

      for t := 0 to return_age-1 do
      begin
        if tune_doneS1[t] then
              All_doneS1 := All_doneS1 + 1;
      end;
      if All_doneS1 = (return_age) then
        All_doneS1 := 1
      else
        All_doneS1 := 0;

      for t := return_age to ny-1 do
      begin
          if tune_done_R[t] then
              All_doneRt := All_doneRt + 1;
      end;
      if All_doneRt = (ny-return_age) then
         All_doneRt := 1
      else
         All_doneRt := 0;

      for t := return_age to ny-1 do
      begin
            if tune_done_lamb[t] then
                All_done_lamb := All_done_lamb + 1;
      end;
      if All_done_lamb = (ny-return_age) then
           All_done_lamb := 1
      else
           All_done_lamb := 0;

      application.ProcessMessages;
      FormSRsimu.EditGibbsCurrIter.text := inttostr(i);

      if i > 40000 then
      begin
         ifail := true;
         exit;
      end;
      inc(i);
 end; //Next iter

 ifail := false;
 FormSRsimu.Label_Iter.Caption := 'Burn in';
 FormSRsimu.LabelCurrEst.Caption := 'Estimates';
 niters := 0;
 Adapting := true;
 for i := 0 to niter-1 do
 begin
   GibbsOneIter(i, ny, return_age, RickerAVarying,
        a, beta, sig2E, sig2C, tau2a, tau2R, tau2lamb, mu_lambda,
        Et, Cobs, a_t, Ct, St, Rt, lambda_t,
        sd_jumpR, p_jumpR, SdS1, P_JumpS1,
        SumRateS1, meanRAS1, tune_doneS1, Adapting,
        SumRateR, meanRAR, tune_done_R,
        sd_jump_lamb, p_jump_lamb,SumRate_lamb,
        meanRA_lamb,tune_done_lamb,ifail);

  if ifail then break;

  if i >= NBurnIn then
  begin
   if (i mod n_thin) = 0 then
   begin
     //a = log(alpha), beta, sig2E, sig2C, tau2a, tau2R, tau2Lamb
     // Time-varying Ricker 'a'
      if RickerAVarying then
      begin
          sumx := 0;
          for t := 0 to ny-return_age-1 do
             sumx := a_t[t] + sumx;
          theta_matrix[niters,0] := exp(sumx/(ny-return_age));
      end
      else
          theta_matrix[niters,0] := exp(a_t[0]);

      // Calculate predicted recruitment for year ny+1
      //   from St[ny-2] in year ny-1
      at_pred := Random_Normal(a_t[ny-return_age-1], sqrt(tau2a), MRNGRandom);
      logR := (at_pred + ln(St[ny-2]) - beta * St[ny-2]);
      Rpred := exp(Random_Normal(logR, sqrt(tau2R), MRNGRandom));    
      Recruit_Pred[niters,0] := Rpred;

      theta_matrix[niters,1] := beta;
      theta_matrix[niters,2] := sig2E;
      theta_matrix[niters,3] := sig2C;
      theta_matrix[niters,4] := tau2a;
      theta_matrix[niters,5] := tau2R;
      theta_matrix[niters,6] := tau2Lamb;
      if FormSRsimu.RadioButtonLambdaPriorHM.Checked then
        theta_matrix[niters,7] := mu_lambda;

      if IsNan(a) or IsInfinite(a) or
         IsNan(beta) or IsInfinite(beta) or
         IsNan(sig2E) or IsInfinite(sig2E) or
         IsNan(tau2R) or IsInfinite(tau2R) or
         IsNan(St[0]) or IsInfinite(St[0])
      then
      begin
         ifail := true;
      end;
      a_last_matrix[niters,0] := a_t[ny-return_age-1];

      inc(niters);
   end;
  end;
 end; //Next iter
end;

procedure Ricker_pred(const para, Ct :TVector; var St :TVector);
var
  t :integer;
  ny :integer;
begin
  //first year
  //other years
  ny := length(Ct);
  setlength(St, ny);
  St[0] := para[4];
  for t := 1 to ny-1 do
  begin
     St[t] := GeneralSR(St[t-1], 0, para[0], para[1]) - Ct[t];
  end;
end;

procedure Gibbs_ExpStateSpaceSR_Est(simu, return_age_i :integer; const E, Cobs, R, Robs :TVector;
                          ainit, binit, sig2init, tau2init :double; var ParaEst: TVector);
var
  i, nrow, np0, np, NRef,//ny,
  NGibbs, nGibbsRows,
  NBurnIn, n_thin :integer;
  Para :TVector;

  NamesPara : TStrVector;
  theta_matrix, a_last_matrix, Recruit_Pred,
  ElevenNum, ElevenNum_a, ElevenNum_Rpred :TMatrix;
  aestimate, bestimate, hMSY, SMSY :double;
  PropAneg, PropBneg :double;
  ifail :boolean;
  header, res_str :string;
begin
     return_age := return_age_i;

     np0 := 7; //a = log(alpha), beta, sig2E, sig2C, tau2a, tau2R, tau2Lamb
     if FormSRsimu.RadioButtonLambdaPriorHM.Checked then
     begin
       np0 := 8;
     end;
     NRef := 0; //reference points

     c_IG_Prior := StrToFloat(FormSRsimu.EditIG_const_sig2.text);

     np := np0 + NRef + 1;

     //********************************************************
     // Simulations
     //********************************************************
     NGibbs := strtoint(FormSRsimu.Edit_N_Gibbs_Iter.text);
     NBurnIn :=  strtoint(FormSRsimu.EditBurn_In.text);
     n_thin :=  strtoint(FormSRsimu.EditN_thinning.text);
     nrow := NGibbs - NBurnIn;
     nGibbsRows := trunc(nrow/n_thin);

     setlength(para, np);
     setlength(theta_matrix, nGibbsRows, np);
     setlength(a_last_matrix, nGibbsRows, 1);
     setlength(Recruit_Pred, nGibbsRows, 1);

     MRandSeed(12345);
     //************************************************
     // initialize parameters
     //************************************************

     para[0] := ainit;// a = ln(alpha)
     para[1] := binit;//bhat;

     para[2] := 0.1; //sig2E;
     para[3] := 0.1; //sig2C;
     para[4] := 0.1; //tau2a;
     para[5] := 0.1; //tau2R;
     para[6] := 0.1; //tau2Lamb;

     // Conditional on S1 in year 1
     Gibbs_sampling(NGibbs, NBurnIn, n_thin, return_age,
         para, E, Cobs, R, theta_matrix,
           a_last_matrix, Recruit_Pred, ifail);

     Summary(theta_matrix, ElevenNum);
     Summary(a_last_matrix, ElevenNum_a);
     Summary(Recruit_Pred, ElevenNum_Rpred); // Recruitment prediction

     // Calculate hMSY and SMSY from posterior draws of parameters
     //  Using posterior mean estimates
     //  Need adjust lognormal error
     if frm_CLIM2.CheckBoxCalSmsybyMCMC.checked then
         SaveSummSimuSMSY(header, theta_matrix, a_last_matrix,
               hMSY, SMSY, PropAneg, PropBneg)
     else
     begin
          // calculate SMSY using mean posterior estimates
        aestimate := (ElevenNum_a[0, 0]);
        bestimate :=  ElevenNum[6,1];
        SMSY_hMSY('LambertW', Exp(aestimate), bestimate, hMSY, SMSY);
     end;
     setlength(NamesPara, np);
     NamesPara[0] := 'Ricker a';
     NamesPara[1] := 'beta';
     NamesPara[2] := 'sig2E';
     NamesPara[3] := 'sig2C';
     NamesPara[4] := 'tau2a';
     NamesPara[5] := 'tau2R';
     NamesPara[6] := 'tau2Lamb';

     if FormSRsimu.RadioButtonLambdaPriorHM.Checked then
     begin
       NamesPara[7] := 'mu_Lamb';
     end;

     // Traceplotting
     if FormSRsimu.CheckBoxPlotHistory.Checked then
     begin
        TracePlot('est', FormSRsimu, NamesPara, theta_matrix);
     end;

     // Time-varying Ricker 'a'
     //   Ricker 'a' of the last year a[t]
     //   Posterior mean of a_t[T]
     ParaEst[0] := (ElevenNum_a[0, 0]);

     ParaEst[1] := ElevenNum[6,1]; //beta
     ParaEst[2] := ElevenNum[6,2]; //sig2E
     ParaEst[3] := ElevenNum[6,5]; //tau2R
     ParaEst[4] := hMSY;
     ParaEst[5] := SMSY;

     // ****************************************
     // Special notes: Zhenming Su March 18, 2019
     // ****************************************
     // Use the posterior mean of the MCMC draws
     //  of R_pred as an estimate of predicted
     //  recruitment for year T+1.
     //  Adjustment for lognormal error is done
     //  automatically
     ParaEst[6] := ElevenNum_Rpred[0,0]; //mean

     // Show parameter values in the log window
     setlength(NamesPara, 12);
     NamesPara[0] := 'Ricka';
     NamesPara[1] := 'beta ';
     NamesPara[2] := 'sig2E';
     NamesPara[3] := 'tau2R';
     NamesPara[4] := 'hMSY ';
     NamesPara[5] := 'SMSY ';
     NamesPara[6] := 'RPred';
     NamesPara[7] := 'sig2C';
     NamesPara[8] := 'tau2a';
     NamesPara[9] := 'tau2L';
     NamesPara[10] := '%neg_a  %neg_b';

     if simu=1 then
     begin
       res_str := '';
       for i := 0 to 10 do
         res_str := res_str + NamesPara[i] + #9;
       MemoWrite(res_str);
     end;
     res_str := '';
     for i := 0 to 6 do
       res_str := res_str + Format('%3.2f',[ParaEst[i]])+#9;

     res_str := res_str +
                Format('%3.2f',[ElevenNum[6,3]])+#9+
                Format('%3.2f',[ElevenNum[6,4]])+#9+
                Format('%3.2f',[ElevenNum[6,6]]);
     res_str := res_str +#9+ Format('%5.2f',[PropAneg]) +'% '+
                Format('%5.2f',[PropBneg])+'%';
     MemoWrite(res_str);
end;


end.
