//*********************************************************************************************************************************
//  Code used in the MCMC sampling algorithm for a state-space model proposed by 
//    Su (2023) “Management performance evaluation of state-space models for 
//                 Pacific pink salmon stock-recruitment analysis” 
//                 submitted for peer review and publication
//  
//  (c) Zhenming Su, 2023
//
// This software is based on two previous articles listed below.
// Peterman, R.M., Pyper, B.J., and Grout, J.A. 2000. Comparison of parameter estimation methods for detecting 
//     climate-induced changes in productivity of Pacific salmon (Oncorhynchus spp.). Can. J. Fish. Aquat. Sci. 57: 181–191.
// Dorner, B., Peterman, R.M., and Su, Z. 2009. Evaluation of performance of alternative management models of 
//     Pacific salmon (Oncorhynchus spp.) in the presence of climatic change and outcome uncertainty 
//     using Monte Carlo simulations. Can. J. Fish. Aquat. Sci. 66: 2199–2221. 
//*******************************************************************************************************************************
unit CLIM2;

interface

  private
    { Private declarations }
  public
    { Public declarations }

  End;

var
  frm_CLIM2: Tfrm_CLIM2;

Type
  TManage = (OptEscape = 0, FixHarvest = 1);
  TEstation = (Std_Ricker, RICKAR1, KF, BSS, EXSSM, EIV, KFMCMC, TSRMCMC, EVPI);
  Tachscenario = (consta, step_func, step_ar1, decline, increase);
Var

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //
  // Simulation control parameters
  //
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  nTrials,               // number of simulation trials
  NYrs_Simu: Integer;    // number of years simulated per trial
  YrsEst: Integer;       // make estimates every YrsEst
  iEst: TEstation;       // Estimation Procedures
  Manag_Policy: TManage; // Type of management policy:
      // OptEscape-optimal escapement determined by Ricker parameter estimates

  CLIM_RndSeed: Integer;

  // Set width of data window in number of yrs of data to include in parameter estimation
  //   when nEstWindow > last year of simulation (NYrs_Simu) ALL data included in param.estim.
  //   e.g. nEstWindow = 15 includes only the most recent 15 years of data
  nEstWindow: Integer; // old: nWindow
  nMinEstYr: Integer;

  nPre_Manage: Integer; 

  HR_Pre_Manag: double; // = ExpRMin. A fixed harvest rate applied to
  // the first "nPre_Manage" pre-management years
  // of the simulation

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //
  // Parameter Estimates
  //
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  aEst, bEst: TVector;    // array[0..MAXYEARS] of double;
  VRecruitEst: TVector; // array[0..MAXYEARS] of double;  process error var
  hMSYEst, SMSYEst: double;
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //
  // Variables estimated or pre-set in the Manage Procedure
  //
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Catch,     // actual catch
  HR,        // harvest rate
  OptEscapeEst,      // Estimated optimal escapement
  OptE_True, // True optimal escapement calculated with true para. values
  OptEscapeAdj  // optimal escapement adjusted for safety margin and outcome uncertainty
     : TVector;

  Escape_Threshold: double; // = "MinEsc"
  // A conservation threshold used for "Low_Escp_Count"
  // to determine the number of times that the realized
  // escapement drops below the level specified by
  // "Escape_Threshold"

  HR_FixHarvPol: double;
  // Old expr. The proportional harvest rate used in the FIXED HARVEST RATE STRATEGY

  // deleted the variable "MinHarvestRate" that is not used anywhere
  MinWildestHR: double; // A fixed harvest rate applied to
  // the years of the simulation when parameter estimates go wild

  // Implementation error
  ImplErr: Boolean;
  OU: Integer; // Outcome uncertainty type 0,1,2
    // for OU type 1 and 2. Predicted recruitment from an estiamtion model
  SDouErr,
  RecruitPred: double; // used in OU1 and OU2

  // Safety margin: correction applied to optimal policy if RiskAverse set
  RiskAverse: Boolean = false;
  RiskCorrection: double;   // the proportion used to adjust S*

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //
  // CONSERVATION MEASURES
  //
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Accept_aMin, Accept_aMax, Accept_bMin, Accept_bMax: double;

  Low_Escp_Count, // =Lowcount: Count # cases in active management yrs
                  //    when realized escapement (S) is below the set level "Escape_Threshold"

  Para_Fail_count, // tallies no. of cases where the "MinWildestHR"
  // was applied for extreme parameter ests.

  HiLo_a_est, // tallies no. of cases where Accept_aMin<=aEst<=Accept_aMax
  HiLo_b_est, // tallies no. of cases where Accept_bMin<=bEst<=Accept_bMax
  HiLo_VR_est, // process variance < Accept_vRprocMax
  ExtinctionCount, // tallies no. of cases where ann. spawners <= Extinction threshold
  EscTargetFailcount // tallies no. of cases where ann. return <= esc. target
    : double;

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Operating model :
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Simulated Data
  Recruit, // = rec
  RecruitObs, Spawners, SpawnerObs, LnRS, LnRS_Obs: TV_NegLB; // =spa

  // LNRS, LnRS_Obs :TVector;  // = ObsLnRS

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //
  // Operating model : True populations quantities
  // and assumptions
  //
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ProcErrSDDefault: double; // SD of the marginal distribition of a recruitment process with AR1
  SDRecruitProcess: double; // the random error of recruitment process with AR1
  ar1_coeff: double; // recruitment process ar1 coefficient
  proc_err0, proc_err1 : double; // process err in recruitment in years t-1 and t

  SDObsErrSpawn, SDObsErrCatch: double; // measurement or observation error

  equilabund: double;
  // equilabund= the equilibrium number of spawners in an unfished pop
  // this variable is only used to initialize the true population

  SInitMin, SInitMax: double; // Min and max. values in range of Spawner values
  // randomly generated to seed the first "lag" yrs of the simulation

  PopAgeOption: Boolean;
  // Option for determining population structure of the fish stock:
  // false = even aged at maturity
  // True = age structured with age 4, 5, and 6 fish
  Return: TVector; // TMatrix;

  // Fixed age at return: pink
  ReturnAge: Integer; // Maximum age at which fish return for an even aged stock

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // True Ricker "a" and "b"
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  amin, amax, bmin, bmax: double; // max. and min. values allowed

  // for True "a" and "b" parameters
  aTrueInitDefault: double;
  aTrueInit: double; //True initial values of "a" parameter
  aTrue: TV_NegLB;

  bTrue: double;
  bTrueInit: double; // True initial values of "b" parameter

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // True "a" Change Scenarios
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  achscenario: Tachscenario;

  // Step function, achscenario = 2
  aIntial_StepDown: Boolean;
  aCh_Down: Boolean;
  ach_step, aCh_min: double;
  aCh_period: Integer;

  // AR(1) process, achscenario = 3
  phi_ach3: double;

  // old: autocorr. Acf for "a" change process if achscenario = 3
  sd_Ricker_a0: double;
  // IMPORTANT: The standard deviation of the variability in the "aTrue" process
  // z_ach3 :array of TV_NegLB; //old: autovariate

  // used in estimation Procedures: LinReg, BSS
  ndata : Integer;
  LnRSEst, CatchEst, EscapeEst, Robs, RecruitEst, CatchObs: TVector;

  XPlot, YPlot: TVector;
  AppPath: string;

var
  IniFileName: string;
  IniFile: TIniFile;

Procedure LoadParFile(FileName: string);
Procedure True_Ricker_a(itrial :integer; NYrs_Simu: Integer; var a_true_mat :TMatrix);
Procedure SystemInit;

implementation

uses
  RickerAR1, Kalman_Filter_MLE, CLIM2_Options, StateSpace_SR, Ev, EKF_SR,
  State_Space_Simu_SR, StateSpaceExpanded_SR, StateSpaceExpanded_SR_AR1, KF_SR,
  StdSR, OutPutUnit;

{$R *.dfm}

{ ------------------------------------------------------
  Procedure: Gen_aTrue
  Author:    zhenming
  Date:      12-Feb-2003
  Revised:   18-Oct-2018
  Revised:   6-Aug-2022
  Arguments: NYrs_Simu :integer
  Result:    aTrue
  ------------------------------------------------------}
Procedure True_Ricker_a(itrial :integer; NYrs_Simu: Integer; var a_true_mat :TMatrix);
var
  iyr, loops: Integer;
  aValue, aConsta, al0, alH, alL, z_ar1, z_ar0, sd_ar1, b_slope : double;
  drift, a_ub, b_ub, a_lb, b_lb, at_ub, at_lb : double;
//  a_spec :TVector;
Begin
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Scenarios for the changes in 'a'
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // achscenario=1 (constant)
  // achscenario=2 (step function)
  // achscenario=3 (step_ar1)
  // achscenario=4 (decline)
  //             5 (increase)
  // achscenario =(consta, step_func, step+ar1, decline, increase);

  // Zero out aTRUE
  for iyr := -Lag + 1 to NYrs_Simu do
  Begin
     aTrue[iyr] := 0.0;
  End;

  // Initialize the a_True scenarios
  aValue := 0;
  aConsta := 0;
  If (achscenario = consta) then
  Begin
    aConsta := RandG(aTrueInit, sd_Ricker_a0);
  End;

  // AR(1) process sd for STEP-AR1
  sd_ar1 := 0.0;
  if ((achscenario = step_func) or (achscenario = step_ar1)) then
  Begin
    // aIntial_StepDown = true:
    // aTrue = maximum value for 1 + n*Period <= iyr <= (Period/2) + n*Period
    // = minimum value for (Period/2 + 1) + n*Period <= iyr <= (Period) + n*Period
    if aIntial_StepDown then
      aCh_Down := True
    else
      aCh_Down := false;

    if aCh_Down then
    Begin   // aTrue changs to min value
      aValue := RandG(aTrueInit-ach_step, sd_Ricker_a0); 
    End
    else
    Begin   // aTrue changs to max value
      aValue := RandG(aTrueInit+ach_step, sd_Ricker_a0); 
    End;

    // calculate noise sd of the ar1 process: sd_ar1
    // phi_ach3 -- AR(1) process coefficient
    if (achscenario = step_ar1) then
        sd_ar1 := sd_Ricker_a0 * Sqrt(1 - phi_ach3 * phi_ach3);

  End;

  if (achscenario = decline) then
  begin
    // decline
    drift := -(2.5-1)/NYrs_Simu;
    a_ub := 3;
    b_ub := -(3-1.5)/NYrs_Simu;
    a_lb := 1.5;
    b_lb := -(2-0.5)/NYrs_Simu;
  end;
  if (achscenario = increase) then
  begin
    // increasing 
    drift := (2.5-1)/NYrs_Simu;
    a_ub := 1.5;
    b_ub := (3-1.5)/NYrs_Simu;
    a_lb := 0.5;
    b_lb := (2-0.5)/NYrs_Simu;
  End;

  // (1) Loops for generating the negative array elements of the
  //     time series of aTrue values for 'lag' years leading up to iyr=1

  // Initialize AR1
  z_ar0 := 0;
  If (achscenario = step_ar1) THEN
  Begin
    for loops := 1 to 100 do
    Begin
      z_ar1 := phi_ach3 * z_ar0 + RandG(0, sd_ar1); 
      z_ar0 := z_ar1;
    End;
  End;

  For loops := -Lag + 1 to 1 do
  Begin
    // Constant 'a'
    If achscenario = consta THEN
    Begin
      aTrue[loops] := aConsta;
    End;

    // Step function
    If achscenario = step_func THEN
    Begin
      aTrue[loops] := aValue;
    End;

    // Step-ar1 function
    If (achscenario = step_ar1) THEN
    Begin
      z_ar1 := phi_ach3 * z_ar0 + RandG(0, sd_ar1); 
      aTrue[loops] := aValue + z_ar1;
      z_ar0 := z_ar1;
    End;

    If achscenario = decline THEN
    Begin
      aTrue[loops] := Random * (3 - 2) + 2;
    end;
    if (achscenario = increase) then
    begin
      aTrue[loops] := Random * (1.5 - 0.5) + 0.5;
    End;

    // Keep aTrue value between amax and amin; same with bTrue
    If aTrue[loops] > amax THEN
    Begin
      aTrue[loops] := amax;
    End;

    If aTrue[loops] < amin THEN
    Begin
      aTrue[loops] := amin;
    End;
  End;

  // (2) loop over remaining years
  for iyr := 2 to NYrs_Simu do
  Begin
    If achscenario = consta THEN
    Begin
      aTrue[iyr] := aConsta;
    End;

    If achscenario = step_func THEN
    Begin
      if (((iyr - 1) mod round(aCh_period/2) = 0) and ((iyr - 1) > 2)) then
      Begin
        aCh_Down := not aCh_Down;

        if aCh_Down then
        Begin // aTrue changs to min value
          aValue := RandG(aTrueInit-ach_step, sd_Ricker_a0);
        End
        else
        Begin // aTrue changs to max value
          aValue := RandG(aTrueInit+ach_step, sd_Ricker_a0);
        End;
      End;
      aTrue[iyr] := aValue;
    End; // If;

    // Step + ar1 function
    If (achscenario = step_ar1) THEN
    Begin
      if (((iyr - 1) mod round(aCh_period/2) = 0) and ((iyr - 1) > 2)) then
      Begin
        aCh_Down := not aCh_Down;

        if aCh_Down then
        Begin // aTrue changs to min value
          aValue := RandG(aTrueInit-ach_step, sd_Ricker_a0);
        End
        else
        Begin // aTrue changs to max value
          aValue := RandG(aTrueInit+ach_step, sd_Ricker_a0);
        End;
      End;

      z_ar1 := phi_ach3 * z_ar0 + RandG(0, sd_ar1); 

      // aTrue value: step + ar1
      aTrue[iyr] := aValue + z_ar1;
      z_ar0 := z_ar1;
    End;

    If (achscenario = decline) or (achscenario = increase) then
    Begin
      aTrue[iyr] := aTrue[iyr-1] + drift + RandG(0, 0.1);
      at_ub := a_ub + b_ub * iyr;
      at_lb := a_lb + b_lb * iyr;
      if (aTrue[iyr] > at_ub) then aTrue[iyr] := at_ub;
      if (aTrue[iyr] < at_lb) then aTrue[iyr] := at_lb;
      if (aTrue[iyr] < 0.5) then aTrue[iyr] := 0.5;
      if (aTrue[iyr] > 3) then aTrue[iyr] := 3;
    End;

    // Keep aTrue value between amax and amin; same with bTrue
    If aTrue[iyr] > amax THEN
    Begin
      aTrue[iyr] := amax;
    End;

    If aTrue[iyr] < amin THEN
    Begin
      aTrue[iyr] := amin;
    End;
  End; // iyr

End; // Procedure aChange;


Procedure True_Pop_Even_aged_Init(itrial: Integer);
// This routine is the 'operating' model, calculating the dynamics of the
//   true population, whose underlying model will be estimated by the parameter
//   estimation routines
var
  t : Integer;
  logSR, prod_total, ObsErrEscape: double;
  equi_abund :double;
Begin
  // Set up and zero out vectors for each trial
  // Recruitment Returns for each return year
  //   for age-structured pop, this is the sum of recruits from different broodyears
  //   Return[1..iyr] compared to Spawners[-(k-1)..iyr] based on brood year
  setlength(Return, NYrs_Simu + 1 + ReturnAge);
  ZeroArray(Return);

  // The following variables are based on brood years
  //  ranged in (-lag + 1 to iyr)
  for t := -Lag + 1 to NYrs_Simu do
  Begin
    Spawners[t] := 0.0;
    SpawnerObs[t] := 0.0;
    Recruit[t] := 0.0;
    RecruitObs[t] := 0.0;
    LnRS[t] := 0.0;
    LnRS_Obs[t] := 0.0;
  End;

  // Initialize the true population
  //   for the first 'lag' brood years of
  //   spawners, recruits, and LNRS
  Lag := ReturnAge;
  equi_abund := equilabund;
  if (achscenario = decline) then equi_abund := 2.5;
  if (achscenario = increase) then equi_abund := 1;

  proc_err0 := 0;    // zero for each trial
  For t := -Lag + 1 to 0 do
  Begin
    // Spawners[-Lag + 1 .. iyr]:
    //   the actual spawners in millions, is calculated by
    //   multiplying equilabund by sinit (measured as a proportion
    //   of equilabund with values (0-1);
    Spawners[t] := equi_abund *
                 (Random * (SInitMax - SInitMin) + SInitMin);
    proc_err1 := proc_err0 * ar1_coeff + RandG(0, SDRecruitProcess);

    prod_total := aTrue[t] + proc_err1;
    if prod_total >= 4 then
       prod_total := 4;
       
    logSR := (prod_total - bTrue * Spawners[t]);

    //logSR := (aTrue[t] - bTrue * Spawners[t] + proc_err1);
    Recruit[t] := Spawners[t] * Exp(logSR);
    proc_err0 := proc_err1;

    ObsErrEscape := RandG(0, SDObsErrSpawn);
    // Random_Normal(0, sqrt(ObsErrSpawnV), MRNGRandom);
    SpawnerObs[t] := Spawners[t] * Exp(ObsErrEscape);

    if t + ReturnAge >= 1 then
    Begin
      Return[t + ReturnAge] := Recruit[t];
      LnRS[t] := Ln(Recruit[t] / Spawners[t]);

      // This LnRS_Obs is what is used by the parameter estimation routines
      // LnRS_Obs[lagyear] := Ln(RecruitObs[lagyear] / SpawnerObs[lagyear]);
    End;
  End; // Next lagyear;
End; // Procedure TruePop;

Procedure True_Pop_Even_Aged(itrial, t :Integer; NewEscapement, errOU, Return_Pred :double);
var
  logSR, ErrSobs, prod_total: double;
  Str1: String;
Begin
  // t = management year
  // NewEscapement = newly set actual escapement for year t
  // Return_Pred = predicted returns for year t; to be compared with true Return[t]

  // Update spawners, SpawnerObs, RecruitObs, LnRS_Obs for each brood year
  Spawners[t] := NewEscapement;  // Newly set spawner level
  ErrSobs := RandG(0, SDObsErrSpawn);  
  SpawnerObs[t] := Spawners[t] * Exp(ErrSobs);
  RecruitObs[t - ReturnAge] := SpawnerObs[t] + CatchObs[t];

  // This ObsLnRS is what is used by the TSR and KF estimation routines
  //   along with SpawnerObs
  LnRS_Obs[t - ReturnAge] :=
    Ln(RecruitObs[t - ReturnAge] / SpawnerObs[t - ReturnAge]);

  // Ricker stock-recruitment curve with AR1 error
  // calculate the number of recruits for this spawning escapement (of brood year iyr)
  // AR1 process error for next R based on Spawners[t] with Ricker model
  proc_err1 := proc_err0 * ar1_coeff + RandG(0, SDRecruitProcess); 
  prod_total := aTrue[t] + proc_err1;
  if prod_total >= 4 then
     prod_total := 4;
  logSR := (prod_total - bTrue * Spawners[t]);
  Recruit[t] := Spawners[t] * Exp(logSR);
  
  // True LogRS[t+ReturnAge]: only used on detailed printouts for checking
    // cases where there is more than one age at maturity and you use the first
    // one or two for estimating the 'a' parameter
  LnRS[t] := Ln(Recruit[t]/Spawners[t]);

  proc_err0 := proc_err1;

  if t + ReturnAge <= NYrs_Simu then
  Begin
    // Returns from brood year t
    Return[t + ReturnAge] := Recruit[t];
  End;

End;


{ -----------------------------------------------------------------------------
  Procedure: Management: management based on optimal escapement policy
  Author:    zhenming
  Date:      1-sep-2022
  Arguments: iyr :integer
  Result:    Escapement - new escapement for year t
             ErrOU = OU error multiplier
  ----------------------------------------------------------------------------- }
Procedure Management(t, itrial :Integer; Var Escapement, errOU: double);
var
  aTruePrime :double;
  aestimate, bestimate, aprime :double;
  VTrueRecruit, VRecruit :double;
  hMSY :double;
  // Outcome Uncertainty
  slopeImplErr, IntercImplErr,
  f_ou, g_ou, h_ou, RActual, errNormalOU: double;
  ErrCobs :double;
  onefish :double;
  a_estOK, b_estOK, VarR_estOK :boolean;
Begin
  // *********************************************************************************
  // COMPUTE 'true optimal escapement' for comparison in all routines.
  // *********************************************************************************

  VTrueRecruit := SDRecruitProcess * SDRecruitProcess

  If iEst = EVPI then
  Begin
    // for EVPI: estimates equal true parameter values;
    aestimate := aTrue[t];
    bestimate := bTrue;
    VRecruit  := VTrueRecruit;
  End
  Else
  Begin
    // for all other param. est. routines:
    //   estimates are most up-to-date estimates
    //    determined from previous years' data (no estimates for 'fixed harvest rate' scenario)
    // Su notes: aEst, bEst are derived from all previous years' data
    // aEst = (x[t],y[t], t = 1,2,3,...,t-1)
    aestimate := aEst[t];
    bestimate := bEst[t];

    // Adjust SMSY using process error variance
    if frm_CLIM2.CheckBoxAdjSMSY.Checked then
       VRecruit := VRecruitEst[t]
    else
       VRecruit := 0;
  End;

  // *********************************************************************************
  // Compute True optimal escapement
  //     using mean stock-recruit curve from lognormal error (sigma-adjusted)
  //     adjusted Ricker a value
  // *********************************************************************************
  aTrueprime := aTrue[t] + VTrueRecruit/2;

  SMSY_hMSY('LambertW', Exp(aTrueprime), bTrue, hMSY, OptE_True[t]);

  // *********************************************************************************
  // Generate a random deviate for implementation error (OU) EVEN WHEN
  // this error is NOT used. Puts all simulations on equal basis.
  // *********************************************************************************
  errNormalOU := 0;
  If ImplErr = FALSE Then
    errNormalOU := RandG(0, SDouErr);

  If ImplErr = True Then
  Begin
      if OU = 0 then // CLIM1: Randall et al. 2002
      begin
        errNormalOU := RandG(0, 0.2);
      End
      else if OU = 1 then // CLIM2 OU1: Dorner et al. 2009
      begin
        errNormalOU := RandG(0, 0.4);
      End
      else
      begin  // CLIM2 OU = 2, Dorner et al. 2009
        errNormalOU := RandG(0, 0.6);
      End
  End;

  // ---------------------------------------------
  // FOR OPTIMAL ESCAPEMENT STRATEGY
  // ---------------------------------------------
  // Estimated Optimal escapement:

  // Checking for parameter estimates out side the acceptiable ranges
  //  In these ranges, optimal escapement is not estimable
  a_estOK := true; b_estOK := true; VarR_estOK := true;
  If (aestimate <= Accept_aMin) or (aestimate >= Accept_aMax) then
  begin
      a_estOK := false;
      HiLo_a_est := HiLo_a_est + 1;
  End;
  if (bestimate <= Accept_bMin) or (bestimate >= Accept_bMax) then
  begin
      b_estOK := false;
      HiLo_b_est := HiLo_b_est + 1;
  End;
  if (VRecruit < 0) or (VRecruit >= 3) then
  begin
      VarR_estOK := false;
      HiLo_VR_est := HiLo_VR_est + 1;
  End;

  // When OptE cannot be calculated because of the above parameter failure, i.e.,
  //   parameter estimates have wild values, such as
  //   b_estimates are negative (density independence; exponential growth),
  // Using a small fixed HR value "MinWildestHR" to set the new escapement
  //   rather than applying an optimal escapement 
  If (not a_estOK) or (not b_estOK) or (not VarR_estOK) then
  Begin
    // Counts the number of years of all parameter failures;
    Para_Fail_count := Para_Fail_count + 1;

    // For wild parameter estimates,
    //   use a fallback fixed harvest rate policy to set the escapement;
    Escapement := Return[t] * (1 - MinWildestHR);
    OptEscapeEst[t] := Escapement;
  End;

  If (a_estOK) and (b_estOK) and (VarR_estOK) then
  begin
    // optimal escapement can be estimated from model estimates

      // for mle methods
      If (iEst = EVPI) or (iEst = Std_Ricker) or (iEst = RICKAR1) or (iEst = KF)
        or (iEst = EIV) or (iEst = KFMCMC) or (iEst = TSRMCMC) then
      Begin
        // Calcul. optimal escapement using lognormal-error-adjusted
        // (i.e. average or mean) stock-recruit curve;
        aprime := aestimate + VRecruit/2;
        hMSY := 0;
        SMSYEst := 0;
        SMSY_hMSY('LambertW', Exp(aprime), bestimate, hMSY, SMSYEst);
      End
      else
      begin
        // for MCMC sampling methods
        if not frm_CLIM2.CheckBoxCalSmsybyMCMC.Checked then
        begin
          // Adjust aprime and calculate SMSY using mean posterior estimates
          aprime := aestimate + VRecruit/2;
          SMSY_hMSY('LambertW', Exp(aprime), bestimate, hMSY, SMSYEst);
        End;
      End;

    OptEscapeEst[t] := SMSYEst;

    If Return[t] <= OptEscapeEst[t] Then
    Begin
      EscTargetFailcount := EscTargetFailcount + 1;
    End; 
  End;

  // Apply risk correction (safety margin) if selected;
  If RiskAverse Then
    // optimal escapement adjusted by a safety margin
    OptEscapeAdj[t] := OptEscapeEst[t] * (1 + RiskCorrection)
  else
    OptEscapeAdj[t] := OptEscapeEst[t];

  // -------------------------------------------------------
  // Applying OutCOME UNCERTAINTY/IMPLEMENTATION ERROR
  // -------------------------------------------------------
  // March 2018 Revised by Zhenming Su
  // Used to model situation where imperfect management and control over fisheries
  // when optimal escapement policy used.

  // Including CLIM1 implementation for Impl Error for Peterman et al. 2000
  // Doner et al. 2009 extended this to two types of outcome uncertainty (OU)
  //    onefish := 1.0e-05;
  errOU := 1.0;
  onefish := 1.0e-05;   
  If ImplErr = True Then
  Begin
    RActual := Return[t];
    if OU = 0 then // CLIM1
    begin
      // Peterman et al. 2000
      // Harvesting_CLIM1_ImplementationError (CLIM1_IMPL_ERR or LINEAR_IMPL_ERR)
      // Implementation error is a linear function of the number of returning recruits
      //  times a normally distributed random error term.
      // It is applied as a multiplicative error term to the target escapement:
      // actual_escapement = target_escapement * N(1.0, stdev) * (intercept + slope * recruits)
      IntercImplErr := 0.5;
      slopeImplErr  := 0.3;
      //sd = 0.2;
      errOU := (IntercImplErr + slopeImplErr * RActual) *
                      (1 + errNormalOU) // multiplicative error;
    End
    else if OU = 1 then  // CLIM2 OU
    begin
      f_ou := 0.5;
      h_ou := 0.4;
      //sd := 0.4;
      errOU := Exp(f_ou + h_ou * (Ln(RActual+onefish) - Ln(RecruitPred+onefish)) + errNormalOU);
    End
    else
    begin  //OU = 2
      f_ou := 0.1;
      g_ou := 0.3;
      h_ou := 0.4;
      //sd = 0.6;
      errOU := Exp(f_ou + g_ou * Ln(RActual+onefish) + h_ou *
                 (Ln(RActual+onefish) - Ln(RecruitPred+onefish)) + errNormalOU);
    End
  End;

  Escapement := OptEscapeAdj[t] * errOU;

  If Escapement <= onefish*2 Then    
  begin
      Escapement := 2*onefish;
      ExtinctionCount := ExtinctionCount+1;
  End;

  If Return[t] <= Escapement Then
  Begin
      Escapement := Return[t];
      //EscTargetFailcount := EscTargetFailcount + 1;  
  End;

  // --------------------------------------------------------------------------------------------
  // CONSERVATION MEASURE:
  //   Count number of cases in active management yrs
  //   when realized escapement (Escapement) is below set level
  If (t > nPre_Manage) and (Escapement < 0.175) Then
      Low_Escp_Count := Low_Escp_Count + 1;

  // -------------------------------------------------------
  // Compute Catch, Exploitation rate used
  // -------------------------------------------------------
  Catch[t] := Return[t] - Escapement;
  HR[t] := Catch[t]/Return[t];

  // Observed catch
  ErrCobs := RandG(0, SDObsErrCatch);
  CatchObs[t] := Catch[t] * Exp(ErrCobs);
End; // Procedure;

Procedure SetLength_Dyn_Arrays();
Begin
  If PopAgeOption = True THEN
  Begin
    setlength(EstRec, NYrs_Simu + 1);
    ZeroArray(EstRec);
  End;

  setlength(aEst, NYrs_Simu + 1);
  ZeroArray(aEst);
  setlength(bEst, NYrs_Simu + 1);
  ZeroArray(bEst);
  setlength(VRecruitEst, NYrs_Simu + 1);
  ZeroArray(VRecruitEst);

  setlength(Return, NYrs_Simu + 1 + ReturnAge);
  ZeroArray(Return);

  setlength(Catch, NYrs_Simu + 1);
  ZeroArray(Catch);
  setlength(HR, NYrs_Simu + 1);
  ZeroArray(HR);
  setlength(OptE_True, NYrs_Simu + 1);
  ZeroArray(OptE_True);
  setlength(OptEscapeEst, NYrs_Simu + 1);
  ZeroArray(OptEscapeEst);
  setlength(OptEscapeAdj, NYrs_Simu + 1);
  ZeroArray(OptEscapeAdj);
  setlength(CatchObs, NYrs_Simu + 1);
  ZeroArray(CatchObs);
End;

Procedure Performance_Measures(itrial, nPerf, yrPerf_Start: Integer;
  var res_vect: TVector; var F_Out_Summ_EachTrial: Textfile);
var
  Bias, vCatch, vSpawn: TVector;
  // ------------------------------------------------------------
  //
  // Bias measures :
  // mean "raw bias", "% bias", "absolute bias" and "sqr bias"
  // averaged over "nperf" years and ntrials
  //
  // ------------------------------------------------------------
  // ************** "a" *******************
  // for differences between aTrue(iyr) and aEst(iyr-1)
  // METHOD 1 -- Compare aTrue[ iyr] which affects broodyear [iyr]
  // with aEst[iyr - 1] which is used to set esc. for broodyear [iyr]
  // but is based on data up to brood year [iyr-EstAge-1].
  // Starting at yr 4 (when nPre_Manage < 4)
  // because least squares will not have any estimates until yr 3.
  avg_a_bias, // avg. bias (aestimated -aTrue), avged over nPerf years
  avg_a_pcBias, // avg. % bias (aestimated -aTrue)/aTrue, avged over nPerf years
  avg_a_absBias, // avg. ABS(Bias. aestimated -aTrue), avged over nPerf years
  ssq_a_Bias,
  // squared difference between aEst and aTrue, summed over nPerf years
  mse_a_Bias,
  // avg. squared difference between aEst and aTrue, avged over nPerf years

  // Ricker "b"
  avg_b_bias, // avg.bias (bestimated -bTrue), avged over nPerf years
  avg_b_pcBias, // avg. % bias (bestimated -bTrue)/bTrue, avged over nPerf years
  avg_b_absBias, // avg. ABS(Bias. bestimated -bTrue), avged over nPerf years
  ssq_b_Bias,
  // squared difference between bEst and bTrue, summed over nPerf years
  mse_b_Bias,
  // avg. squared difference between bEst and bTrue, avged over nPerf years

  // ************ "a2" ************************
  // METHOD 2 -- Compare aEst[iyr-1] with the 'aTrue' that
  // it is trying to estimate, i.e., aTrue[iyr-EstAge-1].
  // Not done here for 'b' because it remains constant.
  // Brood year (iyr) starting in yr. 3 because
  // least squares will not have any estimates until then.
  // for differences between aTrue(iyr-EstAge-1) and aEst(iyr-1)
  avg_a2_Bias, avg_a2_pcBIAS, avg_a2_absBias, ssq_a2_Bias, mse_a2_Bias,

  // "Escapement measure 1"
  // for differences between OptE and OptE_True
  avg_esc_Bias, // avg. bias: (OptE-OptE_True), avged over nPerf years
  avg_esc_pc_BIAS,
  // avg. % bias: (OptE-OptE_True)/OptE_True)*100, avged over nPerf years
  avg_abs_esc_Bias, // avg. ABS((OptE-OptE_True), avged oover nPerf years
  ssq_esc_Bias,
  // squared difference between OptE and OptE_True, summed over years
  mse_esc_Bias,
  // avg. squared difference between  OptE and OptE_True, avged over nPerf years

  // "Escapement measure 2"
  // for differences between actual esc. and OptE_True
  avg_esc2_bias, avg_esc2_pc_BIAS, avg_abs_esc2_Bias, ssq_esc2_Bias,
    mse2_esc_Bias: double;

  // --------------------------------------------------------
  // Total catch, log(catch), spawner abundance
  // sum of sqr(catch), sqr(spawners)
  // over "nperf" years for each trial
  // --------------------------------------------------------
  Cum_C, // cumulative catch over each nperf years' simulation;
  avg_C,
  Cum_lnc, // cumulative natural log catch over each nperf years' simulation;
  avg_spwn, // avrg cum. spawners over each nperf years' simulation;
  avg_HR // avg. harvest rate over years
    : double;

  // --------------------------------------------------------------
  // Variance and CV over years in catch and spawner abundance
  // --------------------------------------------------------------
  VAR_c, CV_c, VAR_spwn, CV_spwn, Q1_spwn: double;

  s_Summ_Stats: string;
  iyr, i: Integer;
Begin
  if nPerf >= 1 then
  Begin
    // ---------------------------------------
    // Bias statistics of paramameter estimates;
    // ---------------------------------------
    setlength(Bias, nPerf);

    // METHOD 1 -- Compare aTrue[iyr] which affects broodyear [iyr]
    // with aEst[iyr] which is used to set escapement for broodyear [iyr] but is
    // based on data up to brood year [iyr-EstAge].
    // Starting at yr 4 (when nPre_Manage < 4) because least squares
    // will not have any estimates until yr 3.
    // Do this for all param. estim. methods to put them all on equal footing.

    // Bias of Ricker "a[y]": yearly bias values of a_t
    // for differences between aTrue(iyr) and aEst(iyr)
    // a(Bias) = aEst[iyr] - aTrue[iyr];

    // aEst[i]: is the estimate based on data upto iyr-1
    for i := yrPerf_Start to NYrs_Simu do
    Begin
      Bias[i - yrPerf_Start] := (aEst[i] - aTrue[i]);
    End; // for
    avg_a_absBias := mean(absl(Bias));
    ssq_a_Bias := SumOfSquares(Bias);
    // mean squared error
    mse_a_Bias := ssq_a_Bias / nPerf;

    for i := yrPerf_Start to NYrs_Simu do
    Begin
      Bias[i - yrPerf_Start] := 100 * Bias[i - yrPerf_Start]/aTrue[i];
    End; // for
    avg_a_pcBias := mean(Bias);

    // bias of "b"
    // for differences between bTrue(iyr) and bEst(iyr-1)
    // vBias := pc_Bias1(bEst, bTrue, yrPerf_Start, iyr, 1);
    for i := yrPerf_Start to NYrs_Simu do
    Begin
      Bias[i - yrPerf_Start] := (bEst[i] - bTrue);
    End; // for
    avg_b_absBias := mean(absl(Bias));
    ssq_b_Bias := SumOfSquares(Bias);
    mse_b_Bias := ssq_b_Bias / nPerf;

    for i := yrPerf_Start to NYrs_Simu do
    Begin
      Bias[i - yrPerf_Start] := 100 * Bias[i - yrPerf_Start]/bTrue;
    End; // for
    avg_b_pcBias := mean(Bias);

    // METHOD 2 -- Compare aEst[iyr] with the 'aTrue' that it is trying
    // to estimate, i.e., aTrue[iyr-EstAge].
    // Not done here for 'b' because it remains constant.
    // Brood year (iyr) starting in yr. 3 because least squares
    // will not have any estimates until then.
    // Do this for all param. estim. methods to put them all on equal footing.
    // a2(Bias) =  aEst[iyr] - aTrue[iyr - EstAge];  (CLIM1)

    //  a2_Bias =  aEst[iyr] - aTrue[iyr - EstAge];
    for i := yrPerf_Start to NYrs_Simu do
    Begin
      Bias[i - yrPerf_Start] := (aEst[i] - aTrue[i-EstAge]);
    End; // for
    avg_a2_absBias := mean(absl(Bias));
    ssq_a2_Bias := SumOfSquares(Bias);
    mse_a2_Bias := ssq_a2_Bias / nPerf;

    for i := yrPerf_Start to NYrs_Simu do
    Begin
      Bias[i - yrPerf_Start] :=
                100 * Bias[i - yrPerf_Start]/aTrue[i-EstAge];
    End; // for
    avg_a2_pcBIAS := mean(Bias);

    // For differences between estimated OptE and OptE_True in iyr;
    // esc(Bias) =  OptE[iyr] - OptE_True[iyr];
    Bias := Bias_fun(OptEscapeEst, OptE_True, yrPerf_Start, NYrs_Simu);
    avg_abs_esc_Bias := mean(absl(Bias));
    ssq_esc_Bias := SumOfSquares(Bias);
    mse_esc_Bias := ssq_esc_Bias / nPerf;

    for i := yrPerf_Start to NYrs_Simu do
    Begin
      Bias[i - yrPerf_Start] := 100 * Bias[i - yrPerf_Start]/OptE_True[i];
    End; // for
    avg_esc_pc_BIAS := mean(Bias);

    // For differences between actual escapement (Spawners[iyr]) and OptE_True in iyr;
    // esc2(Bias) := Spawners[iyr] - OptE_True[iyr];
    Bias := Bias3(Spawners, OptE_True, yrPerf_Start, NYrs_Simu);

    avg_abs_esc2_Bias := mean(absl(Bias));
    ssq_esc2_Bias := SumOfSquares(Bias);
    mse2_esc_Bias := ssq_esc2_Bias / nPerf;

    for i := yrPerf_Start to NYrs_Simu do
    Begin
      Bias[i - yrPerf_Start] := 100 * Bias[i - yrPerf_Start] / OptE_True[i];
    End; // for
    avg_esc2_pc_BIAS := mean(Bias);

    // ---------------------------------------------
    // cumulative performance indices
    // for catch, harvest rate and #spawners
    // summed across nPerf years of simulations
    // ---------------------------------------------
    // catch
    setlength(vCatch, nPerf);
    // CopyVector(vCatch, catch, yrPerf_Start, NYrs_Simu);
    for i := yrPerf_Start to NYrs_Simu do
    Begin
      vCatch[i - yrPerf_Start] := Catch[i];
    End; // for

    Cum_C := Sum(vCatch);

    avg_C := mean(vCatch);
    VAR_c := variance(vCatch);

    // calc. interannual variability in catches for each trial using CV;
    If (VAR_c <= 0.0) Or (Cum_C <= 0) THEN
    Begin
      CV_c := 0.0;
    End
    Else
    Begin
      CV_c := Sqrt(VAR_c) / (Cum_C / nPerf);
    End; // If;

    Cum_lnc := Sum(Log(vCatch));

    // harvest
    avg_HR := SumRange(HR, yrPerf_Start, NYrs_Simu) / nPerf;

    // spawners
    setlength(vSpawn, nPerf);
    Copy(Spawners, vSpawn, yrPerf_Start, NYrs_Simu);
    avg_spwn := mean(vSpawn);
    VAR_spwn := variance(vSpawn);
    Q1_spwn := Quantile(length(vSpawn), 0.25, vSpawn);

    // variance and CV over time in spawner abundance;
    // calc. interannual variability in spawners for each trial using CV;
    If (VAR_spwn <= 0.0) Or (avg_spwn <= 0) THEN
      CV_spwn := 0.0
    Else
      CV_spwn := Sqrt(VAR_spwn) / (avg_spwn);

  
End;

{ -----------------------------------------------------------------------------
  Procedure: Clim2
  Author:    zhenming
  Date:      12-Feb-2022
  ----------------------------------------------------------------------------- }
Procedure CLIM();
// Date last revised:
// This version does the fully simulation and parameter estimation.
// Program to determine the optimal parameter estimation methods and
// management algorithms in the context of climate change
// and its effects on parameters of productivity of fish populations.
var
  t1, t2, t, i, j, k, iyr, itrial, npp,
    N_YrTrials, FirstCall : Integer;

  // Number of years of data for calculating performance measures following nPre_Manage.
  // Harvest rate policy for these years are based on updated parameter harvest rates
  YrEstStart, yrPerf_Start, nPerf: Integer;

  // ------------------------------------------------
  // Yearly quantities
  // ------------------------------------------------
  Escapement, err_OU: double;

  ParaEst: TVector;
  wt: TVector;
  spwnindex: TVector;
  ParameterEst: TVector;

  // Summary stats
  statStr: string;
  Str, est_method, s_Annual_Est, s_Summ_Stats, s_Simu_Files_Prefix: String;

  SumResults: TMatrix;
  avg: TVector;

  // Output file of detailed output from simulations
  Ann_Est_Out: Boolean;
  Fname, tempStr, Appl_Header: string;
  F: Textfile;
  F_Out_Annual_Est,     // "_Ann_Est" detailed
  F_Out_Summ_EachTrial, // "*_Sum_stat"
  F_OUT_EST: Textfile;

  a_true_tmp :double;
  a_str :string;
  a_spec :TVector;
  a_true_mat :TMatrix;
  Fdata :Textfile;
  
  res_str, end_time : string;
  NamesPara: array of string;
Begin

  try

    // NYrs_Simu and nPerf below must be AT LEAST 3.
    nPerf := NYrs_Simu - nPre_Manage;
    YrsEst := StrToInt(frm_CLIM2.EditUpdateEst.text);

    // The start year of calculating performance indices
    yrPerf_Start := nPre_Manage + 1;
    nMinEstYr := nPre_Manage;

    bTrue := bTrueInit;

    aTrueInit := aTrueInitDefault;
    equilabund := aTrueInit / bTrueInit;

    ar1_coeff := strtofloat(frm_CLIM2.Editar1_coeff.text);

    // process error sd
    SDRecruitProcess := ProcErrSDDefault * Sqrt(1 - ar1_coeff * ar1_coeff);

    try

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // Initialize the estimators
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      If (iEst = Std_Ricker) or (iEst = RICKAR1) or (iEst = KF) or (iEst = EIV) or
           (iEst = BSS) or (iEst = EXSSM) or (iEst = KFMCMC) or (iEst = TSRMCMC)
      Then
      Begin
        setlength(spwnindex, NYrs_Simu + 1);
           ZeroArray(spwnindex);
        setlength(wt, NYrs_Simu + 1);
        wt := vector(NYrs_Simu + 1, 1.0);
        setlength(ParameterEst, 7);
           ZeroArray(ParameterEst);
      End;

      case iEst of
        Std_Ricker:
          est_method := 'TSR';
        RICKAR1:
          est_method := 'RICKAR1';
        KF:
          est_method := 'KF';
        BSS:
          Begin
            est_method := 'BSS';
            FormSRsimu.RadioGroup_Custom_Model.ItemIndex := 2;
          End;
        EXSSM:
          Begin
            est_method := 'ExSSM';

            if FormSRsimu.RadioGroup_Custom_Model.ItemIndex = 4 then
            begin
              npp := 9; // 7
              setlength(NamesPara, npp);
              if FormSRsimu.RadioButtonLambdaPriorHM.Checked then
              Begin
                npp := 9;
                setlength(NamesPara, npp);
              End;

              NamesPara[0] := 'alpha';
              NamesPara[1] := 'beta';
              NamesPara[2] := 'sig2E';
              NamesPara[3] := 'tau2R';
              NamesPara[4] := 'hMSYEst';
              NamesPara[5] := 'SMSYEst';
              NamesPara[6] := 'sig2C';
              NamesPara[7] := 'tau2a';
              NamesPara[8] := 'tau2Lamb';

              { if FormSRsimu.RadioButtonLambdaPriorHM.Checked then
                Begin
                NamesPara[7] := 'mu_Lamb';
                End; }
            End
            else
            begin
              npp := 13;
              setlength(NamesPara, npp);
              NamesPara[0] := 'alpha';
              NamesPara[1] := 'beta';
              NamesPara[2] := 'sig2E';
              NamesPara[3] := 'tau2R';
              NamesPara[4] := 'hMSYEst';
              NamesPara[5] := 'SMSYEst';
              NamesPara[6] := 'sig2C';
              NamesPara[7] := 'tau2a';
              NamesPara[8] := 'tau2Lamb';
              // NamesPara[9] := 'x3';
              // NamesPara[10] := 'x4';
              // NamesPara[11] := 'x5';
              // NamesPara[12] := 'x6';
            End;

            res_str := '';
            for i := 0 to npp - 1 do
              res_str := res_str + NamesPara[i] + #9;

            { res_str := res_str + '--' +#9;
              for i := npp to npp+NYrs_Simu-1 do
              res_str := res_str  +'a[' + inttostr(i-npp+1) + ']'+ #9; }
            // MemoWrite(res_str);
          End;

        EIV:
          Begin
            est_method := 'EV';
            FormSRsimu.RadioGroup_Custom_Model.ItemIndex := 1;
          End;
        KFMCMC:
          Begin
            est_method := 'KFMCMC';
            FormSRsimu.RadioGroup_Custom_Model.ItemIndex := 6;
          End;
        TSRMCMC:
          Begin
            est_method := 'TSRMCMC';
            FormSRsimu.RadioGroup_Custom_Model.ItemIndex := 7;
          End;
        EVPI:
          est_method := 'EVPI';
      else
        est_method := 'NA';
      End;

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // Loop over Monte Carlo trials for an estimator
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      For itrial := 1 To nTrials do
      Begin

        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        // Cumu. quantities over years
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        // For each Monte Carlo trial reset all the initial conditions
        ExtinctionCount := 0;
        // Number of cases where OptE < Escape_Threshold;
        Low_Escp_Count := 0;
        // Number of cases where ann. return <= esc. target;
        EscTargetFailcount := 0;
        // Number of cases where the MinWildestHR was applied for extreme parameter ests.
        Para_Fail_count := 0;
        // Number of cases where 0<=baEst <=4;
        HiLo_a_est := 0;
        // Number of cases where 0<= bEst <=5;
        HiLo_b_est := 0;
        HiLo_VR_est := 0;

        SetLength_Dyn_Arrays();

        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // (1) Generate aTrue for for all years
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        RandSeed := CLIM_RndSeed + 10*itrial;  

        True_Ricker_a(itrial, NYrs_Simu, a_true_mat);

        RandSeed := CLIM_RndSeed + itrial; 
        setlength(SumResults, itrial, 20);
        ZeroArray(SumResults[itrial - 1]);

        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // (2) Initialize true population Spawner-Recruit data
        //      for the first lag years
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        True_Pop_Even_aged_Init(itrial);

        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // Loop over years
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        FirstCall := 0;
        For iyr := 1 To NYrs_Simu do
        Begin

          // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          // (1) make estimation
          // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          // Assessment model estimation
          If (iEst = Std_Ricker) or (iEst = RICKAR1) or (iEst = KF) or
               (iEst = BSS) or (iEst = EXSSM) or (iEst = KFMCMC) or
                 (iEst = TSRMCMC) or (iEst = EIV) then
          Begin
            // Set the number of data points for estimation routines;
            // Set up starting and finishing years for estimation window of relevant data
            If NYrs_Simu <= nEstWindow then
            begin
              //use ALL data
              YrEstStart := 1;
              ndata := iyr - 1;
            End
            Else
            begin
              // use only nEstWindow data points
              If iyr <= nEstWindow then
              begin
                YrEstStart := 1;
                ndata := iyr-1;
              End
              else
              begin
                YrEstStart := iyr - (nEstWindow) + 1;
                ndata := nEstWindow;
              End;
            End;

            // This section deals with age-structured fish populations
            //    for least squares method,
            // a spawner index was used to remove the estimation age (EstAge)
            //    lag which creates dimensioning problems.
            if (iEst = Std_Ricker) or (iEst = RICKAR1) or (iEst = KF) then
            Begin
              setlength(spwnindex, ndata);
              setlength(LnRSEst, ndata);
              for t := 0 to ndata - 1 do
              Begin
                spwnindex[t] := SpawnerObs[t+YrEstStart-EstAge]; // SpawnerObs[t - EstAge + 1];
                LnRSEst[t]   := LnRS_Obs[t+YrEstStart-EstAge];   // LnRS_Obs[t - EstAge + 1];
              End;
            End;

            if (iEst = BSS) or (iEst = EXSSM) or
                (iEst = KFMCMC) or (iEst = TSRMCMC) or (iEst = EIV)  then
            Begin
              setlength(EscapeEst, ndata + EstAge);
              setlength(RecruitEst, ndata + EstAge);
              setlength(Robs, ndata + EstAge);
              setlength(CatchEst, ndata + EstAge);
              if iyr = 1 then
              Begin
                for t := 0 to EstAge - 1 do
                Begin
                  EscapeEst[t] := SpawnerObs[t - EstAge + 1];
                  CatchEst[t] := 0;
                  RecruitEst[t] := 0;
                  Robs[t] := 0;
                End;
              End
              else
              Begin
                EscapeEst[iyr] := SpawnerObs[iyr - EstAge + 1];

                CatchEst[iyr] := CatchObs[iyr - 1];
                RecruitEst[iyr] := Return[iyr - 1];
                Robs[iyr] := (CatchObs[iyr - 1] + SpawnerObs[iyr - 1]);
              End;
            End;

            if (iEst = Std_Ricker) then
            Begin // Standard Ricker model
              If (iyr >= nPre_Manage - max(YrsEst,EstAge)) then
              Begin // only do estimation for years > nPre_Manage - YrsEst //EstAge+1
                if (iyr mod YrsEst) = 0 then
                Begin // estimating every YrsEst year
                  LeastSq(spwnindex, LnRSEst, wt, length(spwnindex), YrEstStart-1, ParameterEst);
                  aEst[iyr] := ParameterEst[0];
                  bEst[iyr] := ParameterEst[1];
                  VRecruitEst[iyr] := ParameterEst[2];
                End
                else
                Begin
                  aEst[iyr] := ParameterEst[0];
                  bEst[iyr] := ParameterEst[1];
                  VRecruitEst[iyr] := ParameterEst[2];
                End;
                If iyr > nPre_Manage then
                   RecruitPred := PredictTSR(SpawnerObs[iyr-EstAge], aEst[iyr-EstAge],
                                             bEst[iyr-EstAge], VRecruitEst[iyr-EstAge]);
              End
              else
              Begin
                // zero parameter estimates so that for the first EstAge = 2 yrs of
                // simulation they're output as 0 from least squares.
                aEst[iyr] := 0;
                bEst[iyr] := 0;
                VRecruitEst[iyr] := 0;
                RecruitPred := 0;
              End; // If;
            End; // std Ricker

            if (iEst = KF) then
            Begin
              If (iyr >= nPre_Manage - max(YrsEst,EstAge)) then
              Begin
                if (iyr mod YrsEst) = 0 then  // estimating every YrsEst year
                Begin
                  // KalmanFilter_Grid(nd, spwnindex, LnRSEst, BV);
                  Kalman_Filter(length(spwnindex), spwnindex, LnRSEst, ParameterEst);
                  aEst[iyr] := ParameterEst[1];
                  bEst[iyr] := ParameterEst[2];
                  VRecruitEst[iyr] := ParameterEst[3];
                End
                else
                Begin
                  aEst[iyr] := ParameterEst[1];
                  bEst[iyr] := ParameterEst[2];
                  VRecruitEst[iyr] := ParameterEst[3];
                End;
                If iyr > nPre_Manage then
                     RecruitPred := PredictKF(SpawnerObs[iyr-EstAge], aEst[iyr-EstAge],
                                          bEst[iyr-EstAge], VRecruitEst[iyr-EstAge]);

              End
              else
              Begin
                // zero parameter estimates so that for the first 2 yrs of
                // simulation they're output as 0 from least squares.
                aEst[iyr] := 0;
                bEst[iyr] := 0;
                VRecruitEst[iyr] := 0;
                RecruitPred := 0;
              End; // If;
            End; // If

            if (iEst = BSS) or (iEst = EXSSM) or (iEst = EIV) or (iEst = KFMCMC)
              or (iEst = TSRMCMC) then
            Begin
              if (iyr >= nPre_Manage - max(YrsEst,EstAge)) then // nPre_Manage   iyr>=EstAge+7
              Begin
                // E(t), C(t), R(t), Robs(t+k)
                if ((iyr) mod YrsEst) = 0 then  // estimating every YrsEst year
                Begin
                  inc(FirstCall);
                  if iEst = BSS then
                    Gibbs_StateSpaceSR_Est(ReturnAge, EscapeEst, CatchEst,
                      RecruitEst, Robs, aTrueInit, bTrueInit, SDObsErrSpawn*SDObsErrSpawn,
                      SDRecruitProcess * SDRecruitProcess, ParameterEst)
                  else if iEst = EXSSM then
                  begin
                    if FormSRsimu.RadioGroup_Custom_Model.ItemIndex = 3 then
                      Gibbs_ExpStateSpaceSR_AR1_Est(FirstCall, ReturnAge, EscapeEst,
                        CatchEst, RecruitEst, Robs, aTrueInit, bTrueInit,
                        SDObsErrSpawn*SDObsErrSpawn, SDRecruitProcess * SDRecruitProcess, ParameterEst)
                    else
                      Gibbs_ExpStateSpaceSR_Est(FirstCall, ReturnAge, EscapeEst, CatchEst,
                        RecruitEst, Robs, aTrueInit, bTrueInit, SDObsErrSpawn*SDObsErrSpawn,
                        SDRecruitProcess * SDRecruitProcess, ParameterEst)
                  End
                  else if iEst = EIV then
                    EVOptim_clim(ReturnAge, EscapeEst, CatchEst, RecruitEst,
                      Robs, aTrueInit, bTrueInit, SDObsErrSpawn*SDObsErrSpawn,
                      SDRecruitProcess * SDRecruitProcess, ParameterEst)
                  else if iEst = KFMCMC then
                    Gibbs_KF_SR_Est(ReturnAge, EscapeEst, CatchEst, RecruitEst,
                      Robs, aTrueInit, bTrueInit, SDObsErrSpawn*SDObsErrSpawn,
                      SDRecruitProcess * SDRecruitProcess, ParameterEst)
                  else
                    Gibbs_TraditionSR_est(ReturnAge, EscapeEst, CatchEst,
                      RecruitEst, Robs, aTrueInit, bTrueInit, SDObsErrSpawn*SDObsErrSpawn,
                      SDRecruitProcess * SDRecruitProcess, ParameterEst);

                  if (ParameterEst[0] = 1) then
                  Begin
                    ParameterEst[0] := (aTrueInit);
                    ParameterEst[1] := bTrueInit;
                    VRecruitEst[iyr] := SDRecruitProcess * SDRecruitProcess;
                  End;
                  aEst[iyr] := ParameterEst[0];
                  bEst[iyr] := ParameterEst[1];
                  VRecruitEst[iyr] := ParameterEst[3];
                  hMSYEst := ParameterEst[4];
                  SMSYEst := ParameterEst[5];
                  RecruitPred := ParameterEst[6];

                  {If (iyr > nPre_Manage) and (iEst = TSRMCMC)  then
                     RecruitPred := PredictTSR(EscapeEst[iyr+1-EstAge], aEst[iyr-EstAge],
                                               bEst[iyr-EstAge], VRecruitEst[iyr-EstAge]);
                  If (iyr > nPre_Manage) and (iEst = KFMCMC)  then
                     RecruitPred := PredictKF(EscapeEst[iyr+1-EstAge], aEst[iyr-EstAge],
                                               bEst[iyr-EstAge], VRecruitEst[iyr-EstAge]); }
                End
                else
                Begin
                  // Non estimation years
                  aEst[iyr] := ParameterEst[0];
                  bEst[iyr] := ParameterEst[1];
                  VRecruitEst[iyr] := ParameterEst[3];
                  hMSYEst := ParameterEst[4];
                  SMSYEst := ParameterEst[5];
                  RecruitPred := ParameterEst[6];
                  
                  {If (iyr > nPre_Manage) and (iEst = TSRMCMC)  then
                     RecruitPred := PredictTSR(EscapeEst[iyr+1-EstAge], aEst[iyr-EstAge],
                                            bEst[iyr-EstAge], VRecruitEst[iyr-EstAge]);
                  If (iyr > nPre_Manage) and (iEst = KFMCMC)  then
                     RecruitPred := PredictKF(EscapeEst[iyr+1-EstAge], aEst[iyr-EstAge],
                                            bEst[iyr-EstAge], VRecruitEst[iyr-EstAge]); }
                End;
              End
              else
              Begin
                // zero parameter estimates so that for the first k yrs of
                // simulation they're output as 0 from least squares.
                aEst[iyr] := 0;
                bEst[iyr] := 0;
                VRecruitEst[iyr] := 0;
                RecruitPred := 0;
              End; // If;
            End;
          End;

          // EVPI (expected value of perfect information);
          If iEst = EVPI THEN
          Begin
            aEst[iyr] := aTrue[iyr];
            bEst[iyr] := bTrue;

            // Still include ProcErrSD (variation around Ricker curve) in EVPI calculation,
            //    because you need to correct the calcul. of estimated optimal escapement
            //    in the Ricker curve for the +sigmaSquared/2 term (i.e. expected value =
            //                                                       mean + sigmSquared/2).
            VRecruitEst[iyr] := SDRecruitProcess * SDRecruitProcess;
            RecruitPred := Return[iyr];
            If (iyr > nPre_Manage) then
               RecruitPred := PredictTSR(Spawners[iyr-EstAge], aTrue[iyr-EstAge],
                                         bTrue, VRecruitEst[iyr-EstAge]);
            // Zero probability measures of wild parameter estimates,
            //   which are not possible for EVPI;
            HiLo_a_est := 0;
            HiLo_b_est := 0;
            HiLo_VR_est := 0;
          End; // If;
          // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          //
          // (4) Conduct harvesting and management
          //
          // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          If (iyr <= max(ReturnAge, nPre_Manage)) then
             PreManagement(iyr, itrial, Escapement)
          else
             Management(iyr, itrial, Escapement, err_OU);  

          // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          //
          // (5) Generate the true population Spawner-recruit data
          // for brood year "iyr" and also the corresponding
          // annual return and LnRS data for the return year
          // "iyr + ReturnAge"
          //
          // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          True_Pop_Even_Aged(itrial, iyr, Escapement, err_OU, RecruitPred);

        // --------------------
        End; // Next iyr;
        // --------------------

        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // Calculate performance indices for each trial
        // take into account only most recent "nPerf" yrs
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        Performance_Measures(itrial, nPerf, yrPerf_Start,
          SumResults[itrial - 1], F_Out_Summ_EachTrial);

      End; // Next itrial;

    finally
      Close(F_Out_Summ_EachTrial);
    End;
  finally
  End;
End; // CLIM2;

End.
