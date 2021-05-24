unit HVproj;

interface


uses
  Dpconst, HVLink, HVdefs, HVData, GBDefs;

type
THVProj = Class (TObject)
private
  p                     : byte;
  (* variables used to control flow of calcs*)
  _ts,_ARTYrIdx         : integer;
  bBalanceSexActs,
  bHasNewInfections,
  bHasHIV,
  bDoStockFlow          : boolean;


  (*main population variables*)
  _Adults               : HV_AdultsArray;

  _AddlRiskGrpARTCov    : HV_SexByRiskByYear;
  _RiskGrpARTCov,
  _HV_RiskGrpARTPerc    : HV_AdultsArray;
  _DP_RiskGrpARTPerc    : HV_DPRiskGrpARTPerc;

  _AdultsTS             : HV_AdultsTS;//handle to the _Adult array structure prior to updates

  _TotalAdultsHIV,
  _TotalAdultsART,
  _Populations          : HV_SexByYear;

  (* stock flow *)
  _TotalPops,
  _TotalEntrants,
  _TotalExits           : HV_TotalPops;

  (* epidemiological/model variables *)
  _EpidemicStYr         : integer;
  _MonthsInPrimaryStage,
  _YearsInPrimaryStage  : double;
  _TransMultM,
  _TransHIVF,
  _TransMultSTI,
  _CondomEff,
  _InitialPulse,
  _ARTSurviving,
  _TransMultMSM         : double;
  _Infectiousness       : HV_InfectiousnessAr;
  _STIPrev              : HV_STIPrev;
  _MalesCircum          : HV_TYearMatrix;
  _RedWhenCircum        : HV_RedWhenCircum;

  (*behavioral*)
  _CondomPercent        : HV_CondomPercent;
  _SexActs              : HV_SexActs;
  _NumPart              : HV_NumPart;
  _Behavior             : HV_Behavior;
  _IncRecruitment       : HV_IncRecruitment;
  _BehaviorChangeRate,
  _RiskGroupPercent     : HV_RiskGroup;
  _PercMarried          : HV_PercMarried;
  _AgeFirstSex          : HV_SexByYear;
  _PerIDUSharing        : HV_TYearMatrix;

  (* impact matrix *)
  _AdjCoverage          : HV_RN_TInterventions;

  (* vaccines *)
  _TNType               : byte;
  _Targeting            : byte;
  _Unvaccinated         : HV_TYearMatrix;
  _NewVaccinations      : HV_SexByRiskByYear;
  _Veffect              : HV_RN_Veffect;
  _VaccineEffectiveness : HV_RN_VaccineEffectiveness;
  _VacCoverage          : HV_RN_VacCoverage;
  _Vaccinated           : HV_Vaccinated;
  _VaccineBehavEffect   : HV_RN_VaccineBehavEffect;
  _TotalNewVaccinations : HV_TYearMatrix;
  _TotalVaccinated      : HV_SexByYear;

  (* new PrEP method mix functionality *)
  _PrEPPercent          : HV_RN_PrEPCoverage;
  _PrEPEff              : HV_PrEPEff;

  (* new entrants and migration *)
  _Pop15rate            : HV_Pop15Rate;
  _MigrNum,
  _MigrRate             : HV_YearSexHIV;
  _Migrant              : HV_Migrant;
  _TotalNewAdults,
  _TotalMigrants        : HV_TYearMatrix;

  (* force of infection *)
  _ForceInf             : HV_SexByYear;
  _rMultAll,
  _CalcPrevalence       : HV_SexByRiskByYear;
  _MultART              : HV_MultART;
  _MultNoART            : HV_HIVStatus;
  _PulseByStage,
  _NewInfByStage        : HV_NewInfByStage;

  (* new infections *)
  _NewInfections        : HV_NewInfections;
  _NewInfectionsT       : HV_NewInfectionsT;
  _NewHIV               : array [DP_BothSexes .. DP_Female] of double;
  _Incidence,
  _IncSexRatio,
  _TotalNewInfection    : HV_TYearMatrix;

  (* mortality *)
  _NonAIDSDeathRate     : HV_SexByYear;
  _AIDSDeaths           : HV_AIDSDeaths;
  _AIDSDeathsART        : HV_SexByRiskByYear;
  _TotalAIDSDeaths,
  _TotalARTDeaths       : HV_TYearMatrix;

  (* exit rates *)
  _RateofAging          : HV_SexByYear;
  _RateofHIVAging       : HV_RateOfHIVAging;
  _RateofARTAging       : HV_RateofARTAging;
  _ExitRate             : HV_ExitRate;
  _ExitRateART          : HV_ExitRateART;

  (* CD4, ART related *)
  _CD4_lower_limit      : HV_CD4_lower_limit;
  _lambda, _mu,
  _ScaleMuART           : HV_HIVStatusBySex;
  _alpha                : HV_CD4ARTBySex;
  _ProgressRate         : HV_ProgressRate;
  _StageExits,
  _StageProgressors     : HV_Stage;

  (* ART *)
  _CD4ThreshHoldAdults  : double;
  _h_set,_h_cov         : integer;
  _ARTReceive           : HV_TYearMatrix;
  _NewlyOnART           : HV_NewlyOnART;
  _NewlyEligibleART     : HV_NewlyEligibleART;

  _CD4Coverage          : array[2..10,0..100] of double;

   HV                   : THVData;

//  OutFile             : TextFile;

  (* called from Create *)
  procedure Transfer_Local;
//  procedure Init_WithValuesFromFP;
  procedure Init_Veffect;
  procedure Init_RiskGroupPercent;
  procedure Init_DynamicYearArrays;
  procedure Init_Variables;
  procedure Init_Disinhibition;
  procedure Init_IntCD4Vars;
  procedure Init_GetFirstYearART;
  procedure Init_AdjCoverage;
  procedure Init_IDUBehavior;

//  (* called from PreOuterLoopCalcs *)
//  procedure Calc_NewEntrants(t:integer);
//  procedure Calc_InitAgingRates;

  (* called from PreInnerLoopCalcs *)
  procedure Calc_BehavioralCoverage(t:integer);
  procedure Calc_SetARTEligibility(t: integer);
  procedure Calc_CD4Trans(t: integer);
  procedure Calc_NonAIDSMortAIM(t : integer);
  procedure Calc_NewEntrantsAIM(t:integer);
  procedure Calc_MigrNum(t:integer);
  procedure Calc_SumPopGroups(t:integer);

  (* called from InnerLoopCalcs *)
  function  Calc_HasHIV(t:integer):boolean;
  procedure Calc_PopulationTotals(t:integer);
  procedure Calc_MigrRate(t:integer);
  procedure Calc_HIVNegAgingAIM(t,_ts : integer);
  procedure Calc_HIVPosAgingAIM(t,_ts: integer);
  procedure Calc_HIVExitRates(t:integer);
  procedure Calc_HIVExitsProgressors(t:integer);
  //procedure Calc_ARTCoverage(t:integer); falls away
  procedure Calc_VaccinationTotals(t:integer);
  procedure Calc_NewVaccinations(t:integer);
  procedure Calc_ProgressNoRiskHIVNeg(t,_ts:integer);
//  procedure Calc_CD4TransAIM(t,t1,_ts: integer);
  procedure Calc_ProgressNoRiskHIVPos(t:integer);
  procedure Calc_Prevalence(t,_ts:integer);
  function  Calc_NewInfections(t:byte):boolean;
  procedure Calc_DistributeNewRecruits(t:integer);
  procedure Calc_ProgressHIVneg(t,_ts:integer);
  procedure Calc_ProgressHIVPrimary(t:integer);
  procedure Calc_ProgressHIVpos(t,_ts:integer);
  procedure Calc_ProgressART(t,_ts:integer);
  procedure Calc_AllocateART(t: integer);
  procedure Calc_AddNewInfections(t,TS:integer);
  procedure Calc_AddNewInfectionsAIM(t:integer);
  procedure Calc_InitialPulse(t,_ts:integer);
  procedure Calc_InitialPulseIDU(t,_ts:integer);
  procedure Calc_TrackAIDSDeaths(t,_ts:integer);
  procedure Calc_TrackNewInfections(t:integer);
  procedure Calc_TotalExits(t,_ts:integer);
  procedure Calc_TotalEntrants(t,_ts:integer);

  (* called from PostInnerLoopCalcs *)
  procedure Calc_IncSexRatio(t: integer);
  procedure Calc_TotalPops(t,t1:integer);
  procedure Calc_HIVMortTotals(t:integer);
  procedure Calc_HIVPopTotals(t:integer);

  (* called from PostOuterLoopCalcs *)
  procedure Calc_Incidence;
  procedure Calc_PercentPop;
  procedure Calc_ARTCoverageDisplay;
  procedure Calc_WriteStockFlow;

  (* called from Destroy *)
  procedure ReverseTransfer_Local;

   (*Additional ART coverage for specific risk groups: Modeling Consortium project   *)
  procedure Calc_AdditionalARTCov(t:integer);
//  function HV_GetBirthInterval(p : byte) : double;

public

  constructor Create(Proj : byte);
  //procedure PreOuterLoopCalcs; not need anymore
  procedure PreInnerLoopCalcs(t : integer);
  procedure InnerLoopCalcs(t,t1 : integer);
  procedure PostInnerLoopCalcs(t: integer);
  procedure PostOuterLoopCalcs;

  (*Additional ART coverage for specific risk groups: Modeling Consortium project   *)
  function GetRiskGrpARTCov(s,r,h,v,t:byte): double;
  function GetDP_RiskGrpARTPerc(s,h,t:byte): double;

  //number on people treatment by sex,
  procedure CalcNumOnART(t, t1 : integer; var numOnART : GB_TMaleFemale);

  procedure SetCD4Coverage(t, c: byte; value : double);
  destructor Destroy; override;
end;

implementation

uses
  SysUtils, Math, GBConst, GBFlexCel, GBData, GBMDUtil, DPData, HVConst,
  HVUtil, DPUtil;

Constructor THVProj.Create(Proj : byte);
begin
  p := Proj;
  HV := GetHVData(proj);

  Init_DynamicYearArrays;

  //IDU behav variables must me passed to module variables with Transfer_Local
  Init_IDUBehavior;

  Transfer_Local;

  (* toggle 'stockflow' monitoring functionality: Calc_TotalPops, Calc_TotalExits, Calc_Entrants
     note that transitions from behavior change are not taken into account *)
  bDoStockFlow  := False;

  (* Due to the use of t-1, the index of _EpidemicStYr relative to the first year of the projection,
     cannot be in first year if t loop starts at 2
     The use of a year delay for the Primary stage adds an extra year *)
//  if GetGBFinalYear(p) >= (GetGBCalcYear(p) + _EpidemicStYrIdx - GetGBCalcYearIdx(p)) then
  if (GetGBCalcYear(proj) <= _EpidemicStYr) and (GetGBFinalYear(proj) >= _EpidemicStYr) then
    if ((_EpidemicStYr - GetGBCalcYear(proj) + 1) < 3) then
    begin
      _EpidemicStYr := GetGBCalcYear(proj) + 2;
      HV.SetEpidemicStYr(_EpidemicStYr);
    end;

  _ts                   := 2*HV_TimeStep;
  bHasNewInfections     := False;
  (* Exit rate from primary stage *)
  _YearsInPrimaryStage  := _MonthsInPrimaryStage/12;

  Init_Veffect;
  Init_RiskGroupPercent;
  Init_Variables;
  Init_Disinhibition;
  Init_IntCD4Vars;
  Init_GetFirstYearART;
  Init_AdjCoverage;
  //Calc_BehavioralCoverage;
end;

Destructor THVProj.Destroy;
begin
  ReverseTransfer_Local;
end;

procedure THVProj.SetCD4Coverage(t, c: byte; value : double);
begin
_CD4Coverage[c,t]:=value;
end;

procedure THVProj.Transfer_Local;
var
 t,s,r,i,m  : byte;
 sum        : double;
begin
  if GetGBInCalcStateMode(p) then
    HV_RetrieveCalcStateData(HV);

  (* Transfer data to Local Arrays *)
  _EpidemicStYr := HV.GetEpidemicStYr;
  _TransMultSTI := HV.GetTransMultSTI;
  _TNType       := HV_GetRNTypeOfVaccine(p);
  _CondomEff    := HV.GetCondomEff/100;
  _InitialPulse := HV.GetInitialPulse;
  _TransMultMSM := HV.GetTransMultMSM;
  _ARTSurviving := 0.95; //GetHVARTSurviving(p);
  _Targeting    := HV_GetRNVaccineTargeting(p);
  _TransMultM   := HV.GetTransMultM;
  _TransHIVF    := HV.GetTrans;

  (* [CDP] this one should stay *)
  _MonthsInPrimaryStage := HV.GetMonthsInPrimaryStage;

  (* [CDP] balancing sex acts *)
  bBalanceSexActs := HV.GetBalanceSexActs;

  for i := HV_RN_AmongVacc to HV_RN_AmongAdults do
    _VaccineBehavEffect[i] := HV_GetRNVaccineBehavEffect(p,i)/100;

  for i := HV_Susceptibility to HV_Infect do
    _RedWhenCircum[i] := HV.GetRedWhenCircum(i)/100;

  for i := HV_PrimaryInf to HV_SympART do
    _Infectiousness[i] := HV.GetInfectiousness(i);

  for i := HV_RN_Efficacy to HV_RN_Duration do
    _VaccineEffectiveness[i] := HV_GetRNVaccineEffectiveness(p,i);

  for t := GetGBCalcYearIdx(p) to GetGBFinalYearIdx(p) do
  begin
    _PerIDUSharing[t] := HV.GetPerIDUSharing(t)/100;
    _MalesCircum[t] := HV_GetRNCoverage(p,t,HV_RN_MC15_49)/100;
    _ARTReceive[t] := HV.GetARTReceive(t)/100;

    for r := HV_None to HV_IDU_F3 do
      _STIPrev[r,t] := HV.GetSTIPrev(r,t);

    for r := HV_None to HV_MSMHR do
      _CondomPercent[r,t] := HV.GetCondomPercent(r,t);

    (* New PrEP method mix functionality *)
    for s := HV_Male to HV_Female do
      for r := HV_LRH to HV_MSMIDU do
        _PrEPPercent[s,r,t] := HV_GetRNPrEPCoverage(p,s,r,t)/100;

    for r := HV_LRH to HV_MSMIDU do
      for s := HV_Male to HV_Female do
      begin
        sum:=0;
        for m := HV_RN_PrEPOral to HV_RN_PrEPRing  do
          sum := sum + HV_GetRNPrEPEffectiveness(p,HV_RN_Effectiveness,m)/100*
                       HV_GetRNPrEPEffectiveness(p,HV_RN_Adherence,m)/100*
                       HV_GetRNMethodMix(p,s,r,m,t)/100;

        _PrEPEff[s,r,t]:=sum;
      end;

    for r := HV_AllRisk to HV_IDU_F3 do
    begin
      _NumPart[r,t] := HV.GetNumPart(r,t);
      _SexActs[r,t] := HV.GetSexActs(r,t);
    end;

    (* assume that IDU and _MSMIDU has MR and MSM_MR behavior respectively
       see also HV_Proj:CalcIDUBehavior() *)

     (* IDU
        PrEPPercent set in editor *)
    _STIPrev[HV_MSMIDU,t] := _STIPrev[HV_MRH,t];

     (* MSM:IDU *)
    _STIPrev[HV_MSMIDU,t] := _STIPrev[HV_MSMMR,t];

    for s := HV_Bothsexes to HV_Female do
    begin
      _AgeFirstSex[s,t] := HV.GetAgeFirstSex(s,t);
      _ForceInf[s,t] := HV.GetForceInf(s,t);
    end;

    for r := HV_AllRisk to HV_MSM_F3 do
      _VacCoverage[r,t] := HV_GetRNVacCoverage(p,r,t);
  end;

  for r := HV_AllRisk to HV_IDU_F3 do
    _PercMarried[r] := HV.GetPercMarried(r);

  for s := HV_Bothsexes to HV_Female do
    for r := HV_AllRisk to HV_MSMIDU do
      _IncRecruitment[s,r] := HV.GetIncRecruitment(s,r);

  for r := HV_AllRisk to HV_MSM_F3 do
  begin
    _Behavior[r,HV_AvgDur] := HV.GetBehavior(r,HV_AvgDur);
    if(not(_Behavior[r,HV_AvgDur]=0)) then
      _BehaviorChangeRate[r]:=1/_Behavior[r,HV_AvgDur];

    if(_BehaviorChangeRate[r]<1/100) then
        _BehaviorChangeRate[r]:=0;

    _Behavior[r,HV_PercPop] := HV.GetBehavior(r,HV_PercPop);
    (* CDP. Sometimes _Behavior[..,HV_PercPop] don't add to 100.
       Very small percentages slip into, for example, r=HV_IDU.*)
    if(_Behavior[r,HV_PercPop]<0.009) then
       _Behavior[r,HV_PercPop]:=0;
  end;
end;

procedure THVProj.ReverseTransfer_Local;
var
  t,s,r,h,v : byte;
begin
  for v := HV_RN_ALLVacc to HV_RN_NoProt do
  begin
    HV_SetRNVeffect(p,HV_RN_TakeAction,v,_Veffect[HV_RN_TakeAction,v]);
    HV_SetRNVeffect(p,HV_RN_DegreeAction,v,_Veffect[HV_RN_DegreeAction,v]);
  end;

  for t := GetGBCalcYearIdx(p) to GetGBFinalYearIdx(p) do
    for s := HV_BothSexes to HV_Female do
      for h := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
      begin
        HV.SetNewlyOnART(s,h,t,_NewlyOnART[s,h,t]);
        HV.SetNewlyEligibleART(s,h,t,_NewlyEligibleART[s,h,t]);
      end;

  for t := GetGBCalcYearIdx(p) to GetGBFinalYearIdx(p) do
  begin
    HV.SetTotalNewVaccinations(t,_TotalNewVaccinations[t]);
    HV.SetUnvaccinated(t,_Unvaccinated[t]);
    HV.SetIncidence(t,_Incidence[t]);
    HV.SetTotalAIDSDeaths(t,_TotalAIDSDeaths[t]);

    for s := HV_Bothsexes to HV_Female do
    begin
      HV.SetPopulations(s,t,_Populations[s,t]);
      HV.SetTotalVaccinated(s,t,_TotalVaccinated[s,t]);
      HV.SetNonAIDSDeathRate(s,t,_NonAIDSDeathRate[s,t]);
      HV.SetRateofAging(s,t,_RateofAging[s,t]);
      HV.SetTotalAdultsHIV(s,t,_TotalAdultsHIV[s,t]);
      (* [CDP] *)
      HV.SetTotalAdultsART(s,t,_TotalAdultsart[s,t]);

      for r := HV_AllRisk to HV_MSMIDU do
      begin
        HV.SetNewVaccinations(s,r,t,_NewVaccinations[s,r,t]);
        HV.SetCalcPrevalence(s,r,t,_CalcPrevalence[s,r,t]);
        HV.SetAIDSDeathsART(s,r,t,_AIDSDeathsART[s,r,t]);
        HV.SetrMultAll(s,r,t,_rMultAll[s,r,t]);
        (* [CDP] Exit Rate now depends on CD4 stage *)
        for h := HV_Primary to HV_CD4_LT50 do
            HV.SetExitRate(s,r,h,t,_ExitRate[s,r,h,t]);

        (* [CDP] not in h loop! *)
        for v := HV_RN_ALLVacc to HV_RN_NoProt do
          begin
            HV.SetNewInfections(s,r,v,t,_NewInfections[s,r,v,t]);
            HV.SetAIDSDeaths(s,r,v,t,_AIDSDeaths[s,r,v,t]);
          end; {v loop}

        (* [CDP] restructred *)
        for h := HV_Negative to HV_AllHIV do
          HV.SetVaccinated(s,r,h,t,_Vaccinated[s,r,h,t]);

        (* [CDP] *)
        for h := HV_Negative to HV_AllHIV do
          for v := HV_RN_ALLVacc to HV_RN_NoProt do
            HV.SetAdults(s,r,h,v,t,_Adults[s,r,h,v,t]);

      end;{r loop}
    end;{s loop}
  end;{t loop}

  for r := HV_AllRisk to HV_MSM_F3 do
    HV.SetRiskGroupPercent(r,_RiskGroupPercent[r]);

  for s := HV_Bothsexes to HV_Female do
    for r := HV_AllRisk to HV_MSMIDU do
      for h := HV_Negative to HV_AllHIV do
        for v := HV_RN_AllVacc to HV_RN_NoProt do
          Finalize(_AdultsTS[s,r,h,v]);

  if GetGBTempProj > 1 then
    HV_CaptureCalcStateData(HV);
end;

procedure THVProj.Init_Disinhibition;
var
  t,r,s,nr  : byte;
  df        : single;
begin
  if _VaccineBehavEffect[HV_RN_AmongVacc] > _VaccineBehavEffect[HV_RN_AmongAdults] then
    df := _VaccineBehavEffect[HV_RN_AmongVacc]
  else
    df := _VaccineBehavEffect[HV_RN_AmongAdults];

  for t := GetGBCalcYearIdx(p) to HV.GetFinalIndex do
  begin
    for r := HV_LRH to HV_MSMIDU do
    begin
      (* determines which button was selected in the editor, either single or
         different coverage targets *)
      if HV_GetRNVaccineCovType(p) = HV_RN_Single then
        nr := HV_AllRisk
      else
        nr := r;
      if r + HV_Total <= HV_IDU_F3 then
        if (_Behavior[r,HV_PercPop] + _Behavior[r + HV_Total,HV_PercPop]) > 0 then
          _CondomPercent[r,t] := _CondomPercent[r,t] + (_CondomPercent[r,1] - _CondomPercent[r,t]) *
                               (_VacCoverage[nr,t] * _Behavior[r,HV_PercPop] +
                                _VacCoverage[nr + HV_Total,t] * _Behavior[r + HV_Total,HV_PercPop])/
                               (_Behavior[r,HV_PercPop] + _Behavior[r + HV_Total,HV_PercPop])/100 * df;

      for s := HV_Male to HV_Female do
      begin
        if s = HV_Male then
        begin
          _NumPart[r,t] := _NumPart[r,t] + (_NumPart[r,1] - _NumPart[r,t]) *
                             _VacCoverage[nr,t]/100 * df;
          _SexActs[r,t] := _SexActs[r,t] + (_SexActs[r,1] - _SexActs[r,t]) *
                             _VacCoverage[nr,t]/100 * df;
        end
        else if r <= HV_HRH then
        begin
          _NumPart[r + HV_Total,t] := _NumPart[r + HV_Total,t] + (_NumPart[r + HV_Total,1] - _NumPart[r + HV_Total,t]) *
                             _VacCoverage[nr + HV_Total,t]/100 * df;
          _SexActs[r + HV_Total,t] := _SexActs[r + HV_Total,t] + (_SexActs[r + HV_Total,1] - _SexActs[r + HV_Total,t]) *
                             _VacCoverage[nr + HV_Total,t]/100 * df;
        end;
      end;
    end;
  end;
end;

//procedure THVProj.Init_WithValuesFromFP;
//begin
//  if GB_ModuleActive(p, GB_FP) and GetHVlowRiskCondomUseFromFP(p) then
//    HV_ResetValuesFromFP(p);
//end;

procedure THVProj.Init_Veffect;
begin
  _Veffect[HV_RN_TakeAction,HV_RN_Take]       := _VaccineEffectiveness[HV_RN_Efficacy]/100;
  _Veffect[HV_RN_TakeAction,HV_RN_Partial]    := 0;
  _Veffect[HV_RN_TakeAction,HV_RN_NoProt]     := 1-(_VaccineEffectiveness[HV_RN_Efficacy]/100);

  _Veffect[HV_RN_DegreeAction,HV_RN_Take]     := 0;
  _Veffect[HV_RN_DegreeAction,HV_RN_Partial]  := 1;
  _Veffect[HV_RN_DegreeAction,HV_RN_NoProt]   := 0;
end;

procedure THVProj.Init_RiskGroupPercent;
var
  value   : double;
  r       : byte;
begin
  for r := HV_None to HV_IDU_F3 do
  begin
    if not(r in [HV_MSM,HV_Total{HV_AllRisk_F3}]) then
    begin
      value := _Behavior[r,HV_PercPop]/100;
      _RiskGroupPercent[r] := value;
    end
    else
    if r = HV_MSM then
    begin
      _RiskGroupPercent[r] := _Behavior[r,HV_PercPop]/100;
      _RiskGroupPercent[HV_MSM_F3] := 0;
    end;
  end;
end;

procedure THVProj.Init_DynamicYearArrays;
var
  i, j, k, l : byte;
begin
  for i := HV_BothSexes to HV_Female do
    for j := HV_AllRisk to HV_MSMIDU do
      for k := HV_Negative to HV_AllHIV do
        for l := HV_RN_AllVacc to HV_RN_NoProt do
          SetLength(_Adults[i, j, k, l], GetGBFinalYearIdx(p) + 1);

  for i := HV_BothSexes to HV_Female do
    for j := HV_AllRisk to HV_MSMIDU do
      SetLength(_AddlRiskGrpARTCov[i, j], GetGBFinalYearIdx(p) + 1);

  for i := HV_BothSexes to HV_Female do
    for j := HV_AllRisk to HV_MSMIDU do
      for k := HV_Negative to HV_AllHIV do
        for l := HV_RN_AllVacc to HV_RN_NoProt do
          SetLength(_RiskGrpARTCov[i, j, k, l], GetGBFinalYearIdx(p) + 1);

  for i := HV_BothSexes to HV_Female do
    for j := HV_AllRisk to HV_MSMIDU do
      for k := HV_Negative to HV_AllHIV do
        for l := HV_RN_AllVacc to HV_RN_NoProt do
          SetLength(_HV_RiskGrpARTPerc[i, j, k, l], GetGBFinalYearIdx(p) + 1);

  for i := DP_Male to DP_Female do
    for j := DP_CD4_GT500 to DP_CD4_LT50 do
      SetLength(_DP_RiskGrpARTPerc[i, j], GetGBFinalYearIdx(p) + 1);

  for i := HV_BothSexes to HV_Female do
    SetLength(_TotalAdultsHIV[i], GetGBFinalYearIdx(p) + 1);

  for i := HV_BothSexes to HV_Female do
    SetLength(_TotalAdultsART[i], GetGBFinalYearIdx(p) + 1);

  for i := HV_BothSexes to HV_Female do
    SetLength(_Populations[i], GetGBFinalYearIdx(p) + 1);

  for i := HV_BothSexes to HV_Female do
    for j := HV_AllRisk to HV_MSMIDU do
      for k := 1 to 2 do
        SetLength(_TotalPops[i, j, k], GetGBFinalYearIdx(p) + 1);

  for i := HV_BothSexes to HV_Female do
    for j := HV_AllRisk to HV_MSMIDU do
      for k := 1 to 2 do
        SetLength(_TotalEntrants[i, j, k], GetGBFinalYearIdx(p) + 1);

  for i := HV_BothSexes to HV_Female do
    for j := HV_AllRisk to HV_MSMIDU do
      for k := 1 to 2 do
        SetLength(_TotalExits[i, j, k], GetGBFinalYearIdx(p) + 1);

  for i := HV_None to HV_IDU_F3 do
    SetLength(_STIPrev[i], GetGBFinalYearIdx(p) + 1);

  SetLength(_MalesCircum, GetGBFinalYearIdx(p) + 1);

  for i := HV_None to HV_MSMIDU do
    SetLength(_CondomPercent[i], GetGBFinalYearIdx(p) + 1);

  for i := HV_AllRisk to HV_MSM_F3 do
    SetLength(_NumPart[i], GetGBFinalYearIdx(p) + 1);

  for i := HV_AllRisk to HV_MSM_F3 do
    SetLength(_SexActs[i], GetGBFinalYearIdx(p) + 1);

  for i := HV_BothSexes to HV_Female do
    SetLength(_AgeFirstSex[i], GetGBFinalYearIdx(p) + 1);

  SetLength(_PerIDUSharing, GetGBFinalYearIdx(p) + 1);

  for i := 1 to HV_RN_MaxInterventions do
    SetLength(_AdjCoverage[i], GetGBFinalYearIdx(p) + 1);

  SetLength(_Unvaccinated, GetGBFinalYearIdx(p) + 1);

  for i := HV_BothSexes to HV_Female do
    for j := HV_AllRisk to HV_MSMIDU do
      SetLength(_NewVaccinations[i, j], GetGBFinalYearIdx(p) + 1);

  for i := HV_RN_AllRisk to HV_RN_MSM_F do
    SetLength(_VacCoverage[i], GetGBFinalYearIdx(p) + 1);

  for i := HV_BothSexes to HV_Female do
    for j := HV_AllRisk to HV_MSMIDU do
      for k := HV_Negative to HV_AllHIV do
        SetLength(_Vaccinated[i, j, k], GetGBFinalYearIdx(p) + 1);

  SetLength(_TotalNewVaccinations, GetGBFinalYearIdx(p) + 1);

  for i := HV_BothSexes to HV_Female do
    SetLength(_TotalVaccinated[i], GetGBFinalYearIdx(p) + 1);

  for i := HV_Male to HV_Female do
    for j := HV_LRH to HV_MSMIDU do
      SetLength(_PrEPPercent[i, j], GetGBFinalYearIdx(p) + 1);

  for i := HV_Male to HV_Female do
    for j := HV_LRH to HV_MSMIDU do
      SetLength(_PrEPEff[i, j], GetGBFinalYearIdx(p) + 1);

  for i := HV_Male to HV_Female do
    for j := HV_DP_H_No to HV_AllHIV do
      SetLength(_Pop15rate[i, j], GetGBFinalYearIdx(p) + 1);

  SetLength(_TotalNewAdults, GetGBFinalYearIdx(p) + 1);

  SetLength(_MigrNum, GetGBFinalYearIdx(p) + 1);
  SetLength(_MigrRate, GetGBFinalYearIdx(p) + 1);
  SetLength(_TotalMigrants, GetGBFinalYearIdx(p) + 1);

  for i := HV_Male to HV_Female do
    for j := HV_LRH to HV_MSMIDU do
      SetLength(_PrEPEff[i, j], GetGBFinalYearIdx(p) + 1);

  for i := HV_BothSexes to HV_Female do
    SetLength(_ForceInf[i], GetGBFinalYearIdx(p)+ 1);

  for i := HV_BothSexes to HV_Female do
    for j := HV_AllRisk to HV_MSMIDU do
      SetLength(_rMultAll[i, j], GetGBFinalYearIdx(p) + 1);

  for i := HV_BothSexes to HV_Female do
    for j := HV_AllRisk to HV_MSMIDU do
      SetLength(_CalcPrevalence[i, j], GetGBFinalYearIdx(p) + 1);

  for i := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
    SetLength(_MultART[i], GetGBFinalYearIdx(p) + 1);

  for i := HV_BothSexes to HV_Female do
    for j := HV_AllRisk to HV_MSMIDU do
      for k := HV_RN_AllVacc to HV_RN_NoProt do
        SetLength(_NewInfections[i, j, k], GetGBFinalYearIdx(p) + 1);

  SetLength(_Incidence, GetGBFinalYearIdx(p) + 1);
  SetLength(_IncSexRatio, GetGBFinalYearIdx(p) + 1);
  SetLength(_TotalNewInfection, GetGBFinalYearIdx(p) + 1);

  for i := HV_BothSexes to HV_Female do
    SetLength(_NonAIDSDeathRate[i], GetGBFinalYearIdx(p) + 1);

  for i := HV_BothSexes to HV_Female do
    for j := HV_AllRisk to HV_MSMIDU do
      for k := HV_RN_AllVacc to HV_RN_NoProt do
        SetLength(_AIDSDeaths[i, j, k], GetGBFinalYearIdx(p) + 1);

  for i := HV_BothSexes to HV_Female do
    for j := HV_AllRisk to HV_MSMIDU do
      SetLength(_AIDSDeathsART[i, j], GetGBFinalYearIdx(p) + 1);

  SetLength(_TotalAIDSDeaths, GetGBFinalYearIdx(p) + 1);

  SetLength(_TotalARTDeaths, GetGBFinalYearIdx(p) + 1);

  for i := HV_BothSexes to HV_Female do
    SetLength(_RateofAging[i], GetGBFinalYearIdx(p) + 1);

  for i := HV_BothSexes to HV_Female do
    for j := HV_Primary to HV_CD4_LT50 do
      SetLength(_RateofHIVAging[i, j], GetGBFinalYearIdx(p) + 1);

  for i := HV_BothSexes to HV_Female do
    for j := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
      SetLength(_RateofARTAging[i, j], GetGBFinalYearIdx(p) + 1);

  for i := HV_BothSexes to HV_Female do
    for j := HV_AllRisk to HV_MSMIDU do
      for k := HV_Primary to HV_CD4_LT50 do
        SetLength(_ExitRate[i, j, k], GetGBFinalYearIdx(p) + 1);

  for i := HV_BothSexes to HV_Female do
    for j := HV_AllRisk to HV_MSMIDU do
      for k := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
        SetLength(_ExitRateART[i, j, k], GetGBFinalYearIdx(p) + 1);

  for i := HV_BothSexes to HV_Female do
    for j := HV_AllRisk to HV_MSMIDU do
      for k := HV_Primary to HV_CD4_LT50 do
        SetLength(_ProgressRate[i, j, k], GetGBFinalYearIdx(p) + 1);

  SetLength(_ARTReceive, GetGBFinalYearIdx(p) + 1);

  for i := HV_BothSexes to HV_Female do
    for j := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
      SetLength(_NewlyOnART[i, j], GetGBFinalYearIdx(p) + 1);

  for i := HV_BothSexes to HV_Female do
    for j := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
      SetLength(_NewlyEligibleART[i, j], GetGBFinalYearIdx(p) + 1);
end;

procedure THVProj.Init_Variables;
var
  t,s,r,a,h,v,nr  : byte;
  value, value2   : double;
  TotalPop        : HV_MaleFemale;
begin
  TotalPop[HV_Male] := 0;
  TotalPop[HV_Female] := 0;

  for t := GetGBCalcYearIdx(p) to GetGBFinalYearIdx(p) do
    for s := HV_Male to HV_Female do
      for h := HV_DP_H_No to HV_AllHIV do
        _Pop15rate[s, h, t]:=0;

  for s := HV_Male to HV_Female do
  begin
    for h := HV_Primary to HV_CD4_LT50 do
    begin
      _mu[h,s]:=0;
      _lambda[h,s]:=0;
      _ScaleMuART[h,s]:=0;
    end;

    for h := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
      _alpha[h,s]:=0;
  end;

  for s :=  HV_Bothsexes to HV_Female do
    for r := HV_AllRisk to HV_MSMIDU do
      for t := GetGBCalcYearIdx(p) to GetGBFinalYearIdx(p) do
      begin
        (* [CDP] Note range for h loops throughout *)
        for h := HV_Negative to HV_AllHIV do
          for v := HV_RN_AllVacc to HV_RN_NoProt do
          begin
            _Adults[s,r,h,v,t] := 0;
            _AIDSDeaths[s,r,v,t] :=0;
            _AIDSDeathsART[s,r,t] :=0;
          end;
      end;

  if GetGBInCalcStateMode(p) then
  begin
  for s :=  HV_Bothsexes to HV_Female do
    for r := HV_AllRisk to HV_MSMIDU do
      (* [CDP] Note range for h loops throughout *)
      for h := HV_Negative to HV_AllHIV do
        for v := HV_RN_AllVacc to HV_RN_NoProt do
          _Adults[s,r,h,v,1] := HV.GetAdults(s,r,h,v,1)
  end
  else
  begin
    for s := HV_Male to HV_Female do
    for a := HV_DP_A15_19 to HV_DP_A45_49 do
      TotalPop[s] := TotalPop[s] + HV_GetDPPop(p, 1, s, a);

    for s := HV_Male to HV_Female do
      for r := HV_None to HV_MSMIDU do
      begin
        if s = HV_Male then
          nr := r
        else
          nr := r + HV_Total;
        if nr <= HV_IDU_F3 then
        begin
          value := TotalPop[s] *  _Behavior[nr,HV_PercPop]/100;
          _Adults[s,r,HV_Negative,HV_RN_UnV,1] := value;
        end;
      end;

    for s := HV_Male to HV_Female do
      for r := HV_None to HV_MSMIDU do
        for h := HV_Negative to HV_CD4_LT50 do
        begin
          value := 0;
          for v := HV_RN_UnV to HV_RN_NoProt do
          begin
            if not((s = HV_Female) and (r >= HV_MSM)) then
              value := value + _Adults[s,r,h,v,1];
          end;
          _Adults[s,r,h,HV_RN_AllVacc,1] := value;
        end;
  end;

  for s := HV_Male to HV_Female do
  begin
    value := 0;
    for r := HV_None to HV_MSMIDU do
    begin
      for h := HV_Negative to HV_CD4_LT50_ART do
      begin
        for v := HV_RN_UnV to HV_RN_NoProt do
        begin
          value := value + _Adults[s,r,h,v,1];
        end;{v loop}
      end;{h loop}
    end;
      _Populations[s,1] := value;
  end;
  value := _Populations[HV_Male,1] + _Populations[HV_Female,1];
  _Populations[HV_Bothsexes,1] := value;

  for s := HV_Male to HV_Female do
    for r := HV_None to HV_MSMIDU do
      (* [CDP check the role of h in this loop!] *)
      for h := HV_Negative to HV_CD4_LT50 do
      begin
        if not((s = HV_Female) and (r >= HV_MSM)) then
        begin
          (* [CDP change to new CD4 ranges. ART receivers not to be vaccinated] *)
          value := (_Adults[s,r,HV_Negative,HV_RN_UnV,1]
                     + _Adults[s,r,HV_CD4_GT500,HV_RN_UnV,1]
                     + _Adults[s,r,HV_CD4_350_500,HV_RN_UnV,1]
                     + _Adults[s,r,HV_CD4_250_349,HV_RN_UnV,1]
                     + _Adults[s,r,HV_CD4_100_199,HV_RN_UnV,1]
                     + _Adults[s,r,HV_CD4_50_99,HV_RN_UnV,1]
                     + _Adults[s,r,HV_CD4_LT50,HV_RN_UnV,1]);
          if value <> 0 then
          begin
            if HV_GetRNVaccineCovType(p) = HV_RN_Single then
              nr := HV_AllRisk
            else
            begin
              if s = HV_Female then
                nr := r + HV_Total
              else
                nr := r;
            end;
            //if nr <= HV_IDU_F_VC new>HV_IDU_F_RG then
            (* [CDP] new constants in HVCONST, using TNCONST instead *)
            if nr <= HV_IDU_F3 then
              value := ((_VacCoverage[nr,1]/100) * _Adults[s,r,h,HV_RN_AllVacc,1]
                     - _Adults[s,r,h,HV_RN_Take,1]
                     - _Adults[s,r,h,HV_RN_Partial,1]
                     - _Adults[s,r,h,HV_RN_NoProt,1] )
                    * _Adults[s,r,h,HV_RN_UnV,1]
                    / value;
          end;
          if value < 0 then value := 0;
          _NewVaccinations[s,r,1] := value;
        end;
      end;

  for s := HV_Male to HV_Female do
  begin
    value := 0;
    for r := HV_None to HV_MSMIDU do
    begin
      if (s = HV_Female) and (r >= HV_MSM) then
      else
      value := value + _NewVaccinations[s,r,1];
    end;

    _TotalVaccinated[s,1] := value;
  end;
  value := _TotalVaccinated[HV_Male,1] + _TotalVaccinated[HV_Female,1];
  _TotalVaccinated[HV_Bothsexes,1] := value;
  value := _Populations[HV_Bothsexes,1] - value;
  _Unvaccinated[1] := value;

  for s := HV_Male to HV_Female do
  begin
    value2 := 0;
    for r := HV_None to HV_MSMIDU do
    begin
      value := 0;
      for h := HV_Negative to HV_CD4_LT50_ART do
      begin
        for v := HV_RN_UnV to HV_RN_NoProt do
        begin
          if not((s = HV_Female) and (r >= HV_MSM)) then
          value := value + _Adults[s,r,h,v,1];
        end;
      end;{h loop}
      _Adults[s,r,HV_AllHIV,HV_RN_AllVacc,1] := value;
      value2 := value2 + value;
    end;{r loop}
    _Adults[s,HV_AllRisk,HV_AllHIV,HV_RN_AllVacc,1] := value2;
  end; {s loop}

 {* Initialize some of the variables - should be moved elsewhere *}
  for t := GetGBCalcYearIdx(p) to GetGBFinalYearIdx(p) do
    for s := HV_Bothsexes to HV_Female do
      for r := HV_AllRisk to HV_MSMIDU do
        begin
          for h := HV_Primary to HV_CD4_LT50 do
          begin
            _ProgressRate[s,r,h,t] := 0;
            _ExitRate[s,r,h,t] := 0;
          end;
          for h :=  HV_CD4_GT500_ART to HV_CD4_LT50_ART do
            _ExitRateART[s,r,h,t] := 0
        end;{for r,s,t}

  for s := HV_Bothsexes to HV_Female do
    for r := HV_AllRisk to HV_MSMIDU do
      for v := HV_RN_AllVacc to HV_RN_NoProt do
      begin
        _NewInfectionsT[s,r,v] := 0;
        for t := GetGBCalcYearIdx(p) to GetGBFinalYearIdx(p) do
          _NewInfections[s,r,v,t]:= 0;
      end;

  for t := GetGBCalcYearIdx(p) to GetGBFinalYearIdx(p) do
  begin
    _TotalNewInfection[t] := 0;
    _ARTReceive[t] := 0;
    _TotalAIDSDeaths[t] := 0;
    _TotalARTDeaths[t] := 0;
    _TotalNewAdults[t] := 0;
    _TotalMigrants[t] := 0;
  end;

  for s := HV_Bothsexes to HV_Female do
    for r := HV_AllRisk to HV_MSMIDU do
      for h := HV_Negative to HV_AllHIV do
        for v := HV_RN_AllVacc to HV_RN_NoProt do
          SetLength(_AdultsTS[s,r,h,v], (GB_GetNumDynYears(p) * HV_TimeStep) + 1);

  for s := HV_Bothsexes to HV_Female do
    for h := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
      for t := GetGBCalcYearIdx(p) to GetGBFinalYearIdx(p) do
      begin
         _NewlyOnART[s,h,t] := 0;
         _NewlyEligibleART[s,h,t] := 0;
      end;

  for t := GetGBCalcYearIdx(p) to GetGBFinalYearIdx(p) do
    for s := HV_Bothsexes to HV_Female do
       for r := HV_AllRisk to HV_MSMIDU do
       begin
         _TotalPops[s,r,1,t] := 0;
         _TotalPops[s,r,2,t] := 0;
         (* not using dimension 2 on these two *)
         _TotalEntrants[s,r,1,t] := 0;
         _TotalExits[s,r,1,t] := 0;
         _rMultAll[s,r,t] := 1;
       end;

  for s := HV_Male to HV_Female do
    for r := HV_None to HV_MSMIDU do
      for h := HV_Negative to HV_CD4_LT50_ART do
        for v := HV_RN_UnV to HV_RN_NoProt do
          _Migrant[s,r,h,v] := 0;

   for h := HV_CD4_GT500 to HV_CD4_LT50 do
     for s := HV_Male to HV_Female do
        for r := HV_None to HV_MSMIDU do
          for v := HV_RN_UnV to HV_RN_NoProt do
            for t := HV_min_year to GetGBFinalYearIdx(p) do
            begin
              _AddlRiskGrpARTCov[s,r,t] := 0;
              _RiskGrpARTCov[s,r,h,v,t] := 0;
              _HV_RiskGrpARTPerc[s,r,h,v,t] := 0;
            end;

   for h := HV_CD4_GT500 to HV_CD4_LT50 do
     for s := HV_Male to HV_Female do
           for t := HV_min_year to GetGBFinalYearIdx(p) do
              _DP_RiskGrpARTPerc[s,h,t] := 0

end;

procedure THVProj.Calc_NonAIDSMortAIM(t : integer);
var
  s,a : byte;
  AdultDeaths, AdultPop : double;
begin

  for s := HV_Male to HV_Female do
  begin
    _NonAIDSDeathRate[s,t]:=0;
  end;

  for s := HV_Male to HV_Female do
  begin
    AdultPop := HV_GetDPPop15to49(p,t,s);
    AdultDeaths := 0;
    for a := HV_DP_A15_19 to HV_DP_A45_49 do
    begin
      (* [CDP] GetDPDeaths(p, t, a, s) does not include AIDSdeaths! *)
      AdultDeaths := AdultDeaths + HV_GetDPDeaths(p, t, a, s);//  - GetDPAIDSdeaths(p, t, a, s);
      //AdultPop := AdultPop + HV_GetDPPop(p,t,s,a);
    end;

    if(AdultPop>0) then
      _NonAIDSDeathRate[s,t] := AdultDeaths/AdultPop;
  end;
end;

procedure THVProj.Calc_HIVNegAgingAIM(t,_ts : integer);
var
  s,r       : byte;
  numerator,
  denominator : double;
begin
  for s := HV_Male to HV_Female do
  begin
    numerator:=HV_GetDPAIDSAging50(p, t, s, HV_DP_H_No, HV_DP_D_HIVNeg);

    denominator:=0;
    for r := HV_None to HV_MSMIDU do
    //for v := HV_RN_UnV to HV_NoProt do
      denominator:= denominator+_AdultsTS[s,r,HV_Negative,HV_RN_AllVacc,_ts];

    _RateofAging[s,t]:=0;
    if((numerator>0) and (denominator>0)) then
     _RateofAging[s,t]:=numerator/denominator;
  end;
end;

//procedure THVProj.Calc_InitAgingRates;
//var
//  s,h,a       : byte;
//  t           : integer;
//  AdultPop,
//  value       : double;
//begin
//  t := 1;
//  for s := HV_Male to HV_Female do
//  begin
//    AdultPop := 0;
//    for a := HV_DP_A15_19 to HV_DP_A45_49 do
//       AdultPop := AdultPop + HV_GetDPPop(p,t,s,a);
//
//    value := (HV_GetDPPop(p, t, s,HV_DP_A45_49) +
//                          HV_GetDPPop(p, t, s, HV_DP_A50_54)) / 10 / AdultPop;
//    _RateofAging[s,t] := value;
//  end;
//
//  for t := 2 to GetHVFinalIndex(p) do
//    for s := HV_Male to HV_Female do
//      _RateOfAging[s,t] := _RateOfAging[s,1];
//
//  for t := GetGBCalcYearIdx(p) to GetGBFinalIdx(p) do
//    for s := HV_Male to HV_Female do
//    begin
//      for h := HV_Primary to HV_CD4_LT50 do
//        _RateofHIVAging[s,h,t] := 0;
//
//      for h := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
//        _RateofARTAging[s,h,t] := 0;
//    end;
//end;

procedure THVProj.Calc_HIVPosAgingAIM(t,_ts: integer);
var
  s,h,r,d       : byte;
  numerator,
  denominator   : double;
begin
  for s := HV_Male to HV_Female do
  begin
    for h := HV_CD4_GT500 to HV_CD4_LT50 do
    begin
     (* No ART *)
      numerator:=HV_GetDPAIDSAging50(p, t, s, h, HV_DP_D_HIVPos);

      denominator:=0;
      for r := HV_None to HV_MSMIDU do
      //for v := HV_RN_UnV to HV_NoProt do
      begin
        denominator:= denominator+_AdultsTS[s,r,h,HV_RN_AllVacc,_ts];
        if(h=HV_CD4_GT500) then
          denominator:= denominator+_AdultsTS[s,r,HV_Primary,HV_RN_AllVacc,_ts];
      end;

       _RateofHIVAging[s,h,t]:=0;
      if((numerator>0) and (denominator>0)) then
         _RateofHIVAging[s,h,t]:=numerator/denominator;

      (* ART *)
      numerator:=0;
      for d := HV_DP_D_ARTlt6m to HV_DP_D_ARTgt12m do
        numerator:=numerator+HV_GetDPAIDSAging50(p, t, s, h, d);

      denominator:=0;
      for r := HV_None to HV_MSMIDU do
      //for v := HV_RN_UnV to HV_NoProt do
          denominator:= denominator+_AdultsTS[s,r,h+10,HV_RN_AllVacc,_ts];

      _RateofARTAging[s,h+10,t]:=0;
      if((numerator>0) and (denominator>0)) then
        _RateofARTAging[s,h+10,t]:=numerator/denominator;
    end;{h loop}

    (*No ART Primary *)
    _RateofHIVAging[s,HV_Primary,t]:=_RateofHIVAging[s,HV_CD4_GT500,t];
  end;
end;

procedure THVProj.Calc_Incidence;
var
 t            : byte;
 value,
 denominator  : double;
begin
  for t := GetGBCalcYearIdx(p) to HV.GetFinalIndex do
  begin
    denominator:=(_Populations[HV_Bothsexes,t]-(_TotalAdultsHIV[HV_Male,t]+_TotalAdultsHIV[HV_Female,t]));
    value:=0;
    if(denominator>0) then
      value := _TotalNewInfection[t]/denominator;
    _Incidence[t] := value;
  end;
end;

procedure THVProj.Calc_IncSexRatio(t: integer);
var
  r,v               : byte;
  value_m,value_f,
  denominator_m,
  denominator_f     : double;
begin
  (* men *)
  value_m:=0;
  for r := HV_None to HV_MSMIDU do
  for v := HV_RN_UnV to HV_RN_NoProt do
    value_m:=value_m+_NewInfections[HV_Male,r,v,t];

  denominator_m:=_Populations[HV_Male,t]-_TotalAdultsHIV[HV_Male,t];
  if(denominator_m>0) then
    value_m:=value_m/denominator_m;

  (* women *)
  value_f:=0;
  for r := HV_None to HV_MSMIDU do
  for v := HV_RN_UnV to HV_RN_NoProt do
    value_f:=value_f+_NewInfections[HV_Female,r,v,t];

  denominator_f:=_Populations[HV_Female,t]-_TotalAdultsHIV[HV_Female,t];
  if(denominator_f>0) then
    value_f:=value_f/denominator_f;

  (* incidence ratio *)
  _IncSexRatio[t]:=1;
  if(value_m>0) then
    _IncSexRatio[t]:= value_f/value_m;

  HV.SetIncSexRatio(t,_IncSexRatio[t]);
end;

procedure THVProj.Calc_PercentPop;
var
  t,s   : byte;
  value : double;
begin
  for s := HV_Male to HV_Female do
  begin
    for t := GetGBCalcYearIdx(p) to HV.GetFinalIndex do
    begin
      if _Populations[HV_Bothsexes,t] = 0 then
        value := 0
      else
        value := _Populations[s,t]/_Populations[HV_Bothsexes,t];
      HV.SetPercentPop(s,t,value);
    end;
  end;
end;

procedure THVProj.Calc_SumPopGroups(t:integer);
(* Calculate totals over dimensions of _Adults *)
var
  s,r,h,v       : byte;
  value,value2  : double;
begin
  for s := HV_Male to HV_Female do
    for r := HV_None to HV_MSMIDU do
      for h := HV_Negative to HV_CD4_LT50_ART do
      begin
        value := 0;
        for v := HV_RN_UnV to HV_RN_NoProt do
        begin
          if not((s = HV_Female) and (r >= HV_MSM)) then
            value := value + _Adults[s,r,h,v,t];
        end;
        _Adults[s,r,h,HV_RN_AllVacc,t] := value;
      end;

  for s := HV_Male to HV_Female do
  begin
    value2 := 0;
    for r := HV_AllRisk to HV_MSMIDU do
    begin
      value := 0;
      for h := HV_Negative to HV_CD4_LT50_ART do
      begin
        for v := HV_RN_UnV to HV_RN_NoProt do
        begin
          if not((s = HV_Female) and (r >= HV_MSM)) then
            value := value + _Adults[s,r,h,v,t];
        end;{v loop}
      end;{h loop}
      _Adults[s,r,HV_AllHIV,HV_RN_AllVacc,t] := value;
      value2 := value2 + value;
    end;
   _Adults[s,HV_AllRisk,HV_AllHIV,HV_RN_AllVacc,t] := value2;
  end;{s loop}

  for s := HV_Male to HV_Female do
    for r := HV_None to HV_MSMIDU do
      for v := HV_RN_UnV to HV_RN_NoProt do
      begin
        value := 0;
        for h := HV_Negative to HV_CD4_LT50_ART do
        begin
          if not((s = HV_Female) and (r >= HV_MSM)) then
            value := value + _Adults[s,r,h,v,t];
        end;
        _Adults[s,r,HV_AllHIV,v,t] := value;
      end;

  for s := HV_Male to HV_Female do
    for r := HV_None to HV_MSMIDU do
    begin
      value := 0;
      for h := HV_Negative to HV_CD4_LT50_ART do
        value := value + _Adults[s,r,h,HV_RN_AllVacc,t];
      _Adults[s,r,HV_AllHIV,HV_RN_AllVacc,t] := value;
    end;

  (* total CD4 categories over all risk groups for CD4 display *)
  for s := HV_Male to HV_Female do
    for h := HV_Negative to HV_CD4_LT50_ART do
    begin
      value := 0;
      for r := HV_None to HV_MSMIDU do
      begin
        for v := HV_RN_UnV to HV_RN_NoProt do
        begin
          if not((s = HV_Female) and (r >= HV_MSM)) then
            value := value + _Adults[s,r,h,v,t];
        end;
      end;
      _Adults[s,HV_AllRisk,h,HV_RN_AllVacc,t] := value;
    end;

//  (* total newly on art over cd4 and sex *)
//   for s := HV_Male to HV_Female do
//   begin
//     value:=0;
//     for h := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
//        value:=value+_NewlyOnART[s,h,t];
//     _NewlyOnART[s,HV_AllHIV_ART,t]:=value;
//   end;
//   _NewlyOnART[HV_Bothsexes,HV_AllHIV_ART,t]:=_NewlyOnART[HV_Male,HV_AllHIV_ART,t]+_NewlyOnART[HV_Female,HV_AllHIV_ART,t];
//
//   (* total newly eligible for art over cd4 and sex *)
//   for s := HV_Male to HV_Female do
//   begin
//     value:=0;
//     for h := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
//        value:=value+_NewlyEligibleART[s,h,t];
//     _NewlyEligibleART[s,HV_AllHIV_ART,t]:=value;
//   end;
//   _NewlyEligibleART[HV_Bothsexes,HV_AllHIV_ART,t]:=_NewlyEligibleART[HV_Male,HV_AllHIV_ART,t]+_NewlyEligibleART[HV_Female,HV_AllHIV_ART,t];
end;

procedure THVProj.Calc_HIVExitRates(t:integer);
(* calculate  not on ART and on ART exit and progress rates *)
var
  s,r,h,nr  : byte;
  value, numer, denom : double;
begin
  for s := HV_Male to HV_Female do
  begin
    (* Mortality reductions off ART proportional to ART coverage *)
    _ScaleMuART[HV_Primary,s] := 1.0;
    for h := HV_CD4_GT500 to HV_CD4_LT50 do
    begin
      numer := _Adults[s,HV_AllRisk,h   ,HV_RN_AllVacc,t];
      denom := _Adults[s,HV_AllRisk,h+10,HV_RN_AllVacc,t] + numer;
      if denom > 0
      then _ScaleMuART[h,s] := numer / denom
      else _ScaleMuART[h,s] := 1.0;
    end;

    for r := HV_None to HV_MSMIDU do
    begin
      if s = HV_Male then
        nr := r
      else
        nr := r + HV_Total;

      if (r >= HV_MSM) and (s = HV_Female)then
      begin
        for h := HV_Primary to HV_CD4_LT50 do
        begin
          _ExitRate[s,r,h,t] := 0;
          _ProgressRate[s,r,h,t] :=0;
        end;{h loop}
        for h := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
        begin
          _ExitRateART[s,r,h,t] := 0;
        end;{h loop}
       end{if r >= HV_MSM ...}
      else
      begin
        for h := HV_Primary to HV_CD4_LT50 do
        begin
          value := _NonAIDSDeathRate[s,t] +
                   _RateofHIVAging[s,h,t] +
                   _BehaviorChangeRate[nr]+
                   _mu[h,s]*_ScaleMuART[h,s];

          _ExitRate[s,r,h,t] := value;
          _ProgressRate[s,r,h,t] := _lambda[h,s];
        end;{h loop}

        for h := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
        begin
          value := _NonAIDSDeathRate[s,t] +
                   _RateofARTAging[s,h,t] +
                   _BehaviorChangeRate[nr]+_alpha[h,s];

          _ExitRateART[s,r,h,t] := value ;
        end;{h loop}
      end; { else }
    end; { r loop}
  end; {s loop}
end;

procedure THVProj.Calc_HIVExitsProgressors(t:integer);
(* calculate not on ART and on ART exit and progress rates *)
var
  s,r,h,v   : byte;
begin
  (* Keep track of stage exits and progressors *)
  for s := HV_Male to HV_Female do
    for r := HV_None to HV_MSMIDU do
      for v := HV_RN_UnV to HV_RN_NoProt do
      begin
        for h := HV_Primary to HV_CD4_LT50 do
        begin
            _StageExits[s,r,h,v] := _Adults[s,r,h,v,t]*(_ExitRate[s,r,h,t]+_ProgressRate[s,r,h,t]);
            _StageProgressors[s,r,h,v]:=_Adults[s,r,h,v,t]*(_ProgressRate[s,r,h,t]);
        end;{h loop}

        for h := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
            _StageExits[s,r,h,v] := _Adults[s,r,h,v,t]*(_ExitRateART[s,r,h,t]);
      end;{s,r,v loop}
end;

procedure THVProj.Calc_CD4Trans(t: integer);
var
  s,h,a           : byte;
  stageDuration   : double;
  AgePops,
  AgePopsART      : array[HV_Male .. HV_Female, HV_DP_CD4_15_24..HV_DP_CD4_45_54] of double;
  AgePopsTotal,
  AgePopsTotalART : array[HV_Male .. HV_Female] of double;
begin
  {* Age dependent transition paramters should be weighted average over age groups *}
  for s := HV_Male to HV_Female do
  begin
    AgePopsTotal[s]:=0;
    AgePopsTotalART[s]:=0;
    for a := HV_DP_CD4_15_24 to HV_DP_CD4_45_54 do
    begin
      AgePops[s,a]:=0;
      AgePopsART[s,a]:=0;
    end;
  end;

  for s := HV_Male to HV_Female do
  begin
    (* Not on ART
       get total in 4 age groups *)
//    AgePops[s,HV_DP_CD4_15_24]:=HV_GetDPHIV(p, t, s, HV_DP_A15_19)-HV_GetDPOnART(p, t, s, HV_DP_A15_19)+
//                                HV_GetDPHIV(p, t, s, HV_DP_A20_24)-HV_GetDPOnART(p, t, s, HV_DP_A20_24);
//    AgePops[s,HV_DP_CD4_25_34]:=HV_GetDPHIV(p, t, s, HV_DP_A25_29)-HV_GetDPOnART(p, t, s, HV_DP_A25_29)+
//                                HV_GetDPHIV(p, t, s, HV_DP_A30_34)-HV_GetDPOnART(p, t, s, HV_DP_A30_34);
//    AgePops[s,HV_DP_CD4_35_44]:=HV_GetDPHIV(p, t, s, HV_DP_A35_39)-HV_GetDPOnART(p, t, s, HV_DP_A35_39)+
//                                HV_GetDPHIV(p, t, s, HV_DP_A40_44)-HV_GetDPOnART(p, t, s, HV_DP_A40_44);
//    AgePops[s,HV_DP_CD4_45_54]:=HV_GetDPHIV(p, t, s, HV_DP_A45_49)-HV_GetDPOnART(p, t, s, HV_DP_A45_49)+
//                                HV_GetDPHIV(p, t, s, HV_DP_A50_54)-HV_GetDPOnART(p, t, s, HV_DP_A50_54);

    AgePops[s,HV_DP_CD4_15_24]:=HV_GetDPPopSexAgeHIV(p, t, s, HV_DP_A15_19, DP_Pop_H_NoART)+
                                HV_GetDPPopSexAgeHIV(p, t, s, HV_DP_A20_24, DP_Pop_H_NoART);
    AgePops[s,HV_DP_CD4_25_34]:=HV_GetDPPopSexAgeHIV(p, t, s, HV_DP_A25_29, DP_Pop_H_NoART)+
                                HV_GetDPPopSexAgeHIV(p, t, s, HV_DP_A30_34, DP_Pop_H_NoART);
    AgePops[s,HV_DP_CD4_35_44]:=HV_GetDPPopSexAgeHIV(p, t, s, HV_DP_A35_39, DP_Pop_H_NoART)+
                                HV_GetDPPopSexAgeHIV(p, t, s, HV_DP_A40_44, DP_Pop_H_NoART);
    AgePops[s,HV_DP_CD4_45_54]:=HV_GetDPPopSexAgeHIV(p, t, s, HV_DP_A45_49, DP_Pop_H_NoART)+
                                HV_GetDPPopSexAgeHIV(p, t, s, HV_DP_A50_54, DP_Pop_H_NoART);


    (* overall total *)
    AgePopsTotal[s]:= AgePops[s,HV_DP_CD4_15_24]+AgePops[s,HV_DP_CD4_25_34]+
                      AgePops[s,HV_DP_CD4_35_44]+AgePops[s,HV_DP_CD4_45_54];
    AgePopsTotal[s]:=max(AgePopsTotal[s],1); // to avoid division by 0

    (* On ART
       get total in 4 age groups *)
//    AgePopsART[s,HV_DP_CD4_15_24]:=HV_GetDPOnART(p, t, s, HV_DP_A15_19)+
//                                   HV_GetDPOnART(p, t, s, HV_DP_A20_24);
//    AgePopsART[s,HV_DP_CD4_25_34]:=HV_GetDPOnART(p, t, s, HV_DP_A25_29)+
//                                   HV_GetDPOnART(p, t, s, HV_DP_A30_34);
//    AgePopsART[s,HV_DP_CD4_35_44]:=HV_GetDPOnART(p, t, s, HV_DP_A35_39)+
//                                   HV_GetDPOnART(p, t, s, HV_DP_A40_44);
//    AgePopsART[s,HV_DP_CD4_45_54]:=HV_GetDPOnART(p, t, s, HV_DP_A45_49)+
//                                   HV_GetDPOnART(p, t, s, HV_DP_A50_54);

    AgePopsART[s,HV_DP_CD4_15_24]:=HV_GetDPPopSexAgeHIV(p, t, s, HV_DP_A15_19, DP_Pop_H_OnART)+
                                HV_GetDPPopSexAgeHIV(p, t, s, HV_DP_A20_24, DP_Pop_H_OnART);
    AgePopsART[s,HV_DP_CD4_25_34]:=HV_GetDPPopSexAgeHIV(p, t, s, HV_DP_A25_29, DP_Pop_H_OnART)+
                                HV_GetDPPopSexAgeHIV(p, t, s, HV_DP_A30_34, DP_Pop_H_OnART);
    AgePopsART[s,HV_DP_CD4_35_44]:=HV_GetDPPopSexAgeHIV(p, t, s, HV_DP_A35_39, DP_Pop_H_OnART)+
                                HV_GetDPPopSexAgeHIV(p, t, s, HV_DP_A40_44, DP_Pop_H_OnART);
    AgePopsART[s,HV_DP_CD4_45_54]:=HV_GetDPPopSexAgeHIV(p, t, s, HV_DP_A45_49, DP_Pop_H_OnART)+
                                HV_GetDPPopSexAgeHIV(p, t, s, HV_DP_A50_54, DP_Pop_H_OnART);



    (* overall total *)
    AgePopsTotalART[s]:= AgePopsART[s,HV_DP_CD4_15_24]+AgePopsART[s,HV_DP_CD4_25_34]+
                         AgePopsART[s,HV_DP_CD4_35_44]+AgePopsART[s,HV_DP_CD4_45_54];
    AgePopsTotalART[s]:=max(AgePopsTotalART[s],1);
  end;{s loop}

  (* Mortality by CD4 averaged over age groups *)
  for s := HV_Male to HV_Female do
  begin
    for h := HV_CD4_GT500 to HV_CD4_LT50 do
    begin
      _mu[h,s]:=HV_GetDPMortByCD4NoART(p,s,HV_DP_CD4_25_34,h);
      if(AgePopsTotal[s]>0) then
      begin
        _mu[h,s]:=0;
        for a := HV_DP_CD4_15_24 to HV_DP_CD4_45_54 do
        begin
          _mu[h,s]:=_mu[h,s]+
          (AgePops[s,a]/AgePopsTotal[s])*HV_GetDPMortByCD4NoART(p,s,a,h);
        end;{a loop}
      end;{if(AgePopsTotal[s]>0)}

      if t>1 then
       _mu[h, s] := max( _mu[h, s], _NonAIDSDeathRate[s,t-1]);

    end;{h loop}
     _mu[HV_Primary,s]:=_mu[HV_CD4_GT500,s];
  end;{s loop}

  (* CD4 transition by CD4 averaged over age groups *)
  for s := HV_Male to HV_Female do
  begin
    for h := HV_CD4_GT500 to HV_CD4_50_99 do
    begin
      _lambda[h,s]:= HV_GetDPAnnRateProgressLowerCD4(p,s,HV_DP_CD4_25_34,h);
      if(AgePopsTotal[s]>0) then
      begin
        _lambda[h,s]:=0;
        for a := HV_DP_CD4_15_24 to HV_DP_CD4_45_54 do
        begin
          if(h=HV_CD4_GT500) then
          begin
            stageDuration := 1.0 / HV_GetDPAnnRateProgressLowerCD4(p,s,a,h) - _YearsInPrimaryStage;
            _lambda[h,s]:=_lambda[h,s]+
            (AgePops[s,a]/AgePopsTotal[s])*(1-exp(-1*(1.0/stageDuration)/HV_TimeStep))
          end
          else
            _lambda[h,s]:=_lambda[h,s]+
            (AgePops[s,a]/AgePopsTotal[s])*(1-exp(-1*HV_GetDPAnnRateProgressLowerCD4(p,s,a,h)/HV_TimeStep));
        end;{a loop}
         _lambda[h,s]:=_lambda[h,s]*HV_TimeStep;
      end;{if(AgePopsTotal[s]>0)}
    end;{h loop}
    _lambda[HV_Primary,s]:=1/_YearsInPrimaryStage;{this transition is hanlded by delay term}
    _lambda[HV_CD4_LT50,s]:=0;{no transition}
  end;{s loop}


  (* Mortality not on ART by CD4 averaged over age groups *)
  for s := HV_Male to HV_Female do
  begin
    (* note the mapping between h and DP_CD4_GT500,... *)
    for h := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
    begin
      _alpha[h,s]:=HV_GetDPMortByCD4WithARTGT12(p,s,HV_DP_CD4_25_34,h-10);
      if(AgePopsTotalART[s]>0) then
      begin
        _alpha[h,s]:=0;
        for a := HV_DP_CD4_15_24 to HV_DP_CD4_45_54 do
        begin
          _alpha[h,s]:=_alpha[h,s]+(AgePopsART[s,a]/AgePopsTotalART[s]) * HV_GetDPMortByCD4WithARTGT12(p,s,a,h-10)
        end;{a loop}
      end;//if(AgePopsTotalart[s]>0)

      if t>1 then
       _alpha[h, s] := max (_alpha[h, s],  _NonAIDSDeathRate[s,t-1]);

    { Reduce mortality by ART by 50% as viral suppression increases to 95% }
    { Input is percent not virally suppressed, initial default is 0.25 }
    { Fast-Track target is 0.05 }
    _alpha[h,s] := _alpha[h,s]
                   * (1 - 0.5
                 * (HV.GetInfectMultiplierOnART(1) - HV.GetInfectMultiplierOnART(t))
                 / (HV.GetInfectMultiplierOnART(1) - 0.05));


    end;{h loop}
  end;{s loop}
end;

procedure THVProj.Calc_SetARTEligibility(t: integer);
var
  h : byte;
begin
  _h_set := 0;
  _h_cov := 0;
  _CD4ThreshHoldAdults:=HV_GetDPCD4ThreshHoldAdults(p,t);
  for h := HV_CD4_GT500 to HV_CD4_LT50 do
  begin
    if( (_CD4_lower_limit[h]<_CD4ThreshHoldAdults)and(_h_set=0) ) then
    begin
      _h_cov := h;
      _h_set := 1;
    end;{if}
  end;{h loop}
end;

procedure THVProj.Calc_AllocateART(t: integer);
(* Calculate ART enrollment rate based on achieving a specified coverage or number receiving ART
  CDP note ART for HV_None *)
var
  s,r,h,h1,v,sp,s1,r1,a           : byte;
  AllocIter, AllocSafe, AllocMaxIter : integer;
  CurARTCoverage,
  EligibleForARTNum,
  ReceivingARTNum,
  ToReceiveARTNum,
  NeedART,
//  NewARTReceivers,
  ARTMortSum,
  ARTMortSum1,
  AllocatedART,
  sumAdult,
  NewARTReceiversNotKP        : double;
  EligibleForARTStage,
  EligibleForARTPerc,
  MortRate,
  Prop1,
  Prop2,
  NewPatients         : array[HV_Male..HV_Female, HV_None..HV_MSMIDU, HV_CD4_GT500..HV_CD4_LT50] of double;
  NewARTEntrants      : array[HV_Male..HV_Female, HV_None..HV_MSMIDU, HV_RN_UnV..HV_RN_NoProt, HV_CD4_GT500..HV_CD4_LT50] of double;
  TotalNewARTEntrants,
  Entrants, dummy,
  WomenE, TotalWomen  :double;

  KPSet,
  NotKPSet : array [GB_Male..GB_Female] of GB_TByteSet;

  //variables are for John Hopkins University PPL.
  NewARTReceivers,
  receiving, eligible, currCoverage, newCoverage, NewARTPatientsRG : array[HV_Male..HV_Female, HV_None..HV_MSMIDU] of double;
  weight      : array[HV_Male..HV_Female, HV_None..HV_MSMIDU, HV_RN_UnV..HV_RN_NoProt, HV_CD4_GT500..HV_CD4_LT50] of double;

  UnmetNeed : array [HV_None..HV_MSMIDU] of double;
  UnmetNeedNotKP : double;

  ReceivingARTNumByRiskGrp,
  EligibleForARTByRiskGrp,
  ToReceiveARTNumByRiskGrp,
  EligibleForARTNumByRiskGrp      : array[GB_Male..GB_Female, HV_None..HV_MSMIDU] of double;
begin

  (* apply one CD4 coverage for all CD4 categories *)
  if(HV_GetDPARTCoverageSelection(p) = HV_DP_NumOrPercent)
  or((HV_GetDPARTCoverageSelection(p) = HV_DP_PcntByRiskGrp) and (t<=HV.GetGoalsBaseYearIdx)) then
  begin
    (* Number receiving ART by CD4 ranges *)
    ReceivingARTNum := 0;
    (* Number eligible for ART among those in the currently eligible CD4 ranges *)
    EligibleForARTNum := 0;
    for h := HV_CD4_GT500 to HV_CD4_LT50 do
      for s := HV_Male to HV_Female do
        for r := HV_None to HV_MSMIDU do
          begin
            EligibleForARTStage[s,r,h]:=0;
            EligibleForARTPerc[s,r,h]:=0;
            NewPatients[s,r,h]:=0;
            Prop1[s,r,h]:=0;
            Prop2[s,r,h]:=0;
          end;

    TotalWomen := 0;
    if (GetDPData(p).getPopsEligTreat(DP_EligTreatPregnantWomen).Eligible)
      and (GetDPData(p).getPopsEligTreat(DP_EligTreatPregnantWomen).Year <= (GetDPData(p).getFirstYear + t - 1))  then
      for h := HV_CD4_GT500 to HV_CD4_LT50 do
          for r := HV_None to HV_MSMIDU do
            for v := HV_RN_UnV to HV_RN_NoProt do
               TotalWomen :=   TotalWomen + _Adults[HV_Female,r,h,v,t] + _Adults[HV_Female,r,h+10,v,t] ;
    (* Number to receive ART among those in the currently eligible CD4 ranges (i.e. Infected + ART receivers) *)

    for h := HV_CD4_GT500 to HV_CD4_LT50 do
    begin
      for s := HV_Male to HV_Female do
        for r := HV_None to HV_MSMIDU do
          for v := HV_RN_UnV to HV_RN_NoProt do
          begin
            ReceivingARTNum := ReceivingARTNum + _Adults[s,r,h+10,v,t];
//            EPW Following line may be useful
//            ReceivingART[s,r]:= ReceivingART[s,r] + _Adults[s,r,h+10,v,t];

            if _CD4_lower_limit[h] < _CD4ThreshHoldAdults then
            begin
              EligibleForARTNum := EligibleForARTNum + _Adults[s,r,h,v,t];
              EligibleForARTStage[s,r,h]:= EligibleForARTStage[s,r,h] + _Adults[s,r,h,v,t];
            end
            else
            (*Begin: Special pop calculations *)
            begin
              dummy := 0;
              WomenE := 0;
              for sp:= DP_EligTreatTB_HIV to DP_EligTreatOtherPop do
                if (sp = DP_EligTreatTB_HIV) or (sp = DP_EligTreatDiscordantCouples) or  (sp = DP_EligTreatOtherPop) then
                   if (GetDPData(p).getPopsEligTreat(sp).Eligible)
                    and (GetDPData(p).getPopsEligTreat(sp).Year <= (GetDPData(p).getFirstYear + t - 1)) then
                      dummy := dummy + GetDPData(p).getPopsEligTreat(sp).PercentHIV ;


             (*Add sex workers if eligible*)
              if (r = HV_HRH) and (s = HV_Female) and (GetDPData(p).getPopsEligTreat(DP_EligTreatSexWorkers).Eligible)
               and (GetDPData(p).getPopsEligTreat(DP_EligTreatSexWorkers).Year <= (GetDPData(p).getFirstYear + t - 1)) then
                begin
                 EligibleForARTNum := EligibleForARTNum + _Adults[s,r,h,v,t];
                 EligibleForARTStage[s,r,h]:= EligibleForARTStage[s,r,h] + _Adults[s,r,h,v,t];
                 WomenE :=  WomenE + _Adults[s,r,h,v,t];
                end;

            (*Add IDU if eligibile*)
              if ((r = HV_IDU)) and (GetDPData(p).getPopsEligTreat(DP_EligTreatIDU).Eligible)
               and (GetDPData(p).getPopsEligTreat(DP_EligTreatIDU).Year <= (GetDPData(p).getFirstYear + t - 1)) then
               begin
                EligibleForARTNum := EligibleForARTNum + _Adults[s,r,h,v,t];
                EligibleForARTStage[s,r,h]:= EligibleForARTStage[s,r,h] + _Adults[s,r,h,v,t];
                WomenE :=  WomenE + _Adults[s,r,h,v,t];
               end ;


            (*Add MSM if eligibile*)
              if (r = HV_MSM) and (s = HV_Male) and (GetDPData(p).getPopsEligTreat(DP_EligTreatMSM).Eligible)
               and (GetDPData(p).getPopsEligTreat(DP_EligTreatMSM).Year <= (GetDPData(p).getFirstYear + t - 1)) then
               begin
                EligibleForARTNum := EligibleForARTNum + _Adults[s,r,h,v,t];
                EligibleForARTStage[s,r,h]:= EligibleForARTStage[s,r,h] + _Adults[s,r,h,v,t];
               end;

            (*Add SDC + TB + Other populations*)
             if dummy >0 then
               begin
                EligibleForARTNum := EligibleForARTNum + _Adults[s,r,h,v,t] * dummy/100;
                EligibleForARTStage[s,r,h]:= EligibleForARTStage[s,r,h] + _Adults[s,r,h,v,t] * dummy/100;
                WomenE :=  WomenE + _Adults[s,r,h,v,t] * dummy/100;
               end;

            (*Add Pregnant WOmen if eligible*)
            if (s = HV_FeMale) and (GetDPData(p).getPopsEligTreat(DP_EligTreatPregnantWomen).Eligible)
               and (GetDPData(p).getPopsEligTreat(DP_EligTreatPregnantWomen).Year <= (GetDPData(p).getFirstYear + t - 1)) then
               begin
                EligibleForARTNum := EligibleForARTNum
                                     + max(0, (GetDPData(p).getNeedPMTCT(t) * (_Adults[HV_Female,r,h,v,t] - WomenE)/ TotalWomen));
                EligibleForARTStage[s,r,h]:= EligibleForARTStage[s,r,h]
                                            + max(0, (GetDPData(p).getNeedPMTCT(t) * (_Adults[HV_Female,r,h,v,t] - WomenE)/ TotalWomen));
                (*Alternate:
                If = year of eligibility  then
                EligibleForARTNum := EligibleForARTNum + GetDPNeedPMTCT(p,t) * (Adults[s,r,h,v,t]-WomenE)   /TotalWomen
                else
                EligibleForARTNum := EligibleForARTNum + GetDPNeedPMTCT(p,t) * PropFirstPreg * (Adults[s,r,h,v,t]-WomenE) / totalWomen;

                 *)
               end;
            end;
           (*END: Special pop calculations*)
          end;{for v,r,s}
    end;{for h}

   (*begin: Special pop calculations  *)
    for h := HV_CD4_GT500 to HV_CD4_LT50 do
      for s := HV_Male to HV_Female do
        for r := HV_None to HV_MSMIDU do
         begin
          sumAdult:= 0;
          for v := HV_RN_UnV to HV_RN_NoProt do
            sumAdult := sumAdult + _Adults[s,r,h,v,t];
          if sumAdult > 0 then
            EligibleForARTPerc[s,r,h] := EligibleForARTStage[s,r,h] / sumAdult;

//            _Adults[s,r,h,v,t] * EligibleForARTPerc[s,r,h]/
         end;
       (*END: Special pop calculations *)

    (* Need for ART: sum of those on ART and those not on ART but eligible for ART *)
    //NeedART:=EligibleForARTNum+ReceivingARTNum;

    //old way :
    //ToReceiveARTNum := NeedART*ARTCoverage;
    //new way: EPW New values need to come by risk group/sex
    ToReceiveARTNum := GetDPData(p).getPop15to49ART(t);
                    //   HV.GetARTCoverageByRG(s, rg, t)
    //alternative way:
    //ToReceiveARTNum := HV_DPGetNumberOnART(p,t,HV_Bothsexes,HV_DP_A45_49);

      (* New ART receivers *)
      for s := HV_Male to HV_Female do
        for r := HV_None to HV_MSMIDU do
        begin
          NewARTReceivers[s,r] := ToReceiveARTNum - ReceivingARTNum;
           (* can't initiate more than eligible *)
          NewARTReceivers[s,r]:=min(NewARTReceivers[s, r],0.95*EligibleForARTNum);//[CDP note stability issue]
          //EPW come back here to change 0.95 to 0.99?
           (* can't remove ART patients! *)
          NewARTReceivers[s,r]:=max(NewARTReceivers[s, r],0);
        end;

    (* Calculate Option 1 : allocate according to number eligible
       Same proportion starting treatment in each CD4 category *)
    for h := HV_CD4_GT500 to HV_CD4_LT50 do
      for s := HV_Male to HV_Female do
        for r := HV_None to HV_MSMIDU do
          if EligibleForARTNum = 0
            then Prop1[s,r,h] := 0
          else
            Prop1[s,r,h] := NewARTReceivers[s, r] / EligibleForARTNum;

    (* Calculate Option 2 : allocation is weighted by mortality
       First, calculate mortality among those eligible for treatment *)
    for h := HV_CD4_GT500 to HV_CD4_LT50 do
      for s := HV_Male to HV_Female do
        for r := HV_None to HV_MSMIDU do
        begin
          ARTMortSum := _mu[h,s] * EligibleForARTStage[s,r,h];
          ARTMortSum1:=0;
          for v := HV_RN_UnV to HV_RN_NoProt do
          begin
            //ARTMortSum1:= ARTMortSum1 + _Adults[s,r,h,v,t] ;
            ARTMortSum1:= ARTMortSum1 + _Adults[s,r,h,v,t] *  EligibleForARTPerc[s,r,h];
          end;{for v}

          if ARTMortSum1 = 0 then
            MortRate[s,r,h] := 0
          else
            MortRate[s,r,h] := ARTMortSum / ARTMortSum1;
        end;{for s,r,h}

    {* Second, calculate weighted average of number eligible and mortality *}
    {* Allocates new ART recipients proportional to expected mortality in each
       eligible population. If the ART allocation exceeds the need in some
       populations, we loop to distribute the excess among the others. *}
    ARTMortSum := 0;
    AllocatedART := 0;
    for h := HV_CD4_GT500 to HV_CD4_LT50 do
      for s := HV_Male to HV_Female do
        for r := HV_None to HV_MSMIDU do
        begin
          NewPatients[s,r,h] := 0;
          ARTMortSum := ARTMortSum + EligibleForARTStage[s,r,h] * MortRate[s,r,h];
        end;
    if ARTMortSum > 0 then
    begin
      {* AllocMaxIter is the maximum number of passes required to allocate ART.
         This is based on a worst-case scenario in which one group is saturated
         at each pass. In general we will need far fewer passes *}
      AllocMaxIter := (HV_CD4_LT50 - HV_CD4_GT500 + 1) * (HV_Female - HV_Male + 1) * (HV_MSMIDU - HV_None + 1);
      for AllocIter := 1 to AllocMaxIter do
      begin
        AllocSafe := 0;
        for h := HV_CD4_GT500 to HV_CD4_LT50 do
          for s := HV_Male to HV_Female do
            for r := HV_None to HV_MSMIDU do
            begin
              NewPatients[s,r,h] := NewPatients[s,r,h]
                + (NewARTReceivers[s, r] - AllocatedART)
                  * (EligibleForARTStage[s,r,h] - NewPatients[s,r,h]) * MortRate[s,r,h] / ARTMortSum;

              {* Truncate if allocated ART exceeds need, otherwise count the
                 group as safely allocated *}
              if NewPatients[s,r,h] > EligibleForARTStage[s,r,h]
              then NewPatients[s,r,h] := EligibleForARTStage[s,r,h]
              else inc(AllocSafe);
            end;
        ARTMortSum := 0;
        for h := HV_CD4_GT500 to HV_CD4_LT50 do
          for s := HV_Male to HV_Female do
            for r := HV_None to HV_MSMIDU do
            begin
              AllocatedART := AllocatedART + NewPatients[s,r,h];
              ARTMortSum := ARTMortSum + (EligibleForARTStage[s,r,h] - NewPatients[s,r,h]) * MortRate[s,r,h];
            end;
        {* Exit the loop early if all groups have been safely allocated *}
        if AllocSafe = AllocMaxIter then break;
      end;
    end; {if(HV_GetDPARTCoverageSelection(p) = HV_NumOrPercent)}

    (* Third, calculate proportion starting *)
    for h := HV_CD4_GT500 to HV_CD4_LT50 do
     for s := HV_Male to HV_Female do
        for r := HV_None to HV_MSMIDU do
          if (EligibleForARTStage[s,r,h] > 0) then
          begin
            if (EligibleForARTStage[s,r,h] <= 0) then
              Prop2[s,r,h] := 0
            else
              Prop2[s,r,h] := NewPatients[s,r,h] / EligibleForARTStage[s,r,h];
          end;

      (* Finally, average methods one and two *)
      for h := HV_CD4_GT500 to HV_CD4_LT50 do
       for s := HV_Male to HV_Female do
        for r := HV_None to HV_MSMIDU do
         if (EligibleForARTStage[s,r,h] > 0) then
          NewPatients[s,r,h] := EligibleForARTStage[s,r,h] *
            (  Prop1[s,r,h] * HV_GetDPnewARTAllocElig(p)
             + Prop2[s,r,h] * HV_GetDPNewARTAllocMort(p));
      //[CDP note that ART for each CD4 category can be received from AIM]
//         for h := HV_CD4_GT500 to HV_CD4_LT50 do
//          if (_TN_CD4_lower_limit[h] < _CD4ThreshHoldAdults) then
//              NewPatients1[h] :=GetDPNewARTCD4(p, t, h-1);
//
//          for h := HV_CD4_GT500 to HV_CD4_LT50 do
//          if (_TN_CD4_lower_limit[h] < _CD4ThreshHoldAdults) then
//              NewPatients2[h] :=(Prop1[h] + Prop2[h]) / 2;

  end// if overall CD4 coverage is chosen

  //YutengEPW code Start
  else if (HV_GetDPARTCoverageSelection(p) = HV_DP_PcntByRiskGrp)
  and (t>HV.GetGoalsBaseYearIdx) then
  begin

//    NumOrPercentSet
//    PcntByRiskGrpSet

    EligibleForARTNum := 0;
    for s := HV_Male to HV_Female do
      for r := HV_None to HV_MSMIDU do
      begin
        ReceivingARTNumByRiskGrp[s, r] := 0;
        (* Number eligible for ART among those in the currently eligible CD4 ranges *)
        EligibleForARTByRiskGrp[s, r] := 0;
        EligibleForARTNumByRiskGrp[s, r] := 0;
        for h := HV_CD4_GT500 to HV_CD4_LT50 do
        begin
          EligibleForARTStage[s,r,h]:=0;
          EligibleForARTPerc[s,r,h]:=0;
          NewPatients[s,r,h]:=0;
          Prop1[s,r,h]:=0;
          Prop2[s,r,h]:=0;
        end;
      end;


    (*Allocate MSM, IDU, FSW first via new option*)
    TotalWomen := 0;
    if (GetDPData(p).getPopsEligTreat(DP_EligTreatPregnantWomen).Eligible)
      and (GetDPData(p).getPopsEligTreat(DP_EligTreatPregnantWomen).Year <= (GetDPData(p).getFirstYear + t - 1))  then
      for h := HV_CD4_GT500 to HV_CD4_LT50 do
          for r := HV_None to HV_MSMIDU do //r in riskSet[s]
            for v := HV_RN_UnV to HV_RN_NoProt do
               TotalWomen :=   TotalWomen + _Adults[HV_Female,r,h,v,t] + _Adults[HV_Female,r,h+10,v,t] ;
    (* Number to receive ART among those in the currently eligible CD4 ranges (i.e. Infected + ART receivers) *)

    for h := HV_CD4_GT500 to HV_CD4_LT50 do
    begin
      for s := HV_Male to HV_Female do
        for r := HV_None to HV_MSMIDU do  //r in riskSet[s]
          for v := HV_RN_UnV to HV_RN_NoProt do
          begin
            ReceivingARTNumByRiskGrp[s,r]:= ReceivingARTNumByRiskGrp[s,r] + _Adults[s,r,h+10,v,t];

            if _CD4_lower_limit[h] < _CD4ThreshHoldAdults then
            begin
              EligibleForARTNum := EligibleForARTNum + _Adults[s,r,h,v,t];
              EligibleForARTStage[s,r,h]:= EligibleForARTStage[s,r,h] + _Adults[s,r,h,v,t];
              EligibleForARTNumByRiskGrp[s, r] :=  EligibleForARTNumByRiskGrp[s,r] + _Adults[s,r,h,v,t];
            end
            else
            (*Begin: Special pop calculations *)
            begin
              dummy := 0;
              WomenE := 0;
              for sp:= DP_EligTreatTB_HIV to DP_EligTreatOtherPop do
                if (sp = DP_EligTreatTB_HIV) or (sp = DP_EligTreatDiscordantCouples) or  (sp = DP_EligTreatOtherPop) then
                   if (GetDPData(p).getPopsEligTreat(sp).Eligible)
                    and (GetDPData(p).getPopsEligTreat(sp).Year <= (GetDPData(p).getFirstYear + t - 1)) then
                      dummy := dummy + GetDPData(p).getPopsEligTreat(sp).PercentHIV ;


             (*Add sex workers if eligible*)
              if (r = HV_HRH) and (s = HV_Female) and (GetDPData(p).getPopsEligTreat(DP_EligTreatSexWorkers).Eligible)
               and (GetDPData(p).getPopsEligTreat(DP_EligTreatSexWorkers).Year <= (GetDPData(p).getFirstYear + t - 1)) then
                begin
                 EligibleForARTNum := EligibleForARTNum + _Adults[s,r,h,v,t];
                 EligibleForARTStage[s,r,h]:= EligibleForARTStage[s,r,h] + _Adults[s,r,h,v,t];
                 EligibleForARTNumByRiskGrp[s, r] :=  EligibleForARTNumByRiskGrp[s,r] + _Adults[s,r,h,v,t];
                 WomenE :=  WomenE + _Adults[s,r,h,v,t];
                end;

            (*Add IDU if eligibile*)
              if ((r = HV_IDU)) and (GetDPData(p).getPopsEligTreat(DP_EligTreatIDU).Eligible)
               and (GetDPData(p).getPopsEligTreat(DP_EligTreatIDU).Year <= (GetDPData(p).getFirstYear + t - 1)) then
               begin
                EligibleForARTNum := EligibleForARTNum + _Adults[s,r,h,v,t];
                EligibleForARTStage[s,r,h]:= EligibleForARTStage[s,r,h] + _Adults[s,r,h,v,t];
                EligibleForARTNumByRiskGrp[s, r] :=  EligibleForARTNumByRiskGrp[s,r] + _Adults[s,r,h,v,t];
                WomenE :=  WomenE + _Adults[s,r,h,v,t];
               end ;


            (*Add MSM if eligibile*)
              if (r = HV_MSM) and (s = HV_Male) and (GetDPData(p).getPopsEligTreat(DP_EligTreatMSM).Eligible)
               and (GetDPData(p).getPopsEligTreat(DP_EligTreatMSM).Year <= (GetDPData(p).getFirstYear + t - 1)) then
               begin
                EligibleForARTNum := EligibleForARTNum + _Adults[s,r,h,v,t];
                EligibleForARTStage[s,r,h]:= EligibleForARTStage[s,r,h] + _Adults[s,r,h,v,t];
                EligibleForARTNumByRiskGrp[s, r] :=  EligibleForARTNumByRiskGrp[s,r] + _Adults[s,r,h,v,t];
               end;

            (*Add SDC + TB + Other populations*)
             if dummy >0 then
               begin
                EligibleForARTNum := EligibleForARTNum + _Adults[s,r,h,v,t] * dummy/100;
                EligibleForARTStage[s,r,h]:= EligibleForARTStage[s,r,h] + _Adults[s,r,h,v,t] * dummy/100;
                EligibleForARTNumByRiskGrp[s, r] :=  EligibleForARTNumByRiskGrp[s,r] + _Adults[s,r,h,v,t]* dummy/100;
                WomenE :=  WomenE + _Adults[s,r,h,v,t] * dummy/100;
               end;

            (*Add Pregnant WOmen if eligible*)
            if (s = HV_FeMale) and (GetDPData(p).getPopsEligTreat(DP_EligTreatPregnantWomen).Eligible)
               and (GetDPData(p).getPopsEligTreat(DP_EligTreatPregnantWomen).Year <= (GetDPData(p).getFirstYear + t - 1)) then
               begin
                EligibleForARTNum := EligibleForARTNum
                                     + max(0, (GetDPData(p).getNeedPMTCT(t) * (_Adults[HV_Female,r,h,v,t] - WomenE)/ TotalWomen));
                EligibleForARTStage[s,r,h]:= EligibleForARTStage[s,r,h]
                                            + max(0, (GetDPData(p).getNeedPMTCT(t) * (_Adults[HV_Female,r,h,v,t] - WomenE)/ TotalWomen));
                EligibleForARTNumByRiskGrp[s,r]:= EligibleForARTNumByRiskGrp[s,r]
                                            + max(0, (GetDPData(p).getNeedPMTCT(t) * (_Adults[HV_Female,r,h,v,t] - WomenE)/ TotalWomen));
                (*Alternate:
                If = year of eligibility  then
                EligibleForARTNum := EligibleForARTNum + GetDPNeedPMTCT(p,t) * (Adults[s,r,h,v,t]-WomenE)   /TotalWomen
                else
                EligibleForARTNum := EligibleForARTNum + GetDPNeedPMTCT(p,t) * PropFirstPreg * (Adults[s,r,h,v,t]-WomenE) / totalWomen;

                 *)
               end;
            end;
           (*END: Special pop calculations*)
          end;{for v,r,s}
    end;{for h}

   (*begin: Special pop calculations  *)
    for h := HV_CD4_GT500 to HV_CD4_LT50 do
      for s := HV_Male to HV_Female do
        for r := HV_None to HV_MSMIDU do
         begin
           sumAdult:= 0;
           for v := HV_RN_UnV to HV_RN_NoProt do
             sumAdult := sumAdult + _Adults[s,r,h,v,t];
           if sumAdult > 0 then
             EligibleForARTPerc[s,r,h] := EligibleForARTStage[s,r,h] / sumAdult;

//            _Adults[s,r,h,v,t] * EligibleForARTPerc[s,r,h]/
         end;
       (*END: Special pop calculations *)

    (* Need for ART: sum of those on ART and those not on ART but eligible for ART *)
    //NeedART:=EligibleForARTNum+ReceivingARTNum;

    for s := HV_Male to HV_Female do
      for r := HV_None to HV_MSMIDU do
         NewARTReceivers[s, r] := 0;



   KPSet[GB_Male]   := [HV_IDU, HV_MSM];
   KPSet[GB_Female] := [HV_HRH, HV_IDU];

   notKPSet[GB_Male]   := [HV_None .. HV_MSMIDU] - KPSet[GB_Male];
   notKPSet[GB_Female] := [HV_None .. HV_IDU] - KPSet[GB_Female];

    for s := HV_Male to HV_Female do
    begin
      for r in KPSet[s] do //r in riskSet[s]
//      for r := HV_None to HV_MSMIDU do //r in riskSet[s]
      begin

        //old way :
        //ToReceiveARTNum := NeedART*ARTCoverage;
        //new way: EPW New values need to come by risk group/sex
  //      ToReceiveARTNumByRiskGrp[s, r] := GetDPData(p).getPop15to49ART(t);  old


        sumAdult := 0;
        for h in [HV_Primary..HV_CD4_LT50, HV_CD4_GT500_ART..HV_CD4_LT50_ART] do
          for v := HV_RN_UnV to HV_RN_NoProt do
            sumAdult := sumAdult+ _Adults[s, r, h, v, t];


        ToReceiveARTNumByRiskGrp[s, r] := sumAdult*HV.GetARTInputCoverageByRG(s, r, t)/100;

        //alternative way:
        //ToReceiveARTNum := HV_DPGetNumberOnART(p,t,HV_Bothsexes,HV_DP_A45_49);

          (* New ART receivers *)
         NewARTReceivers[s, r] := ToReceiveARTNumByRiskGrp[s, r] - ReceivingARTNumByRiskGrp[s, r];
          (* can't initiate more than eligible, note: code will only work for
             universal eligibility *)

         NewARTReceivers[s, r]:=min(NewARTReceivers[s, r],0.999*EligibleForARTNumByRiskGrp[s,r]);//[CDP note stability issue]
          (* can't remove ART patients! *)
         NewARTReceivers[s, r]:=max(NewARTReceivers[s, r],0);
      end;

      //initial value is the to receive ART total
      NewARTReceiversNotKP := GetNumberOnART(GetDPData(p), t, s, DP_A45_49);

      //procduces to receive not KP
      for r in KPSet[s] do
        NewARTReceiversNotKP := NewARTReceiversNotKP - ToReceiveARTNumByRiskGrp[s, r];

      //produces the new art receivers not KP
      for r in notKPSet[s] do
        NewARTReceiversNotKP := NewARTReceiversNotKP - ReceivingARTNumByRiskGrp[s,r];


       UnmetNeedNotKP := 0;
       for r in notKPSet[s] do // not FSW, MSM, PWID
       begin
         UnmetNeed[r] := 0.0;
         for h in [HV_Primary..HV_CD4_LT50] do// not +10
           for v := HV_RN_UnV to HV_RN_NoProt do
           begin
             UnmetNeed[r] :=  UnmetNeed[r] + _Adults[s, r, h, v, t];
             UnmetNeedNotKP := UnmetNeedNotKP + _Adults[s, r, h, v, t];
           end;
       end;

      for r in notKPSet[s] do //r in riskSet[s]
      begin

        NewARTReceivers[s, r] := NewARTReceiversNotKP*(UnmetNeed[r]/UnmetNeedNotKP);

        NewARTReceivers[s, r]:=min(NewARTReceivers[s, r],0.999*EligibleForARTNumByRiskGrp[s,r]);//[CDP note stability issue]
        (* can't remove ART patients! *)
        NewARTReceivers[s, r]:=max(NewARTReceivers[s, r],0);
      end;

    end;

    for s := HV_Male to HV_Female do
      for r := HV_None to HV_MSMIDU do //r in riskSet[s]
      begin

         if r= HV_MSM then
          NewARTReceivers := NewARTReceivers;

    (* Calculate Option 1 : allocate according to number eligible
       Same proportion starting treatment in each CD4 category *)
         for h := HV_CD4_GT500 to HV_CD4_LT50 do
           if EligibleForARTNumByRiskGrp[s,r] = 0 then
             Prop1[s,r,h] := 0
           else
             Prop1[s,r,h] := NewARTReceivers[s, r] / EligibleForARTNumByRiskGrp[s,r];

    (* Calculate Option 2 : allocation is weighted by mortality
       First, calculate mortality among those eligible for treatment *)
          for h := HV_CD4_GT500 to HV_CD4_LT50 do
          begin
            ARTMortSum := _mu[h,s] * EligibleForARTStage[s,r,h];
            ARTMortSum1 := 0;
            for v := HV_RN_UnV to HV_RN_NoProt do
            begin
              //ARTMortSum1:= ARTMortSum1 + _Adults[s,r,h,v,t] ;
              ARTMortSum1:= ARTMortSum1 + _Adults[s,r,h,v,t] *  EligibleForARTPerc[s,r,h];
            end;{for v}

            if ARTMortSum1 = 0 then
              MortRate[s,r,h] := 0
            else
              MortRate[s,r,h] := ARTMortSum / ARTMortSum1;
          end;{for s,r,h}

    {* Second, calculate weighted average of number eligible and mortality *}
    {* Allocates new ART recipients proportional to expected mortality in each
       eligible population. If the ART allocation exceeds the need in some
       populations, we loop to distribute the excess among the others. *}
          ARTMortSum := 0;
          AllocatedART := 0;
          for h := HV_CD4_GT500 to HV_CD4_LT50 do
          begin
            NewPatients[s,r,h] := 0;
            ARTMortSum := ARTMortSum + EligibleForARTStage[s,r,h] * MortRate[s,r,h];
          end;
          if ARTMortSum > 0 then
          begin
            {* AllocMaxIter is the maximum number of passes required to allocate ART.
               This is based on a worst-case scenario in which one group is saturated
               at each pass. In general we will need far fewer passes *}
            AllocMaxIter := (HV_CD4_LT50 - HV_CD4_GT500 + 1) ;
            for AllocIter := 1 to AllocMaxIter do
            begin
              AllocSafe := 0;
              for h := HV_CD4_GT500 to HV_CD4_LT50 do
              begin
                NewPatients[s,r,h] := NewPatients[s,r,h]
                  + (NewARTReceivers[s, r] - AllocatedART)
                    * (EligibleForARTStage[s,r,h] - NewPatients[s,r,h]) * MortRate[s,r,h] / ARTMortSum;

                {* Truncate if allocated ART exceeds need, otherwise count the
                   group as safely allocated *}
                if NewPatients[s,r,h] > EligibleForARTStage[s,r,h]
                then NewPatients[s,r,h] := EligibleForARTStage[s,r,h]
                else inc(AllocSafe);
              end;

              ARTMortSum := 0;
              for h := HV_CD4_GT500 to HV_CD4_LT50 do
              begin
                AllocatedART := AllocatedART + NewPatients[s,r,h];
                ARTMortSum := ARTMortSum + (EligibleForARTStage[s,r,h] - NewPatients[s,r,h]) * MortRate[s,r,h];
              end;
              {* Exit the loop early if all groups have been safely allocated *}
              if (AllocSafe = AllocMaxIter) or (ARTMortSum=0) then break;
            end;
          end; {if(HV_GetDPARTCoverageSelection(p) = HV_DP_NumOrPercent)}

          (* Third, calculate proportion starting *)
          for h := HV_CD4_GT500 to HV_CD4_LT50 do
            if (EligibleForARTStage[s,r,h] > 0) then
            begin
              if (EligibleForARTStage[s,r,h] <= 0) then
                Prop2[s,r,h] := 0
              else
                Prop2[s,r,h] := NewPatients[s,r,h] / EligibleForARTStage[s,r,h];
            end;

            (* Finally, average methods one and two *)
          for h := HV_CD4_GT500 to HV_CD4_LT50 do
             if (EligibleForARTStage[s,r,h] > 0) then
              NewPatients[s,r,h] := EligibleForARTStage[s,r,h] *
                (  Prop1[s,r,h] * HV_GetDPnewARTAllocElig(p)
                 + Prop2[s,r,h] * HV_GetDPNewARTAllocMort(p));

      end;//Sex/Risk for loop


    (*Allocate rest via Number or percent algorithm*)

  end;//end of else percent option

  (* ART entrants are proportionally allocated to each CD4 bin so as to achieve
  a specified number of HIV+ cases receiving ART *)
  for h := HV_CD4_GT500 to HV_CD4_LT50 do
    for s := HV_Male to HV_Female do
      for r := HV_None to HV_MSMIDU do
        for v := HV_RN_UnV to HV_RN_NoProt do
          NewARTEntrants[s,r,v,h] := 0;

  TotalNewARTEntrants:=0;
  ReceivingARTNum:=0;

  for h := HV_CD4_GT500 to HV_CD4_LT50 do
  begin
    for s := HV_Male to HV_Female do
      for r := HV_None to HV_MSMIDU do
        for v := HV_RN_UnV to HV_RN_NoProt do
        begin
            NewARTEntrants[s,r,v,h] := 0;

            if( (HV_GetDPARTCoverageSelection(p) in [HV_DP_NumOrPercent, HV_DP_PcntByRiskGrp])
            and (EligibleForARTStage[s,r,h]>0) ) then
               NewARTEntrants[s,r,v,h] := (_Adults[s,r,h,v,t] * EligibleForARTPerc[s,r,h]/EligibleForARTStage[s,r,h]) * NewPatients[s,r,h]
                             (* if allocating by CD4 coverage then calc is much easier... *)
            else if (HV_GetDPARTCoverageSelection(p) = HV_DP_CD4Percent) then
            begin

              EligibleForARTNum:=_Adults[s,r,h,v,t];//eligible
              NeedART:=EligibleForARTNum+_Adults[s,r,h+10,v,t];

              CurARTCoverage:=0;
              if(NeedART>0) then
                CurARTCoverage:=_Adults[s,r,h+10,v,t]/NeedART;

              (*Changing coverage by risk group plus CD4 count *)

              if _HV_RiskGrpARTPerc[s,r,h,v,t] > CurARTCoverage then
                NewARTEntrants[s,r,v,h]:=(_HV_RiskGrpARTPerc[s,r,h,v,t]-CurARTCoverage)*NeedART;

              // if(HV_GetDPCD4Coverage(p,t,h)/100>CurARTCoverage) then
               // NewARTEntrants[s,r,v,h]:=(HV_GetDPCD4Coverage(p,t,h)/100-CurARTCoverage)*NeedART;

              (* can't initiate more than eligible *)
              NewARTEntrants[s,r,v,h]:=min(NewARTEntrants[s,r,v,h],EligibleForARTNum);
              (* can't remove ART patients! *)
              NewARTEntrants[s,r,v,h]:=max(NewARTEntrants[s,r,v,h],0);
            end;{else if(HV_GetDPARTCoverageSelection(p)= HV_DP_CD4Percent)}

            TotalNewARTEntrants:=TotalNewARTEntrants+NewARTEntrants[s,r,v,h];
            ReceivingARTNum:=ReceivingARTNum+_Adults[s,r,h+10,v,t];
        end;{for v,r,s}
  end;{for h}


//  if  (HV_GetDPARTCoverageSelection(p)= HV_ARTPercRiskGrp)
//  and (t>HV.GetGoalsBaseYearIdx) then
//  begin
//
//    for h := HV_CD4_GT500 to HV_CD4_LT50 do
//    begin
//      for s := HV_Male to HV_Female do
//        for r := HV_None to HV_MSMIDU do
//          for v := HV_RN_UnV to HV_RN_NoProt do
//          begin
//            EligibleForARTNum := _Adults[s,r,h,v,t];//eligible
//            NeedART := EligibleForARTNum+_Adults[s,r,h+10,v,t];
//
//            CurARTCoverage := 0;
//            if(NeedART>0) then
//              CurARTCoverage := _Adults[s,r,h+10,v,t]/NeedART;
//
//            (*Changing coverage by risk group plus CD4 count *)
//         // ratio := (prop1[s,r,h] prop2[s,r,h])
//         // newPeopleOnART[s,r,t] =  (HV.GetARTCoverageByRG(s,r,t)- HV.GetCoverageFromXLS(s, r, t))*PopEligable[s,r,t]
//         // NewARTEntrants[s,r,h,t]  = (prop1[s,r,h]+prop2[s,r,h])*newPeopleOnART[s,r,t]
//
//            _HV_RiskGrpARTPerc[s,r,h,v,t] := HV.GetCoverageFromXLS(s, r, t)/100 ;//
//            if _HV_RiskGrpARTPerc[s,r,h,v,t] > CurARTCoverage then
//              NewARTEntrants[s,r,v,h]:=(_HV_RiskGrpARTPerc[s,r,h,v,t]-CurARTCoverage)*NeedART;
//
//            // if(HV_GetDPCD4Coverage(p,t,h)/100>CurARTCoverage) then
//             // NewARTEntrants[s,r,v,h]:=(HV_GetDPCD4Coverage(p,t,h)/100-CurARTCoverage)*NeedART;
//
//            (* can't initiate more than eligible *)
//            NewARTEntrants[s,r,v,h]:=min(NewARTEntrants[s,r,v,h],EligibleForARTNum);
//            (* can't remove ART patients! *)
//            NewARTEntrants[s,r,v,h]:=max(NewARTEntrants[s,r,v,h],0);
//
//            TotalNewARTEntrants:=TotalNewARTEntrants+NewARTEntrants[s,r,v,h];
//            ReceivingARTNum:=ReceivingARTNum+_Adults[s,r,h+10,v,t];
//          end;
//    end;
//  end;{else if(HV_GetDPARTCoverageSelection(p)= HV_ARTPercRiskGrp)}

//  if( (HV_GetDPARTCoverageSelection(p) = HV_NumOrPercent) and (EligibleForARTStage[s,r,h]>0) )
//  or((HV_GetDPARTCoverageSelection(p) = HV_ARTPercRiskGrp)
//     and (EligibleForARTStage[s,r,h]>0)
//     and (t<=HV.GetGoalsBaseYearIdx))then
//	begin
//		for h := HV_CD4_GT500 to HV_CD4_LT50 do
//		begin
//		  for s := HV_Male to HV_Female do
//		    for r := HV_None to HV_MSMIDU do
//          for v := HV_RN_UnV to HV_RN_NoProt do
//          begin
//            NewARTEntrants[s,r,v,h] := (_Adults[s,r,h,v,t] * EligibleForARTPerc[s,r,h]/EligibleForARTStage[s,r,h]) * NewPatients[s,r,h];
//            TotalNewARTEntrants:=TotalNewARTEntrants+NewARTEntrants[s,r,v,h];
//            ReceivingARTNum:=ReceivingARTNum+_Adults[s,r,h+10,v,t];
//          end;
//    end;
//  end
//  (* if allocating by CD4 coverage then calc is much easier... *)
//  else if(HV_GetDPARTCoverageSelection(p)= HV_CD4Category) then
//  begin
//
//	  for h := HV_CD4_GT500 to HV_CD4_LT50 do
//	  begin
//		  for s := HV_Male to HV_Female do
//		    for r := HV_None to HV_MSMIDU do
//			    for v := HV_RN_UnV to HV_RN_NoProt do
//          begin
//            EligibleForARTNum:=_Adults[s,r,h,v,t];//eligible
//            NeedART:=EligibleForARTNum+_Adults[s,r,h+10,v,t];
//
//            CurARTCoverage:=0;
//            if(NeedART>0) then
//            CurARTCoverage:=_Adults[s,r,h+10,v,t]/NeedART;
//
//            (*Changing coverage by risk group plus CD4 count *)
//
//            if _HV_RiskGrpARTPerc[s,r,h,v,t] > CurARTCoverage then
//            NewARTEntrants[s,r,v,h]:=(_HV_RiskGrpARTPerc[s,r,h,v,t]-CurARTCoverage)*NeedART;
//
//            // if(HV_GetDPCD4Coverage(p,t,h)/100>CurARTCoverage) then
//             // NewARTEntrants[s,r,v,h]:=(HV_GetDPCD4Coverage(p,t,h)/100-CurARTCoverage)*NeedART;
//
//            (* can't initiate more than eligible *)
//            NewARTEntrants[s,r,v,h]:=min(NewARTEntrants[s,r,v,h],EligibleForARTNum);
//            (* can't remove ART patients! *)
//            NewARTEntrants[s,r,v,h]:=max(NewARTEntrants[s,r,v,h],0);
//
//            TotalNewARTEntrants:=TotalNewARTEntrants+NewARTEntrants[s,r,v,h];
//            ReceivingARTNum:=ReceivingARTNum+_Adults[s,r,h+10,v,t];
//			    end;
//	  end;
//  end{else if(HV_GetDPARTCoverageSelection(p)= HV_CD4Category)}
//  else if(HV_GetDPARTCoverageSelection(p)= HV_ARTPercRiskGrp)
//  and (t>HV.GetGoalsBaseYearIdx) then
//  begin
//	  for h := HV_CD4_GT500 to HV_CD4_LT50 do
//	  begin
//		  for s := HV_Male to HV_Female do
//		    for r := HV_None to HV_MSMIDU do
//			    for v := HV_RN_UnV to HV_RN_NoProt do
//			    begin
//			      EligibleForARTNum := _Adults[s,r,h,v,t];//eligible
//			      NeedART := EligibleForARTNum+_Adults[s,r,h+10,v,t];
//
//			      CurARTCoverage := 0;
//			      if(NeedART>0) then
//			    	CurARTCoverage := _Adults[s,r,h+10,v,t]/NeedART;
//
//			    		  (*Changing coverage by risk group plus CD4 count *)
//		//                  ratio := (prop1[s,r,h] prop2[s,r,h])
//		//                  newPeopleOnART[s,r,t] =  (HV.GetARTCoverageByRG(s,r,t)- HV.GetCoverageFromXLS(s, r, t))*PopEligable[s,r,t]
//		//                  NewARTEntrants[s,r,h,t]  = (prop1[s,r,h]+prop2[s,r,h])*newPeopleOnART[s,r,t]
//
//			    	_HV_RiskGrpARTPerc[s,r,h,v,t] := HV.GetCoverageFromXLS(s, r, t)/100 ;//
//			    	if _HV_RiskGrpARTPerc[s,r,h,v,t] > CurARTCoverage then
//			    	NewARTEntrants[s,r,v,h]:=(_HV_RiskGrpARTPerc[s,r,h,v,t]-CurARTCoverage)*NeedART;
//
//			    	// if(HV_GetDPCD4Coverage(p,t,h)/100>CurARTCoverage) then
//			    	 // NewARTEntrants[s,r,v,h]:=(HV_GetDPCD4Coverage(p,t,h)/100-CurARTCoverage)*NeedART;
//
//			    	(* can't initiate more than eligible *)
//			    	NewARTEntrants[s,r,v,h]:=min(NewARTEntrants[s,r,v,h],EligibleForARTNum);
//			    	(* can't remove ART patients! *)
//			    	NewARTEntrants[s,r,v,h]:=max(NewARTEntrants[s,r,v,h],0);
//
//			    	TotalNewARTEntrants:=TotalNewARTEntrants+NewARTEntrants[s,r,v,h];
//			    	ReceivingARTNum:=ReceivingARTNum+_Adults[s,r,h+10,v,t];
//			    end;
//    end;
//  end;{else if(HV_GetDPARTCoverageSelection(p)= HV_ARTPercRiskGrp)}



 if (HV_GetDPARTCoverageSelection(p)= HV_DP_CD4Percent) then
 begin
  (*Keep targeted risk groups coverage constant and adjust the rest *)
  for h := HV_CD4_LT50 downto HV_CD4_GT500 do
    for s := HV_Male to HV_Female do
      for r := HV_None to HV_MSMIDU do
        for v := HV_RN_UnV to HV_RN_NoProt do
          if (TotalNewARTEntrants > 0) then
          begin
          if ((s = HV_Female) and (r = HV_HRH) and (GetDPData(p).getPopsEligTreat(DP_EligTreatSexWorkers).Eligible))
               or ((s= HV_Male) and (r = HV_MSM)and (GetDPData(p).getPopsEligTreat(DP_EligTreatMSM).Eligible)) or
              ((r = HV_IDU) and (GetDPData(p).getPopsEligTreat(DP_EligTreatIDU).Eligible)) then
          begin
             (* can't remove ART patients! *)
            Entrants := max(NewARTEntrants[s, r, v, h], 0);
            (* can't initiate more than eligible *)
            Entrants := min(_Adults[s, r, h, v, t], Entrants);

            TotalNewARTEntrants := TotalNewARTEntrants - Entrants;
            ReceivingARTNum := ReceivingARTNum  + Entrants;

            _Adults[s, r, h, v, t] := _Adults[s, r, h, v, t] -
              Entrants; { HIV+ not on ART }
            (* note offset on h *)
            _Adults[s, r, h + 10, v, t] := _Adults[s, r, h + 10, v, t] +
              Entrants; { HIV+ on ART }
            (* note offset on h *)
            _NewlyOnART[s, h + 10, t] := _NewlyOnART[s, h + 10, t] +
              Entrants;
           end;
          end; { for h,v,r,s }

   //Adjustment for rest of the population
  (* Renormalize NewARTEntrants to allocate exact number of New ART required to match total ART from AIM *)


  ToReceiveARTNum := GetDPData(p).getPop15to49ART(t) - ReceivingARTNum;
  // ToReceiveARTNum:= TotalNewARTEntrants ;// -ReceivingARTNum;
  (* can't remove ART receivers *)
  ToReceiveARTNum := max(ToReceiveARTNum, 0);

  for h := HV_CD4_LT50 downto HV_CD4_GT500 do
    for s := HV_Male to HV_Female do
      for r := HV_None to HV_MSMIDU do
        for v := HV_RN_UnV to HV_RN_NoProt do
          if (TotalNewARTEntrants > 0) then
          begin
            (*if not(((s = HV_Female) and (r = HV_HRH) )
               or ((s= HV_Male) and (r = HV_MSM)) or
              ((r = HV_IDU) )) then   *)
            if not(((s = HV_Female) and (r = HV_HRH) and (GetDPData(p).getPopsEligTreat(DP_EligTreatSexWorkers).Eligible))
               or ((s= HV_Male) and (r = HV_MSM)and (GetDPData(p).getPopsEligTreat(DP_EligTreatMSM).Eligible)) or
              ((r = HV_IDU) and (GetDPData(p).getPopsEligTreat(DP_EligTreatIDU).Eligible))) then
            begin
            Entrants := (NewARTEntrants[s, r, v, h] / TotalNewARTEntrants)
              * ToReceiveARTNum;
            (* can't remove ART patients! *)
            Entrants := max(Entrants, 0);
            (* can't initiate more than eligible *)
            Entrants := min(_Adults[s, r, h, v, t], Entrants);

            _Adults[s, r, h, v, t] := _Adults[s, r, h, v, t] -
              Entrants; { HIV+ not on ART }
            (* note offset on h *)
            _Adults[s, r, h + 10, v, t] := _Adults[s, r, h + 10, v, t] +
              Entrants; { HIV+ on ART }
            (* note offset on h *)
            _NewlyOnART[s, h + 10, t] := _NewlyOnART[s, h + 10, t] +
              Entrants;
            end;
          end; { for h,v,r,s }
 end
 else if (HV_GetDPARTCoverageSelection(p)= HV_DP_PcntByRiskGrp)
 and (t > HV.GetGoalsBaseYearIdx) then
 begin


  //numReceivingEverybodyElse := GetDPData(p).getPop15to49ART(t) - NumReceiving;
  //for r in FSW, MSM, IDU
  //  numReceivingEverybodyElse -=  NumHIVPos[r]*covByRisk[r]
  //

    (*Keep targeted risk groups coverage constant and adjust the rest *)
    for h := HV_CD4_LT50 downto HV_CD4_GT500 do
      for s := HV_Male to HV_Female do
        for r := HV_None to HV_MSMIDU do
          for v := HV_RN_UnV to HV_RN_NoProt do
            if (TotalNewARTEntrants > 0) then
            begin
               (* can't remove ART patients! *)
              Entrants := max(NewARTEntrants[s, r, v, h], 0);
              (* can't initiate more than eligible *)
              Entrants := min(_Adults[s, r, h, v, t], Entrants);

              TotalNewARTEntrants := TotalNewARTEntrants - Entrants;
              ReceivingARTNum := ReceivingARTNum  + Entrants;

              _Adults[s, r, h, v, t] := _Adults[s, r, h, v, t] -
                Entrants; { HIV+ not on ART }
              (* note offset on h *)
              _Adults[s, r, h + 10, v, t] := _Adults[s, r, h + 10, v, t] +
                Entrants; { HIV+ on ART }
              (* note offset on h *)
              _NewlyOnART[s, h + 10, t] := _NewlyOnART[s, h + 10, t] +
                Entrants;
            end; { for h,v,r,s }

//   //Adjustment for rest of the population
//  (* Renormalize NewARTEntrants to allocate exact number of New ART required to match total ART from AIM *)
//
//  ToReceiveARTNum := GetDPData(p).getPop15to49ART(t) - ReceivingARTNum;
//  // ToReceiveARTNum:= TotalNewARTEntrants ;// -ReceivingARTNum;
//  (* can't remove ART receivers *)
//  ToReceiveARTNum := max(ToReceiveARTNum, 0);
//
//  for h := HV_CD4_LT50 downto HV_CD4_GT500 do
//    for s := HV_Male to HV_Female do
//      for r := HV_None to HV_MSMIDU do
//        for v := HV_RN_UnV to HV_RN_NoProt do
//          if (TotalNewARTEntrants > 0) then
//          begin
//
//            (*if not(((s = HV_Female) and (r = HV_HRH) )
//               or ((s= HV_Male) and (r = HV_MSM)) or
//              ((r = HV_IDU) )) then   *)
//            Entrants := (NewARTEntrants[s, r, v, h] / TotalNewARTEntrants)
//              * ToReceiveARTNum;
//            (* can't remove ART patients! *)
//            Entrants := max(Entrants, 0);
//            (* can't initiate more than eligible *)
//            Entrants := min(_Adults[s, r, h, v, t], Entrants);
//
//            _Adults[s, r, h, v, t] := _Adults[s, r, h, v, t] -
//              Entrants; { HIV+ not on ART }
//            (* note offset on h *)
//            _Adults[s, r, h + 10, v, t] := _Adults[s, r, h + 10, v, t] +
//              Entrants; { HIV+ on ART }
//            (* note offset on h *)
//            _NewlyOnART[s, h + 10, t] := _NewlyOnART[s, h + 10, t] +
//              Entrants;
//          end; { for h,v,r,s }
 end
 else
 begin

   ToReceiveARTNum := GetDPData(p).getPop15to49ART(t) - ReceivingARTNum;
   // ToReceiveARTNum:= TotalNewARTEntrants ;// -ReceivingARTNum;
   (* can't remove ART receivers *)
   ToReceiveARTNum := max(ToReceiveARTNum, 0);

    for h := HV_CD4_LT50 downto HV_CD4_GT500 do
      for s := HV_Male to HV_Female do
        for r := HV_None to HV_MSMIDU do
          for v := HV_RN_UnV to HV_RN_NoProt do
             if(TotalNewARTEntrants>0) then
             begin

                Entrants:=(NewARTEntrants[s,r,v,h]/TotalNewARTEntrants)*ToReceiveARTNum;
                (* can't remove ART patients! *)
                Entrants:=max(entrants,0);
                (* can't initiate more than eligible *)
                Entrants:=min(_Adults[s,r,h,v,t],Entrants);

                _Adults[s,r,h,v,t] := _Adults[s,r,h,v,t] - Entrants;{HIV+ not on ART}
                (* note offset on h *)
                _Adults[s,r,h+10,v,t] := _Adults[s,r,h+10,v,t] + Entrants;{HIV+ on ART}
                (* note offset on h *)
                _NewlyOnART[s,h+10,t]:= _NewlyOnART[s,h+10,t] + Entrants;
             end;{for h,v,r,s}
 end;
end;

procedure THVProj.Init_IntCD4Vars;
var
  t : integer;
  s, h : byte;
begin
  (* In theory the pulse should be distributed by CD4 stage, but doing so
     leads to different dynamics. For the moment all are put in the Primary stage. *)
  { Could average the CD4 distribution according to incidence by age from AIM }
  for s := HV_Male to HV_Female do
  begin
    _NewInfByStage[s,HV_Primary]   := HV_GetDPDistNewInfectionsCD4(p, s, HV_DP_CD4_25_34, HV_CD4_GT500) / 100;
    _NewInfByStage[s,HV_CD4_GT500] := 0;
    for h := HV_CD4_350_500 to HV_CD4_LT50 do
      _NewInfByStage[s,h] := HV_GetDPDistNewInfectionsCD4(p, s, HV_DP_CD4_25_34, h) / 100;
  end;

  for s := HV_Male to HV_Female do
    for h := HV_Primary to HV_CD4_LT50 do
      _PulseByStage[s,h] := _NewInfByStage[s,h];

  (* Infectiousness multipliers *)
  _MultNoART[HV_Primary] := _Infectiousness[HV_PrimaryInf];
  _MultNoART[HV_CD4_GT500] :=  _Infectiousness[HV_Asymptomatic];
  _MultNoART[HV_CD4_350_500] :=  _Infectiousness[HV_Asymptomatic];
  _MultNoART[HV_CD4_250_349] := _Infectiousness[HV_Asymptomatic];
  _MultNoART[HV_CD4_200_249] := _Infectiousness[HV_Asymptomatic];
  _MultNoART[HV_CD4_100_199] := _Infectiousness[HV_SympNoART];
  _MultNoART[HV_CD4_50_99] := _Infectiousness[HV_SympNoART];
  _MultNoART[HV_CD4_LT50] := _Infectiousness[HV_SympNoART];

  (* Model infectiousness on ART as reduction in asymptomatic infectiousness *)
  for h := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
    for t := GetGBCalcYearIdx(p) to GetGBFinalYearIdx(p) do
      _MultART[h,t] := HV.GetInfectMultiplierOnART(t)*_Infectiousness[HV_Asymptomatic];

  {* Lower limits for categories *}
  _CD4_lower_limit[HV_Primary] := 500;
  _CD4_lower_limit[HV_CD4_GT500] := 500;
  _CD4_lower_limit[HV_CD4_350_500] := 350;
  _CD4_lower_limit[HV_CD4_250_349] := 250;
  _CD4_lower_limit[HV_CD4_200_249] := 200;
  _CD4_lower_limit[HV_CD4_100_199] := 100;
  _CD4_lower_limit[HV_CD4_50_99] := 50;
  _CD4_lower_limit[HV_CD4_LT50] := 0;
end;

//falls away
//procedure THVProj.Calc_ARTCoverage(t:integer);
//var
//  s,r,h,v : byte;
//  value   : double;
//begin
//  _ARTReceive[t]:=0;
//
//  if HV_GetDPHAARTPerNum(p,t) = HV_DP_Percent then
//    _ARTReceive[t] := HV_GetDPHAART(p,t)/100
//  else
//  begin
//    value := 0;
//    for s := HV_Male to HV_Female do
//      for r := HV_None to HV_MSMIDU do
//        for h := HV_CD4_GT500 to HV_CD4_LT50 do
//          for v := HV_RN_UnV to HV_RN_NoProt do
//          begin
//            if (_CD4_lower_limit[h] < _CD4ThreshHoldAdults) then
//              value := value + _Adults[s,r,h,v,t];
//              (* note offset on h *)
//            value := value + _Adults[s,r,h+10,v,t];
//          end;
//    if value = 0 then
//      _ARTReceive[t] := 0
//    else
//      _ARTReceive[t] := (HV_GetDPHAART(p,t))/value;
//  end;
//end;

procedure THVProj.Calc_TrackAIDSDeaths(t,_ts:integer);
var
  s,r,h,v : byte;
  value   : double;
begin
  for s := HV_Male to HV_Female do
  begin
    for r := HV_None to HV_MSMIDU do
    begin
      for v := HV_RN_UnV to HV_RN_NoProt do
      begin
       value:=0;
        if not((s = HV_Female) and (r >= HV_MSM)) then
        begin
          for h := HV_Primary to HV_CD4_LT50 do
                 value := value + (1/HV_TimeStep)*_AdultsTS[s,r,h,v,_ts]*_mu[h,s]*_ScaleMuART[h,s];

          for h := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
                 value := value + (1/HV_TimeStep)*_AdultsTS[s,r,h,v,_ts]*_alpha[h,s];
        end;
          //value := (_Adults[s,r,h,v,t]/_YearsInSymptomaticStage) * (1-_ArtSurviving * _ArtReceiveL[h]);
          _AIDSDeaths[s,r,v,t] := _AIDSDeaths[s,r,v,t]+value;
          _TotalAIDSDeaths[t] := _TotalAIDSDeaths[t]+value;
      end;{v loop}
    end;{r loop}
  end;{s loop}

  //AIDSDeathsTotalART:=0;
  for s := HV_Male to HV_Female do
  begin
    for r := HV_None to HV_MSMIDU do
    begin
      value:=0;
      for v := HV_RN_UnV to HV_RN_NoProt do
      begin
        if not((s = HV_Female) and (r >= HV_MSM)) then
        begin
          for h := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
            value := value + (1/HV_TimeStep)*_AdultsTS[s,r,h,v,_ts]*_alpha[h,s];
        end;{if}
      end; {v loop}
      //value := (_Adults[s,r,h,v,t]/_YearsInSymptomaticStage) * (1-_ArtSurviving * _ArtReceiveL[h]);
      _AIDSDeathsART[s,r,t] := _AIDSDeathsART[s,r,t]+value;
      _TotalARTDeaths[t] := _TotalARTDeaths[t]+value;
    end;{r loop}
  end;{s loop}
end;

procedure THVProj.Calc_TrackNewInfections(t:integer);
var
  s,r,v : byte;
begin
  (* accumulate new infection over year *)
  for s := HV_Male to HV_Female do
    for r := HV_LRH to HV_MSMIDU do
      for v := HV_RN_UnV to HV_RN_NoProt  do
         if not((s = HV_Female) and (r >= HV_MSM)) then
         begin
            _NewInfections[s,r,v,t]:=_NewInfections[s,r,v,t]+(1/HV_TimeStep)*_NewInfectionsT[s,r,v];
            _TotalNewInfection[t] := _TotalNewInfection[t] + (1/HV_TimeStep)*_NewInfectionsT[s,r,v];
         end;

  HV.SetTotalNewInfection(t,_TotalNewInfection[t]);

//  for s := HV_Male to HV_Female do
//    for r := HV_LRH to HV_MSMIDU do
//      for v := HV_RN_UnV to HV_RN_NoProt do
//      begin
//        if not((s = HV_Female) and (r >= HV_MSM)) then
//        _TotalNewInfection[t] := _TotalNewInfection[t] + (1/HV_TimeStep)*_NewInfectionsT[s,r,v];
//      end;
//  SetHVTotalNewInfection(p,t,_TotalNewInfection[t]);
end;

procedure THVProj.Calc_InitialPulse(t,_ts:integer);
var
  s,r,h   : byte;
  Pulse   : array[HV_Male..HV_Female,HV_None..HV_MSMIDU] of double;
begin
  for s := HV_Male to HV_Female do
    for r := HV_None to HV_MSMIDU do
      pulse[s,r] := 0;

  (* Calc initial Pulse non IDU *)
  for s := HV_Male to HV_Female do
  begin
    for r := HV_LRH to HV_MSMIDU do
    begin
      (* Use _AdultsTS because _Adults[s,r,HV_Negative,HV_RN_UnV,_ts] is already updated
         making an assumption that 'more' of the initial pulse goes to high risk groups *)
      pulse[s,r] := _InitialPulse*_AdultsTS[s,r,HV_Negative,HV_RN_UnV,_ts];

      if (r=HV_MRH) then
        pulse[s,r] := 3/5*_InitialPulse*_AdultsTS[s,r,HV_Negative,HV_RN_UnV,_ts];

      if (r=HV_LRH) then
        pulse[s,r] := 1/5*_InitialPulse*_AdultsTS[s,r,HV_Negative,HV_RN_UnV,_ts];

      if ((r=HV_None) or (r=HV_IDU)) then
        pulse[s,r] :=0;

      if ((s = HV_Female) and (r=HV_MSM)) then
        pulse[s,r] :=0;
    end;{r loop}
  end;{s loop}

  (* add initial pulse to HIV+ states and remove from HIV- *)
  for s := HV_Male to HV_Female do
    for r := HV_LRH to HV_MSMIDU do
      for h := HV_Primary to HV_CD4_LT50 do
      begin
        _Adults[s,r,h,HV_RN_UnV,t] := _Adults[s,r,h,HV_RN_UnV,t]+pulse[s,r]*_PulseByStage[s,h];
        _Adults[s,r,HV_Negative,HV_RN_UnV,t] := _Adults[s,r,HV_Negative,HV_RN_UnV,t]-pulse[s,r]*_PulseByStage[s,h];
        _NewHIV[s]:=_NewHIV[s]+pulse[s,r]*_PulseByStage[s,h];
      end;
end;

procedure THVProj.Calc_InitialPulseIDU(t,_ts:integer);
var
  s,h   : byte;
  pulse : double;
begin
 (* Calc initial Pulse IDU *)
  for s := HV_Male to HV_Female do
  begin
    (* HIV may be introduced among IDU relative to rest of population
       Add initial pulse to IDU if are no HIV+ IDU already due to sexual transmission *)
    pulse := _InitialPulse*_Adults[s,HV_IDU,HV_Negative,HV_RN_UnV,t];
    for h := HV_Primary to HV_CD4_LT50 do
    begin
      _Adults[s,HV_IDU,h,HV_RN_UnV,t] := _Adults[s,HV_IDU,h,HV_RN_UnV,t]+pulse*_PulseByStage[s,h];
      _Adults[s,HV_IDU,HV_Negative,HV_RN_UnV,t] := _Adults[s,HV_IDU,HV_Negative,HV_RN_UnV,t]-pulse*_PulseByStage[s,h];
      _NewHIV[s]:=_NewHIV[s]+pulse*_PulseByStage[s,h];
    end{h loop}
  end;{s loop}
end;

function THVProj.Calc_NewInfections(t:byte):boolean;
var
  s,r,h,v                             : byte;
  PrevM, PrevF, PrevM_IDU,
  PrevF_IDU, PrevML, PrevFL,
  PrevB, SusceptibleIDU,
  YearsInAsymptomatic,
  TotalSexActsM, TotalSexActsF,
  TotalSexActs, SexActsRatioM,
  SexActsRatioF, SexActsRatioM_IDU,
  SexActsRatioF_IDU,
  rMultNumerator,rMultDenominator,
  rMultNumeratorAll,
  rMultDenominatorAll,
  rMultEqlb, rMultM, rMultF,
  rMultM_IDU, rMultF_IDU,
  _ForceSexTransmIDU,
  _ForceNeedleTransmIDU,
  value                               : double;
begin
  Calc_NewInfections:=False;

  (* HIVInfection *)
  for s := HV_Bothsexes to HV_Female do
    for r := HV_AllRisk to HV_MSMIDU do
      for v := HV_RN_AllVacc to HV_RN_NoProt do
      begin
        _NewInfectionsT[s,r,v]:=0;
      end;{s,r,v loop}

  (* Calculating multiplier for each sex and risk group *)
  rMultEqlb:=1.7;//2.44 in previous version of Goals}
  rMultNumeratorAll:=0;
  rMultDenominatorAll :=0;
  for s := HV_Male to HV_Female do
    for r := HV_LRH to HV_HRH do
    begin
      (* Don't include IDU prevalence in multipliers. IDU force of infection is external input *)
      rMultNumerator:=0;
      rMultDenominator:= 0;
      for h := HV_Primary to HV_CD4_LT50 do
      begin
        rMultNumerator := rMultNumerator + _MultNoART[h]*_Adults[s,r,h,HV_RN_AllVacc,t];
        rMultDenominator := rMultDenominator +_Adults[s,r,h,HV_RN_AllVacc,t];
        rMultNumeratorAll := rMultNumeratorAll + _MultNoART[h]*_Adults[s,r,h,HV_RN_AllVacc,t];
        rMultDenominatorAll := rMultDenominatorAll +_Adults[s,r,h,HV_RN_AllVacc,t];
      end;{h loop}
      for h := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
      begin
        rMultNumerator := rMultNumerator + _MultART[h,t]*_Adults[s,r,h,HV_RN_AllVacc,t];
        rMultDenominator := rMultDenominator + _Adults[s,r,h,HV_RN_AllVacc,t];
        rMultNumeratorAll := rMultNumeratorAll + _MultART[h,t]*_Adults[s,r,h,HV_RN_AllVacc,t];
        rMultDenominatorAll := rMultDenominatorAll + _Adults[s,r,h,HV_RN_AllVacc,t];
      end;{h loop}

      if( (rMultDenominator>0) and (_CalcPrevalence[s,r,t]>0) ) then
        _rMultAll[s,r,t]:= (rMultNumerator / rMultDenominator) / rMultEqlb *
                           ( 1 - (_VaccineEffectiveness[HV_RN_Infectiousness]/100) *
                            _TotalVaccinated[HV_Bothsexes,t-1] /_Populations[HV_Bothsexes,t-1])
      else
        _rMultAll[s,r,t]:=1;
    end;{r,s loop}

  {* The 'global' multiplier *}
  if(rMultDenominatorAll>0) then
    _rMultAll[HV_Bothsexes,HV_AllRisk,t]:= (rMultNumeratorAll / rMultDenominatorAll) / rMultEqlb *
                                           ( 1 - (_VaccineEffectiveness[HV_RN_Infectiousness]/100) *
                                           _TotalVaccinated[HV_Bothsexes,t-1] /_Populations[HV_Bothsexes,t-1])
  else
    _rMultAll[HV_Bothsexes,HV_AllRisk,t]:=1;

  {* The MSM multiplier. MSM not mixing with other risk groups. *}
  rMultNumeratorAll:=0;
  rMultDenominatorAll :=0;
  for r := HV_MSM to HV_MSMHR do
  begin
    for h := HV_Primary to HV_CD4_LT50 do
    begin
      rMultNumeratorAll := rMultNumeratorAll + _MultNoART[h]*_Adults[HV_Male,r,h,HV_RN_AllVacc,t];
      rMultDenominatorAll := rMultDenominatorAll +_Adults[HV_Male,r,h,HV_RN_AllVacc,t];
    end;{h loop}
    for h := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
    begin
      rMultNumeratorAll := rMultNumeratorAll + _MultART[h,t]*_Adults[HV_Male,r,h,HV_RN_AllVacc,t];
      rMultDenominatorAll := rMultDenominatorAll + _Adults[HV_Male,r,h,HV_RN_AllVacc,t];
    end;{h loop}
  end;{r loop}

  for r := HV_MSM to HV_MSMHR do
  begin
    if( (rMultDenominatorAll>0)) then
      _rMultAll[HV_Male,r,t]:= (rMultNumeratorAll / rMultDenominatorAll) / rMultEqlb *
                               ( 1 - (_VaccineEffectiveness[HV_RN_Infectiousness]/100) *
                               _TotalVaccinated[HV_Male,t-1] /_Populations[HV_Male,t-1])
    else
      _rMultAll[HV_Male,r,t]:=1;
  end;{r loop}

  (* New Infections *)

  PrevML := 0;
  PrevFL := 0;
  PrevM_IDU:=0;
  PrevF_IDU:=0;
  rMultM_IDU:=1;
  rMultF_IDU:=1;

  for r := HV_LRH to HV_HRH do
  begin
    (* Balance the number of sex acts between low risk men and women *)

    (* A factor SexActsRatio is applied to the number of sex acts for men
    in the transmission formula. Sometimes the results look more realistic without this factor...
    {True balancing may require a proportional mixing formulation *)
//    SexActsRatioM:= 1;
//    SexActsRatioF:= 1;
    TotalSexActsF:= 0;
    TotalSexActsM:= 0;
//    TotalSexActs:=0;

    if r = HV_LRH then
    begin
      (* Number of low risk SexActs women are having *)
      TotalSexActsF :=_Adults[HV_Female,HV_LRH,HV_AllHIV,HV_RN_AllVacc,t]*_NumPart[HV_LRH_F3,t]*_SexActs[HV_LRH_F3,t]+
      +_Adults[HV_Female,HV_MRH,HV_AllHIV,HV_RN_AllVacc,t]*_PercMarried[HV_MRH_F3]/100*_NumPart[HV_MRH_F3,t]*_SexActs[HV_MRH_F3,t]
      +_Adults[HV_Female,HV_HRH,HV_AllHIV,HV_RN_AllVacc,t]*_PercMarried[HV_HRH_F3]/100*_NumPart[HV_HRH_F3,t]*_SexActs[HV_HRH_F3,t]
      +_Adults[HV_Female,HV_IDU,HV_AllHIV,HV_RN_AllVacc,t]*_PercMarried[HV_IDU_F3]/100*_NumPart[HV_IDU_F3,t]*_SexActs[HV_IDU_F3,t];
      (* Number of low risk SexActs men are having *)
      TotalSexActsM :=_Adults[HV_Male,HV_LRH,HV_AllHIV,HV_RN_AllVacc,t]*_NumPart[HV_LRH,t]*_SexActs[HV_LRH,t]+
      +_Adults[HV_Male,HV_MRH,HV_AllHIV,HV_RN_AllVacc,t]*_PercMarried[HV_MRH]/100*_NumPart[HV_MRH,t]*_SexActs[HV_MRH,t]
      +_Adults[HV_Male,HV_HRH,HV_AllHIV,HV_RN_AllVacc,t]*_PercMarried[HV_HRH]/100*_NumPart[HV_HRH,t]*_SexActs[HV_HRH,t]
      +_Adults[HV_Male,HV_IDU,HV_AllHIV,HV_RN_AllVacc,t]*_PercMarried[HV_IDU]/100*_NumPart[HV_IDU,t]*_SexActs[HV_IDU,t]
      +_Adults[HV_Male,HV_MSM,HV_AllHIV,HV_RN_AllVacc,t]*_PercMarried[HV_MSM]/100*_NumPart[HV_MSM,t]*_SexActs[HV_MSM,t]
      +_Adults[HV_Male,HV_MSMLR,HV_AllHIV,HV_RN_AllVacc,t]*_PercMarried[HV_MSMLR]/100*_NumPart[HV_MSMLR,t]*_SexActs[HV_MSMLR,t]
      +_Adults[HV_Male,HV_MSMMR,HV_AllHIV,HV_RN_AllVacc,t]*_PercMarried[HV_MSMMR]/100*_NumPart[HV_MSMMR,t]*_SexActs[HV_MSMMR,t]
      +_Adults[HV_Male,HV_MSMHR,HV_AllHIV,HV_RN_AllVacc,t]*_PercMarried[HV_MSMHR]/100*_NumPart[HV_MSMHR,t]*_SexActs[HV_MSMHR,t]
      +_Adults[HV_Male,HV_MSMIDU,HV_AllHIV,HV_RN_AllVacc,t]*_PercMarried[HV_MSMIDU]/100*_NumPart[HV_MSMHR,t]*_SexActs[HV_MSMHR,t];
    end;

    if r = HV_MRH then
      begin
      (* Number of medium risk SexActs women are having *)
      TotalSexActsF :=_Adults[HV_Female,HV_MRH,HV_AllHIV,HV_RN_AllVacc,t]*_NumPart[HV_MRH_F3,t]*_SexActs[HV_MRH_F3,t];
      (* Number of medium risk SexActs men are having *)
      TotalSexActsM :=_Adults[HV_Male,HV_MRH,HV_AllHIV,HV_RN_AllVacc,t]*_NumPart[HV_MRH,t]*_SexActs[HV_MRH,t];
    end;

    if r = HV_HRH then
    begin
      (* Number of high risk SexActs women are having *)
      TotalSexActsF :=_Adults[HV_Female,HV_HRH,HV_AllHIV,HV_RN_AllVacc,t]*_NumPart[HV_HRH_F3,t]*_SexActs[HV_HRH_F3,t];
      (* Number of high risk SexActs men are having *)
      TotalSexActsM :=_Adults[HV_Male,HV_HRH,HV_AllHIV,HV_RN_AllVacc,t]*_NumPart[HV_HRH,t]*_SexActs[HV_HRH,t];
    end;

    if r = HV_LRH then
    begin
      if (_Adults[HV_Male,HV_LRH,HV_AllHIV,HV_RN_AllVacc,t]+
                _Adults[HV_Male,HV_MRH,HV_AllHIV,HV_RN_AllVacc,t]*
                _PercMarried[HV_MRH]/100+
                _Adults[HV_Male,HV_HRH,HV_AllHIV,HV_RN_AllVacc,t]*
                _PercMarried[HV_HRH]/100+
                _Adults[HV_Male,HV_IDU,HV_AllHIV,HV_RN_AllVacc,t]*
                _PercMarried[HV_IDU]/100+
                _Adults[HV_Male,HV_MSM,HV_AllHIV,HV_RN_AllVacc,t]*
                _PercMarried[HV_MSM]/100+
                _Adults[HV_Male,HV_MSMLR,HV_AllHIV,HV_RN_AllVacc,t]*
                _PercMarried[HV_MSMLR]/100+
                _Adults[HV_Male,HV_MSMMR,HV_AllHIV,HV_RN_AllVacc,t]*
                _PercMarried[HV_MSMMR]/100+
                _Adults[HV_Male,HV_MSMHR,HV_AllHIV,HV_RN_AllVacc,t]*
                _PercMarried[HV_MSMHR]/100+
                _Adults[HV_Male,HV_MSMIDU,HV_AllHIV,HV_RN_AllVacc,t]*
                _PercMarried[HV_MSMIDU]/100) = 0 then
        PrevM := 0
      else
        PrevM := (_Adults[HV_Male,HV_LRH,HV_AllHIV,HV_RN_AllVacc,t]*
                _CalcPrevalence[HV_Male,HV_LRH,t]+
                _Adults[HV_Male,HV_MRH,HV_AllHIV,HV_RN_AllVacc,t]*
                _CalcPrevalence[HV_Male,HV_MRH,t]*
                _PercMarried[HV_MRH]/100+
                _Adults[HV_Male,HV_HRH,HV_AllHIV,HV_RN_AllVacc,t]*
                _CalcPrevalence[HV_Male,HV_HRH,t]*
                _PercMarried[HV_HRH]/100+
                _Adults[HV_Male,HV_IDU,HV_AllHIV,HV_RN_AllVacc,t]*
                _CalcPrevalence[HV_Male,HV_IDU,t]*
                _PercMarried[HV_IDU]/100+
                _Adults[HV_Male,HV_MSM,HV_AllHIV,HV_RN_AllVacc,t]*
                _CalcPrevalence[HV_Male,HV_MSM,t]*
                _PercMarried[HV_MSM]/100+
                _Adults[HV_Male,HV_MSMLR,HV_AllHIV,HV_RN_AllVacc,t]*
                _CalcPrevalence[HV_Male,HV_MSMLR,t]*
                _PercMarried[HV_MSMLR]/100+
                _Adults[HV_Male,HV_MSMMR,HV_AllHIV,HV_RN_AllVacc,t]*
                _CalcPrevalence[HV_Male,HV_MSMMR,t]*
                _PercMarried[HV_MSMMR]/100+
                _Adults[HV_Male,HV_MSMHR,HV_AllHIV,HV_RN_AllVacc,t]*
                _CalcPrevalence[HV_Male,HV_MSMHR,t]*
                _PercMarried[HV_MSMHR]/100+
                _Adults[HV_Male,HV_MSMIDU,HV_AllHIV,HV_RN_AllVacc,t]*
                _CalcPrevalence[HV_Male,HV_MSMIDU,t]*
                _PercMarried[HV_MSMIDU]/100
                )/
               (_Adults[HV_Male,HV_LRH,HV_AllHIV,HV_RN_AllVacc,t]+
                _Adults[HV_Male,HV_MRH,HV_AllHIV,HV_RN_AllVacc,t]*
                _PercMarried[HV_MRH]/100+
                _Adults[HV_Male,HV_HRH,HV_AllHIV,HV_RN_AllVacc,t]*
                _PercMarried[HV_HRH]/100+
                _Adults[HV_Male,HV_IDU,HV_AllHIV,HV_RN_AllVacc,t]*
                _PercMarried[HV_IDU]/100+
                _Adults[HV_Male,HV_MSM,HV_AllHIV,HV_RN_AllVacc,t]*
                _PercMarried[HV_MSM]/100+
                _Adults[HV_Male,HV_MSMLR,HV_AllHIV,HV_RN_AllVacc,t]*
                _PercMarried[HV_MSMLR]/100+
                _Adults[HV_Male,HV_MSMMR,HV_AllHIV,HV_RN_AllVacc,t]*
                _PercMarried[HV_MSMMR]/100+
                _Adults[HV_Male,HV_MSMHR,HV_AllHIV,HV_RN_AllVacc,t]*
                _PercMarried[HV_MSMHR]/100+
                _Adults[HV_Male,HV_MSMIDU,HV_AllHIV,HV_RN_AllVacc,t]*
                _PercMarried[HV_MSMIDU]/100);

      if (_Adults[HV_Female,HV_LRH,HV_AllHIV,HV_RN_AllVacc,t]+
                _Adults[HV_Female,HV_MRH,HV_AllHIV,HV_RN_AllVacc,t]*
                _PercMarried[HV_MRH_F3]/100+
                _Adults[HV_Female,HV_HRH,HV_AllHIV,HV_RN_AllVacc,t]*
                _PercMarried[HV_HRH_F3]/100+
                _Adults[HV_Female,HV_IDU,HV_AllHIV,HV_RN_AllVacc,t]*
                _PercMarried[HV_IDU_F3]/100) = 0 then
        PrevF := 0
      else
        PrevF := (_Adults[HV_Female,HV_LRH,HV_AllHIV,HV_RN_AllVacc,t]*
                _CalcPrevalence[HV_Female,HV_LRH,t]+
                _Adults[HV_Female,HV_MRH,HV_AllHIV,HV_RN_AllVacc,t]*
                _CalcPrevalence[HV_Female,HV_MRH,t]*
                _PercMarried[HV_MRH_F3]/100+
                _Adults[HV_Female,HV_HRH,HV_AllHIV,HV_RN_AllVacc,t]*
                _CalcPrevalence[HV_Female,HV_HRH,t]*
                _PercMarried[HV_HRH_F3]/100+
                _Adults[HV_Female,HV_IDU,HV_AllHIV,HV_RN_AllVacc,t]*
                _CalcPrevalence[HV_Female,HV_IDU,t]*
                _PercMarried[HV_IDU_F3]/100)/
               (_Adults[HV_Female,HV_LRH,HV_AllHIV,HV_RN_AllVacc,t]+
                _Adults[HV_Female,HV_MRH,HV_AllHIV,HV_RN_AllVacc,t]*
                _PercMarried[HV_MRH_F3]/100+
                _Adults[HV_Female,HV_HRH,HV_AllHIV,HV_RN_AllVacc,t]*
                _PercMarried[HV_HRH_F3]/100+
                _Adults[HV_Female,HV_IDU,HV_AllHIV,HV_RN_AllVacc,t]*
                _PercMarried[HV_IDU_F3]/100);

      PrevML := PrevM;
      PrevFL := PrevF;

      (* Infectiousness multipliers *)
      {rMultM := (_Adults[HV_Male,HV_LRH,HV_AllHIV,HV_RN_AllVacc,t]*
                _rMultAll[HV_Male,HV_LRH,t]+
                _Adults[HV_Male,HV_MRH,HV_AllHIV,HV_RN_AllVacc,t]*
                _rMultAll[HV_Male,HV_MRH,t]*
                _PercMarried[HV_MRH]/100+
                _Adults[HV_Male,HV_HRH,HV_AllHIV,HV_RN_AllVacc,t]*
                _rMultAll[HV_Male,HV_HRH,t]*
                _PercMarried[HV_HRH]/100+
                _Adults[HV_Male,HV_IDU,HV_AllHIV,HV_RN_AllVacc,t]*
                _rMultAll[HV_Male,HV_IDU,t]*
                _PercMarried[HV_IDU]/100+
                _Adults[HV_Male,HV_MSM,HV_AllHIV,HV_RN_AllVacc,t]*
                _rMultAll[HV_Male,HV_MSM,t]*
                _PercMarried[HV_MSM]/100+
                _Adults[HV_Male,HV_MSMLR,HV_AllHIV,HV_RN_AllVacc,t]*
                _rMultAll[HV_Male,HV_MSMLR,t]*
                _PercMarried[HV_MSMLR]/100+
                _Adults[HV_Male,HV_MSMMR,HV_AllHIV,HV_RN_AllVacc,t]*
                _rMultAll[HV_Male,HV_MSMMR,t]*
                _PercMarried[HV_MSMMR]/100+
                _Adults[HV_Male,HV_MSMHR,HV_AllHIV,HV_RN_AllVacc,t]*
                _rMultAll[HV_Male,HV_MSMHR,t]*
                _PercMarried[HV_MSMHR]/100+
                _Adults[HV_Male,HV_MSMIDU,HV_AllHIV,HV_RN_AllVacc,t]*
                _rMultAll[HV_Male,HV_MSMIDU,t]*
                _PercMarried[HV_MSMIDU]/100
                )/
               (_Adults[HV_Male,HV_LRH,HV_AllHIV,HV_RN_AllVacc,t]+
                _Adults[HV_Male,HV_MRH,HV_AllHIV,HV_RN_AllVacc,t]*
                _PercMarried[HV_MRH]/100+
                _Adults[HV_Male,HV_HRH,HV_AllHIV,HV_RN_AllVacc,t]*
                _PercMarried[HV_HRH]/100+
                _Adults[HV_Male,HV_IDU,HV_AllHIV,HV_RN_AllVacc,t]*
                _PercMarried[HV_IDU]/100+
                _Adults[HV_Male,HV_MSM,HV_AllHIV,HV_RN_AllVacc,t]*
                _PercMarried[HV_MSM]/100+
                _Adults[HV_Male,HV_MSMLR,HV_AllHIV,HV_RN_AllVacc,t]*
                _PercMarried[HV_MSMLR]/100+
                _Adults[HV_Male,HV_MSMMR,HV_AllHIV,HV_RN_AllVacc,t]*
                _PercMarried[HV_MSMMR]/100+
                _Adults[HV_Male,HV_MSMHR,HV_AllHIV,HV_RN_AllVacc,t]*
                _PercMarried[HV_MSMHR]/100+
                _Adults[HV_Male,HV_MSMIDU,HV_AllHIV,HV_RN_AllVacc,t]*
                _PercMarried[HV_MSMIDU]/100);}

      {rMultF := (_Adults[HV_Female,HV_LRH,HV_AllHIV,HV_RN_AllVacc,t]*
                _rMultAll[HV_Female,HV_LRH,t]+
                _Adults[HV_Female,HV_MRH,HV_AllHIV,HV_RN_AllVacc,t]*
                _rMultAll[HV_Female,HV_MRH,t]*
                _PercMarried[HV_MRH_F1]/100+
                _Adults[HV_Female,HV_HRH,HV_AllHIV,HV_RN_AllVacc,t]*
                _rMultAll[HV_Female,HV_HRH,t]*
                _PercMarried[HV_HRH_F1]/100+
                _Adults[HV_Female,HV_IDU,HV_AllHIV,HV_RN_AllVacc,t]*
                _rMultAll[HV_Female,HV_IDU,t]*
                _PercMarried[HV_IDU_F1]/100)/
               (_Adults[HV_Female,HV_LRH,HV_AllHIV,HV_RN_AllVacc,t]+
                _Adults[HV_Female,HV_MRH,HV_AllHIV,HV_RN_AllVacc,t]*
                _PercMarried[HV_MRH_F1]/100+
                _Adults[HV_Female,HV_HRH,HV_AllHIV,HV_RN_AllVacc,t]*
                _PercMarried[HV_HRH_F1]/100+
                _Adults[HV_Female,HV_IDU,HV_AllHIV,HV_RN_AllVacc,t]*
                _PercMarried[HV_IDU_F1]/100);}
    end
    else {medium risk, high risk}
    begin
      PrevM := _CalcPrevalence[HV_Male,r,t];
      PrevF := _CalcPrevalence[HV_Female,r,t];

      (* assume HIV prevalence among LR to be the lowest risk faced *)
      PrevM:=max(PrevML,PrevM);
      PrevF:=max(PrevFL,PrevF);

      {* Infectiousness multipliers *}
      //rMultM:=_rMultAll[HV_Male,r,t];
      //rMultF:=_rMultAll[HV_Female,r,t];
    end;

    (* This can be used to overide the risk-group multipliers with the global multiplier *)
    rMultM:=_rMultAll[HV_Bothsexes,HV_AllRisk,t];
    rMultF:=_rMultAll[HV_Bothsexes,HV_AllRisk,t];

    (* Sex acts acts ratio, toggled by _BalanceSexActs *)
    SexActsRatioM := 1;
    SexActsRatioF := 1;
    (* assume equal adjustment for men and women *)
    TotalSexActs:=(TotalSexActsM+TotalSexActsF)/2;
    if(TotalSexActsM>0) then
        SexActsRatioM := (TotalSexActs / TotalSexActsM);
    if(TotalSexActsF>0) then
        SexActsRatioF := (TotalSexActs / TotalSexActsF);

    (* Toggle balancing *)
    if bBalanceSexActs = FALSE then
    begin
      SexActsRatioM :=  1;{no balancing}
      SexActsRatioF :=  1;{no balancing}
    end;

    (* Assuming IDU have medium risk sexual behavior, and will be exposed to same risk as HV_MRH *)
    if(r=HV_MRH) then
    begin
      PrevM_IDU:=PrevM;
      PrevF_IDU:=PrevF;
      rMultM_IDU:=rMultM;
      rMultF_IDU:=rMultF;
      SexActsRatioM_IDU:=SexActsRatioM;
      SexActsRatioF_IDU:=SexActsRatioF;
    end;

    (* Male IDU Transmission is based on male IDU prevalence, but female IDU is
       also based on male IDU prevalence *)
    if not(r = HV_IDU) then
    begin
      (* men *)
      value := 1-(Power(PrevF*
               (Power((1-_TransHIVF*rMultF*
               ((1 - _MalesCircum[t]) + (1-_RedWhenCircum[HV_Susceptibility]) * _MalesCircum[t])*
               (1+(_TransMultSTI-1)*_STIPrev[r,t])*
               (1-_CondomPercent[r,t]*_CondomEff)*
               (1 - _PrEPPercent[HV_Male,r,t]*_PrEPEff[HV_Male,r,t])),
               _SexActs[r,t]*SexActsRatioM))+
               (1-PrevF),_NumPart[r,t]));


       value := value * _Adults[HV_Male,r,HV_Negative,HV_RN_UnV,t];
      _NewInfectionsT[HV_Male,r,HV_RN_UnV] := value;

      _NewInfectionsT[HV_Male,r,HV_RN_Take] := 0;

      value := 1-(Power(PrevF*
                 (Power((1-_TransHIVF*(1 -
                 (_VaccineEffectiveness[HV_RN_Efficacy]/100))*rMultF*
                 ((1- _MalesCircum[t]) + (1-_RedWhenCircum[HV_Susceptibility]) * _MalesCircum[t])*
                 (1+(_TransMultSTI-1)*
                 _STIPrev[r,t])*
                 (1-_CondomPercent[r,t]*_CondomEff) *
                 (1 - _PrEPPercent[HV_Male,r,t]*_PrEPEff[HV_Male,r,t])),
                 _SexActs[r,t]*SexActsRatioM))+
                 (1-PrevF),_NumPart[r,t]));



      value := value * _Adults[HV_Male,r,HV_Negative,HV_RN_Partial,t];
      _NewInfectionsT[HV_Male,r,HV_RN_Partial] := value;


      value := 1-(Power(PrevF*
                (Power((1-_TransHIVF*rMultF*
                ((1- _MalesCircum[t]) + (1-_RedWhenCircum[HV_Susceptibility]) * _MalesCircum[t])*
                (1+(_TransMultSTI-1)*_STIPrev[r,t])*
                (1-_CondomPercent[r,t]*_CondomEff)*
                (1 - _PrEPPercent[HV_Male,r,t]*_PrEPEff[HV_Male,r,t]) ),
                _SexActs[r,t]*SexActsRatioM))+
                (1-PrevF),_NumPart[r,t]));

      value := value * _Adults[HV_Male,r,HV_Negative,HV_RN_NoProt,t];
      _NewInfectionsT[HV_Male,r,HV_RN_NoProt] := value;

    end;
    (* women *)
    value := 1-(Power(PrevM*
              (Power((1-_TransHIVF * _TransMultM *rMultM*
              ((1- _MalesCircum[t]) + (1-_RedWhenCircum[HV_Infect]) * _MalesCircum[t])*
              (1+(_TransMultSTI-1)*_STIPrev[r + HV_Total,t])*
              (1-_CondomPercent[r,t]*_CondomEff)*
              (1 - _PrEPPercent[HV_Female,r,t]*_PrEPEff[HV_Female,r,t])),
              _SexActs[r + HV_Total,t]*SexActsRatioF))+
              (1-PrevM),_NumPart[r + HV_Total,t]));

    value := value * _Adults[HV_Female,r,HV_Negative,HV_RN_UnV,t];
    _NewInfectionsT[HV_Female,r,HV_RN_UnV] := value;

    _NewInfectionsT[HV_Female,r,HV_RN_Take] := 0;

    value := 1-(Power(PrevM*
              (Power((1-_TransHIVF * _TransMultM *(1 -
              (_VaccineEffectiveness[HV_RN_Efficacy]/100))*rMultM*
              ((1- _MalesCircum[t]) + (1-_RedWhenCircum[HV_Infect]) * _MalesCircum[t])*
              (1+(_TransMultSTI-1)*_STIPrev[r + HV_Total,t])*
              (1-_CondomPercent[r,t]*_CondomEff)*
              (1 - _PrEPPercent[HV_Female,r,t]*_PrEPEff[HV_Female,r,t])),
              _SexActs[r + HV_Total,t]*SexActsRatioF))+
              (1-PrevM),_NumPart[r + HV_Total,t]));

    value := value * _Adults[HV_Female,r,HV_Negative,HV_RN_Partial,t];
    _NewInfectionsT[HV_Female,r,HV_RN_Partial] := value;

    value := 1-(Power(PrevM*
             (Power((1-_TransHIVF * _TransMultM *rMultM*
             ((1- _MalesCircum[t]) + (1-_RedWhenCircum[HV_Infect]) * _MalesCircum[t])*
             (1+(_TransMultSTI-1)*_STIPrev[r + HV_Total,t])*
             (1-_CondomPercent[r,t]*_CondomEff)*
             (1 - _PrEPPercent[HV_Female,r,t]*_PrEPEff[HV_Female,r,t])),
             _SexActs[r + HV_Total,t]*SexActsRatioF))+
             (1-PrevM),_NumPart[r + HV_Total,t]));

    value := value * _Adults[HV_Female,r,HV_Negative,HV_RN_NoProt,t];
      _NewInfectionsT[HV_Female,r,HV_RN_NoProt] := value;
  end; { r loop }

   (* Add new infections due to contaminated blood *)
  if (t > 1) and (HV.GetBloodInfection(t) > 0) then
    for s := HV_Male to HV_Female do
      for v := HV_RN_UnV to HV_RN_NoProt do
        if (_Adults[HV_Male,HV_LRH,HV_Negative,HV_RN_AllVacc,t] + _Adults[HV_Female,HV_LRH,HV_Negative,HV_RN_AllVacc,t]) > 0 then
          _NewInfectionsT[s,HV_LRH,v] := _NewInfectionsT[s,HV_LRH,v]
                                        + HV.GetBloodInfection(t)
                                        * _Adults[s,HV_LRH,HV_Negative,v,t] /
                                        (_Adults[HV_Male,HV_LRH,HV_Negative,HV_RN_AllVacc,t]
                                        + _Adults[HV_Female,HV_LRH,HV_Negative,HV_RN_AllVacc,t]);

  for r := HV_MSM to HV_MSMIDU do
  begin
    (* MSM multiplier calculated among MSM groups only *)
    rMultM:=_rMultAll[HV_Male,r,t];
    //rMultM:=_rMultAll[HV_Bothsexes,HV_AllRisk,t];

    (* assume HIV prevalence among LR men to be the lowest risk faced by MSM *)
    PrevM := max(_CalcPrevalence[HV_Male,HV_LRH,t],_CalcPrevalence[HV_Male,r,t]);
    value := 1-(Power(PrevM*
              (Power((1-_TransHIVF*rMultM*
              ((1- _MalesCircum[t]) + (1-_RedWhenCircum[HV_Infect]) * _MalesCircum[t])*
              _TransMultMSM *
              (1+(_TransMultSTI-1)*_STIPrev[r,t])*
              (1-_CondomPercent[r,t]*_CondomEff)*
              (1 - _PrEPPercent[HV_Male,r,t]*_PrEPEff[HV_Male,r,t])),
              _SexActs[r,t]*(1)))+
              (1-PrevM),_NumPart[r,t]));

    value := value * _Adults[HV_Male,r,HV_Negative,HV_RN_UnV,t];
    _NewInfectionsT[HV_Male,r,HV_RN_UnV] := value;

    _NewInfectionsT[HV_Male,r,HV_RN_Take] := 0;

    value := 1-(Power(PrevM*
               (Power((1-_TransHIVF*(1 -
               (_VaccineEffectiveness[HV_RN_Efficacy]/100))*rMultM* _TransMultMSM *
               ((1- _MalesCircum[t]) + (1-_RedWhenCircum[HV_Infect]) * _MalesCircum[t])*
               (1+(_TransMultSTI-1)*_STIPrev[r,t])*
               (1-_CondomPercent[r,t]*_CondomEff)*
               (1 - _PrEPPercent[HV_Male,r,t]*_PrEPEff[HV_Male,r,t])),
               _SexActs[r,t]*(1)))+
               (1-PrevM),_NumPart[r,t]));

    value := value * _Adults[HV_Male,r,HV_Negative,HV_RN_Partial,t];
    _NewInfectionsT[HV_Male,r,HV_RN_Partial] := value;

    value := 1-(Power(PrevM*
              (Power((1-_TransHIVF*rMultM* _TransMultMSM *
              ((1- _MalesCircum[t]) + (1-_RedWhenCircum[HV_Infect]) * _MalesCircum[t])*
              (1+(_TransMultSTI-1)*_STIPrev[r,t])*
              (1-_CondomPercent[r,t]*_CondomEff)*
              (1 - _PrEPPercent[HV_Male,r,t]*_PrEPEff[HV_Male,r,t])),
              _SexActs[r,t]*(1)))+
              (1-PrevM),_NumPart[r,t]));

    value := value * _Adults[HV_Male,r,HV_Negative,HV_RN_NoProt,t];
    _NewInfectionsT[HV_Male,r,HV_RN_NoProt] := value;
  end; { r loop }

  r := HV_IDU;
  (* Assuming IDU have medium risk sexual behavior, see r=HV_MRH above *)
  PrevM:=PrevM_IDU;
  PrevF:=PrevF_IDU;
  SexActsRatioM:=SexActsRatioM_IDU;
  SexActsRatioF:=SexActsRatioF_IDU;
  (* This can be used to overide the risk-group multipliers with the global multiplier *)
  rMultM:=_rMultAll[HV_Bothsexes,HV_AllRisk,t];
  rMultF:=_rMultAll[HV_Bothsexes,HV_AllRisk,t];

  (* Prevalence of HIV among all IDU, not that MSM
     All IDU face same risk
     IDU contributes to HIV risk for IDU *)
  if _Adults[HV_Male,HV_IDU,HV_AllHIV,HV_RN_AllVacc,t] +
    _Adults[HV_Female,HV_IDU,HV_AllHIV,HV_RN_AllVacc,t] +
    _Adults[HV_Male,HV_MSMIDU,HV_AllHIV,HV_RN_AllVacc,t] = 0
  then
   PrevB := 0
  else
    PrevB := (_CalcPrevalence[HV_Male,HV_IDU,t] *
    _Adults[HV_Male,HV_IDU,HV_AllHIV,HV_RN_AllVacc,t] +
    _CalcPrevalence[HV_Female,HV_IDU,t] *
    _Adults[HV_Female,HV_IDU,HV_AllHIV,HV_RN_AllVacc,t]+
    _CalcPrevalence[HV_Male,HV_MSMIDU,t] *
    _Adults[HV_Male,HV_MSMIDU,HV_AllHIV,HV_RN_AllVacc,t]) /
    (_Adults[HV_Male,HV_IDU,HV_AllHIV,HV_RN_AllVacc,t] +
    _Adults[HV_Female,HV_IDU,HV_AllHIV,HV_RN_AllVacc,t]+
    _Adults[HV_Male,HV_MSMIDU,HV_AllHIV,HV_RN_AllVacc,t]);

  if _PerIDUSharing[t] > 0 then
    PrevB := PrevB / _PerIDUSharing[t];
    (* assume HIV prevalence among LR to be the lowest risk faced by IDU *)
  PrevB:=max(PrevB,min(_CalcPrevalence[HV_Male,HV_LRH,t],_CalcPrevalence[HV_Female,HV_LRH,t]));
    (* end of IDU risk calculation *)

  for s := HV_Male to HV_Female do
  begin
    SusceptibleIDU := _Adults[s,HV_IDU,HV_Negative,HV_RN_UnV,t];
    if SusceptibleIDU > (_Adults[s,HV_IDU,HV_AllHIV,HV_RN_UnV,t] * _PerIDUSharing[t]) -
                        (_Adults[s,HV_IDU,HV_AllHIV,HV_RN_UnV,t] -
                         _Adults[s,HV_IDU,HV_Negative,HV_RN_UnV,t]) then
      SusceptibleIDU := (_Adults[s,HV_IDU,HV_AllHIV,HV_RN_UnV,t] * _PerIDUSharing[t]) -
                           (_Adults[s,HV_IDU,HV_AllHIV,HV_RN_UnV,t] -
                            _Adults[s,HV_IDU,HV_Negative,HV_RN_UnV,t]);

    (* Sexual transmission for IDU HV_RN_UnV
       Assuming IDU have medium risk sexual behavior *)
//      _ForceSexTransmIDU := 0;
    if s = HV_Male then
      _ForceSexTransmIDU := 1-(Power(PrevF*
                 (Power((1-_TransHIVF*rMultF*
                 ((1 - _MalesCircum[t]) + (1-_RedWhenCircum[HV_Susceptibility]) * _MalesCircum[t])*
                 (1+(_TransMultSTI-1)*_STIPrev[HV_IDU,t])*
                 (1-_CondomPercent[HV_IDU,t]*_CondomEff)*
                 (1 - _PrEPPercent[s,HV_IDU,t]*_PrEPEff[s,HV_IDU,t])),
                 _SexActs[HV_IDU,t]*SexActsRatioM))+
                 (1-PrevF),_NumPart[HV_IDU,t]))
     else
     begin
      _ForceSexTransmIDU := 1-(Power(PrevM*
                (Power((1-_TransHIVF * _TransMultM *rMultM*
                ((1- _MalesCircum[t]) + (1-_RedWhenCircum[HV_Infect]) * _MalesCircum[t])*
                (1+(_TransMultSTI-1)*_STIPrev[HV_IDU + HV_Total,t])*
                (1-_CondomPercent[HV_IDU,t]*_CondomEff)*
                (1 - _PrEPPercent[s,HV_IDU,t]*_PrEPEff[s,HV_IDU,t])),
                _SexActs[HV_IDU + HV_Total,t]*SexActsRatioF))+
                (1-PrevM),_NumPart[HV_IDU + HV_Total,t]));
     end;

     _ForceNeedleTransmIDU:=0;
     (* can build in a needle sharing mechanism here *)
     //if(t>=_EpidemicStYrIDU) then
     //   _ForceNeedleTransmIDU:= _NeedleShareRate[t] * PrevB *_TransHIVNeedle
     //else
     //   _ForceNeedleTransmIDU:=0;
     //end;

     value := _ForceInf[s, t] * _rMultAll[HV_Bothsexes,HV_AllRisk,t] * PrevB * SusceptibleIDU *(1 - _PrEPPercent[s,HV_IDU,t]*_PrEPEff[s,HV_IDU,t]) +
             _ForceNeedleTransmIDU*SusceptibleIDU;
     (* all IDU susceptible to sexual transmission *)
     SusceptibleIDU := _Adults[s,HV_IDU,HV_Negative,HV_RN_UnV,t];
     value := value+_ForceSexTransmIDU*SusceptibleIDU;

    if value < 0 then
      value := 0;
    _NewInfectionsT[s,r,HV_RN_UnV] := value;

    _NewInfectionsT[s,r,HV_RN_Take] := 0;

    SusceptibleIDU := _Adults[s,HV_IDU,HV_Negative,HV_RN_Partial,t];
    if SusceptibleIDU > (_Adults[s,HV_IDU,HV_AllHIV,HV_RN_Partial,t] * _PerIDUSharing[t]) -
                        (_Adults[s,HV_IDU,HV_AllHIV,HV_RN_Partial,t] -
                         _Adults[s,HV_IDU,HV_Negative,HV_RN_Partial,t]) then
      SusceptibleIDU := (_Adults[s,HV_IDU,HV_AllHIV,HV_RN_Partial,t] * _PerIDUSharing[t]) -
                           (_Adults[s,HV_IDU,HV_AllHIV,HV_RN_Partial,t] -
                            _Adults[s,HV_IDU,HV_Negative,HV_RN_Partial,t]);

    (* Sexual transmission for IDU HV_RN_Partial *)
//      _ForceSexTransmIDU := 0;
    if s = HV_Male then
      _ForceSexTransmIDU := 1-(Power(PrevF*
                   (Power((1-_TransHIVF*(1 -
                   (_VaccineEffectiveness[HV_RN_Efficacy]/100))*rMultF*
                   ((1- _MalesCircum[t]) + (1-_RedWhenCircum[HV_Susceptibility]) * _MalesCircum[t])*
                   (1+(_TransMultSTI-1)*_STIPrev[HV_IDU,t])*
                   (1-_CondomPercent[HV_IDU,t]*_CondomEff) *
                   (1 - _PrEPPercent[s,HV_IDU,t]*_PrEPEff[s,HV_IDU,t])),
                   _SexActs[HV_IDU,t]*SexActsRatioM))+
                   (1-PrevF),_NumPart[HV_IDU,t]))
     else
     begin
      _ForceSexTransmIDU := 1-(Power(PrevM*
                (Power((1-_TransHIVF * _TransMultM *(1 -
                (_VaccineEffectiveness[HV_RN_Efficacy]/100))*rMultM*
                ((1- _MalesCircum[t]) + (1-_RedWhenCircum[HV_Infect]) * _MalesCircum[t])*
                (1+(_TransMultSTI-1)*_STIPrev[HV_IDU + HV_Total,t])*
                (1-_CondomPercent[HV_IDU,t]*_CondomEff)*
                (1 - _PrEPPercent[s,HV_IDU,t]*_PrEPEff[s,HV_IDU,t])),
                _SexActs[HV_IDU + HV_Total,t]*SexActsRatioF))+
                (1-PrevM),_NumPart[HV_IDU + HV_Total,t]));
     end;

    value := _ForceInf[s, t] * _rMultAll[HV_Bothsexes,HV_AllRisk,t] * PrevB * SusceptibleIDU *(1 - _PrEPPercent[s,HV_IDU,t]*_PrEPEff[s,HV_IDU,t]) * (1 - _VaccineEffectiveness[HV_RN_Efficacy]/100) +
             _ForceNeedleTransmIDU*SusceptibleIDU;
     (* all IDU susceptible to sexual transmission *)
    SusceptibleIDU := _Adults[s,HV_IDU,HV_Negative,HV_RN_Partial,t];
    value:=value+_ForceSexTransmIDU*SusceptibleIDU;

    if value < 0 then
      value := 0;
    _NewInfectionsT[s,r,HV_RN_Partial] := value;

    SusceptibleIDU := _Adults[s,HV_IDU,HV_Negative,HV_RN_NoProt,t];
    if SusceptibleIDU > (_Adults[s,HV_IDU,HV_AllHIV,HV_RN_NoProt,t] * _PerIDUSharing[t]) -
                        (_Adults[s,HV_IDU,HV_AllHIV,HV_RN_NoProt,t] -
                         _Adults[s,HV_IDU,HV_Negative,HV_RN_NoProt,t]) then
      SusceptibleIDU := (_Adults[s,HV_IDU,HV_AllHIV,HV_RN_NoProt,t] * _PerIDUSharing[t]) -
                           (_Adults[s,HV_IDU,HV_AllHIV,HV_RN_NoProt,t] -
                            _Adults[s,HV_IDU,HV_Negative,HV_RN_NoProt,t]);

    (* Sexual transmission for IDU HV_RN_NoProt *)
//      _ForceSexTransmIDU := 0;
    if s = HV_Male then
      _ForceSexTransmIDU := 1-(Power(PrevF*
                  (Power((1-_TransHIVF*rMultF*
                  ((1- _MalesCircum[t]) + (1-_RedWhenCircum[HV_Susceptibility]) * _MalesCircum[t])*
                  (1+(_TransMultSTI-1)*_STIPrev[HV_IDU,t])*
                  (1-_CondomPercent[HV_IDU,t]*_CondomEff)*
                  (1 - _PrEPPercent[s,HV_IDU,t]*_PrEPEff[s,HV_IDU,t])),
                  _SexActs[HV_IDU,t]*SexActsRatioM))+
                  (1-PrevF),_NumPart[HV_IDU,t]))
     else
     begin
       _ForceSexTransmIDU := 1-(Power(PrevM*
               (Power((1-_TransHIVF * _TransMultM *rMultM*
               ((1- _MalesCircum[t]) + (1-_RedWhenCircum[HV_Infect]) * _MalesCircum[t])*
               (1+(_TransMultSTI-1)*_STIPrev[HV_IDU + HV_Total,t])*
               (1-_CondomPercent[HV_IDU,t]*_CondomEff)*
               (1 - _PrEPPercent[s,HV_IDU,t]*_PrEPEff[s,HV_IDU,t])),
               _SexActs[HV_IDU + HV_Total,t]*SexActsRatioF))+
               (1-PrevM),_NumPart[HV_IDU + HV_Total,t]));
     end;

    value := _ForceInf[s, t] * _rMultAll[HV_Bothsexes,HV_AllRisk,t] * PrevB * (1 - _PrEPPercent[s,HV_IDU,t]*_PrEPEff[s,HV_IDU,t]) * SusceptibleIDU +
             _ForceNeedleTransmIDU*SusceptibleIDU;
    (* all IDU susceptible to sexual transmission *)
    SusceptibleIDU := _Adults[s,HV_IDU,HV_Negative,HV_RN_NoProt,t];
    value := value+_ForceSexTransmIDU*SusceptibleIDU;

    if value < 0 then
      value := 0;
    _NewInfectionsT[s,r,HV_RN_NoProt] := value;
  end;{s loop}

  (* add IDU infections for MSM:IDU, assuming same needle sharing risk as other IDU
     note that infection from sexual transmission already added
     same PrevB calculated for all IDU, applies to MSM:IDU *)
  _NewInfectionsT[HV_Male,HV_MSMIDU,HV_RN_Take] := 0;

  SusceptibleIDU := _Adults[HV_Male,HV_MSMIDU,HV_Negative,HV_RN_UnV,t];
  if SusceptibleIDU > (_Adults[HV_Male,HV_MSMIDU,HV_AllHIV,HV_RN_UnV,t] * _PerIDUSharing[t]) -
                      (_Adults[HV_Male,HV_MSMIDU,HV_AllHIV,HV_RN_UnV,t] -
                       _Adults[HV_Male,HV_MSMIDU,HV_Negative,HV_RN_UnV,t])
  then SusceptibleIDU := (_Adults[HV_Male,HV_MSMIDU,HV_AllHIV,HV_RN_UnV,t] * _PerIDUSharing[t]) -
                         (_Adults[HV_Male,HV_MSMIDU,HV_AllHIV,HV_RN_UnV,t] -
                          _Adults[HV_Male,HV_MSMIDU,HV_Negative,HV_RN_UnV,t]);

   _ForceNeedleTransmIDU:=0;
   {can build in a needle sharing mechanism here}

   value := _NewInfectionsT[HV_Male,HV_MSMIDU,HV_RN_UnV] +
            _ForceInf[HV_Male, t] * _rMultAll[HV_Bothsexes,HV_AllRisk,t] * PrevB * SusceptibleIDU * (1 - _PrEPPercent[HV_Male,HV_MSMIDU,t]*_PrEPEff[HV_Male,HV_MSMIDU,t]) +
            _ForceNeedleTransmIDU*SusceptibleIDU;

  if value < 0 then
    value := 0;
  _NewInfectionsT[HV_Male,HV_MSMIDU,HV_RN_UnV] := value;

  SusceptibleIDU := _Adults[HV_Male,HV_MSMIDU,HV_Negative,HV_RN_Partial,t];
  if SusceptibleIDU > (_Adults[HV_Male,HV_MSMIDU,HV_AllHIV,HV_RN_Partial,t] * _PerIDUSharing[t]) -
                      (_Adults[HV_Male,HV_MSMIDU,HV_AllHIV,HV_RN_Partial,t] -
                       _Adults[HV_Male,HV_MSMIDU,HV_Negative,HV_RN_Partial,t]) then
    SusceptibleIDU := (_Adults[HV_Male,HV_MSMIDU,HV_AllHIV,HV_RN_Partial,t] * _PerIDUSharing[t]) -
                         (_Adults[HV_Male,HV_MSMIDU,HV_AllHIV,HV_RN_Partial,t] -
                          _Adults[HV_Male,HV_MSMIDU,HV_Negative,HV_RN_Partial,t]);

  value := _NewInfectionsT[HV_Male,HV_MSMIDU,HV_RN_Partial] +
           _ForceInf[HV_Male, t] * _rMultAll[HV_Bothsexes,HV_AllRisk,t] * PrevB * SusceptibleIDU * (1 - _PrEPPercent[HV_Male,HV_MSMIDU,t]*_PrEPEff[HV_Male,HV_MSMIDU,t]) * (1 - _VaccineEffectiveness[HV_RN_Efficacy]/100) +
           _ForceNeedleTransmIDU*SusceptibleIDU;

  if value < 0 then
    value := 0;
  _NewInfectionsT[HV_Male,HV_MSMIDU,HV_RN_Partial] := value;

  SusceptibleIDU := _Adults[HV_Male,HV_MSMIDU,HV_Negative,HV_RN_NoProt,t];
  if SusceptibleIDU > (_Adults[HV_Male,HV_MSMIDU,HV_AllHIV,HV_RN_NoProt,t] * _PerIDUSharing[t]) -
                      (_Adults[HV_Male,HV_MSMIDU,HV_AllHIV,HV_RN_NoProt,t] -
                       _Adults[HV_Male,HV_MSMIDU,HV_Negative,HV_RN_NoProt,t]) then
    SusceptibleIDU := (_Adults[HV_Male,HV_MSMIDU,HV_AllHIV,HV_RN_NoProt,t] * _PerIDUSharing[t]) -
                         (_Adults[HV_Male,HV_MSMIDU,HV_AllHIV,HV_RN_NoProt,t] -
                          _Adults[HV_Male,HV_MSMIDU,HV_Negative,HV_RN_NoProt,t]);

  value := _NewInfectionsT[HV_Male,HV_MSMIDU,HV_RN_NoProt]+_ForceInf[HV_Male, t] * _rMultAll[HV_Bothsexes,HV_AllRisk,t] * PrevB * SusceptibleIDU * (1 - _PrEPPercent[HV_Male,HV_MSMIDU,t]*_PrEPEff[HV_Male,HV_MSMIDU,t]) +
           _ForceNeedleTransmIDU*SusceptibleIDU;

  if value < 0 then
    value := 0;
  _NewInfectionsT[HV_Male,HV_MSMIDU,HV_RN_NoProt] := value;
  (* end of IDU infection MSM:IDU *)

  value:=0;
  for s := HV_Male to HV_Female do
    for r := HV_LRH to HV_MSMIDU do
    for v := HV_RN_UnV to HV_RN_NoProt do
      value:=value+_NewInfectionsT[s,r,v];

  if(value>0) then
    Calc_NewInfections:=True;
end;

procedure THVProj.Calc_AddNewInfections(t,TS:integer);
var
  s,r,h,v : byte;
begin
  (*add new infections to Goals*)
  for s := HV_Male to HV_Female do
    for r := HV_LRH to HV_MSMIDU do
      for v := HV_RN_UnV to HV_RN_NoProt do
      begin
        if not((s = HV_Female) and (r >= HV_MSM)) then
        begin
            (* remove new infections from HIV- *)
            _Adults[s,r,HV_Negative,v,t] := _Adults[s,r,HV_Negative,v,t]-
                                          (1/TS)*_NewInfectionsT[s,r,v];
            (* add new infection to HIV+ categories *)
            for h := HV_Primary to HV_CD4_LT50 do
            _Adults[s,r,h,v,t] := _Adults[s,r,h,v,t]+
                                          (1/TS)*_NewInfByStage[s,h]*_NewInfectionsT[s,r,v];
         end;
      end;{s,r,v loop}

   (*add newinfections to AIM*)
    for s := HV_Male to HV_Female do
    begin
      //_NewHIV[s]:=0; set to zero in InnerLoopCalcs
      for r := HV_LRH to HV_MSMIDU do
        for v := HV_RN_UnV to HV_RN_NoProt do
          _NewHIV[s]:=_NewHIV[s]+(1/HV_TimeStep)*_NewInfectionsT[s,r,v];
    end;
end;

procedure THVProj.Calc_AddNewInfectionsAIM(t:integer);
var
s:byte;
begin
     for s := DP_Male to DP_Female do
         GetDPData(p).setGoalsNewHIVbySex(s,_NewHIV[s]);
end;

function THVProj.Calc_HasHIV(t:integer):boolean;
var
  s,r,h,v : byte;
  value   : double;
begin
  Calc_HasHIV := false;

  value:=0;
  for s := HV_Male to HV_Female do
    for r := HV_None to HV_MSMIDU do
      for h := HV_Primary to HV_CD4_LT50_ART do
        for v := HV_RN_UnV to HV_RN_NoProt do
          value:=value+_Adults[s,r,h,v,t];

  if(value>0) then
    Calc_HasHIV:=true;
end;

procedure THVProj.Calc_Prevalence(t,_ts:integer);
var
  s,r,h,v       : byte;
  value,value2  : double;
begin
  (* Calc prevalence *)
  for s := HV_Male to HV_Female do
  begin
    for r := HV_None to HV_MSMIDU do
    begin
      if not((s = HV_Female) and (r >= HV_MSM)) then
      begin
        value2 := 0;
        for v := HV_RN_UnV to HV_RN_NoProt do
          for h := HV_Negative to HV_CD4_LT50_ART do
            value2 := value2 + _AdultsTS[s,r,h,v,_ts];{use _AdultsTS}

        if value2 = 0 then
          value := 0
        else
        begin
          value := 0;
          for v := HV_RN_UnV to HV_RN_NoProt do
            for h := HV_Primary to HV_CD4_LT50_ART do
              value := value + _AdultsTS[s,r,h,v,_ts];{use _AdultsTS}
          value := value / value2;
        end;{if value2 = 0}

        if value < 0 then
          value := 0;
        if value > 1 then
          value := 1;
        _CalcPrevalence[s,r,t] := value;
      end; {if s and r statement}
    end;{r loop}
  end;{s loop}
end;

procedure THVProj.Calc_PopulationTotals(t:integer);
var
  s,r,h,v : byte;
  value   : double;
begin
  (*Sum to get total pops*)
  for s := HV_Male to HV_Female do
  begin
    value := 0;
    for r := HV_None to HV_MSMIDU do
    begin
      for h := HV_Negative to HV_CD4_LT50_ART do
      begin
        for v := HV_RN_UnV to HV_RN_NoProt do
        begin
          value := value + _Adults[s,r,h,v,t];
        end;
      end;
    end;
    _Populations[s,t] := value;
  end;
  value := _Populations[HV_Male,t] + _Populations[HV_Female,t];
  _Populations[HV_Bothsexes,t] := value;
end;

procedure THVProj.Calc_VaccinationTotals(t:integer);
var
  s,r,h : byte;
  value : double;
begin
  for s := HV_Male to HV_Female do
  begin
    for r := HV_None to HV_MSMIDU do
    begin
      for h := HV_Negative to HV_CD4_LT50 do
      begin
        value := _Adults[s,r,h,HV_RN_Take,t] +
                 _Adults[s,r,h,HV_RN_Partial,t] +
                 _Adults[s,r,h,HV_RN_NoProt,t];
        _Vaccinated[s,r,h,t] := value;
      end;
    end;
  end;

  for s := HV_Male to HV_Female do
  begin
    value := 0;
    for r := HV_None to HV_MSMIDU do
    begin
      for h := HV_Negative to HV_CD4_LT50 do
      begin
        value := value + _Vaccinated[s,r,h,t];
      end;
    end;
    _TotalVaccinated[s,t] := value;
  end;

  for s := HV_Male to HV_Female do
  begin
    for r := HV_None to HV_MSMIDU do
    begin
      value := 0;
      for h := HV_Negative to HV_CD4_LT50_ART do
          value := value + _Adults[s,r,h,HV_RN_UnV,t];
      _Adults[s,r,HV_AllHIV,HV_RN_UnV,t] := value;
    end;
  end;

  value := _TotalVaccinated[HV_Male,t] +
           _TotalVaccinated[HV_Female,t];
  _TotalVaccinated[HV_Bothsexes,t] := value;

  value := _Populations[HV_Bothsexes,t] - value;
  _Unvaccinated[t] := value;
end;

procedure THVProj.Calc_NewVaccinations(t:integer);
var
  s,r,h,h1,nr           : byte;
  TotalNewVaccinations,
  value,value2          : double;
begin
  TotalNewVaccinations := 0;

  for s := HV_Male to HV_Female do
  begin
    for r := HV_None to HV_MSMIDU do
    begin
      value := 0;

      (* [CDP] h1 := TN_SymptomaticHIV; *)
      h1 := HV_CD4_LT50;
      if HV_GetRNVaccineTargeting(p) = 1 then h1 := HV_Negative;
      value2 := 0;
      for h := HV_Negative to h1 do
      begin
        if not((s = HV_Female) and (r >= HV_MSM)) then
        begin
          if ( _Adults[s,HV_AllRisk,HV_AllHIV,HV_RN_AllVacc,t] = 0) or ( _VaccineEffectiveness[HV_RN_Duration] = 0) then
            value2 := 0
          else
          begin
            if HV_GetRNVaccineCovType(p) = HV_RN_Single then
              nr := HV_AllRisk
            else
            begin
              if s = HV_Female then
                nr := r + HV_Total
              else
                nr := r;
            end;
            if nr <= HV_IDU_F3 then
            begin
              value2 := _VacCoverage[nr,t]/100 *
                         {(1/HV_TimeStep)*} (_Adults[s,r,h,HV_RN_AllVacc,t] +
                            _Pop15rate[s,HV_DP_H_No,t] *
                           _Adults[s,r,h,HV_RN_AllVacc,t] / _Adults[s,HV_AllRisk,HV_AllHIV,HV_RN_AllVacc,t]) -
                           (_Adults[s,r,h,HV_RN_Take,t]+
                            _Adults[s,r,h,HV_RN_Partial,t]+
                            _Adults[s,r,h,HV_RN_NoProt,t]) * ( 1 - (1 / _VaccineEffectiveness[HV_RN_Duration])) ;

             if value2 < 0 then
              value2 := 0;
            end;
          end;
          value := value + value2;
        end;
      end; {h loop}
      _NewVaccinations[s,r,t] := value;
      TotalNewVaccinations := TotalNewVaccinations + value;
    end; {r loop}
  end; { s loop}
  _TotalNewVaccinations[t] := TotalNewVaccinations;
end;

//procedure THVProj.Calc_NewEntrants(t:integer);
//begin
//  (* Calculate population entering adulthood at age 15 *)
//  _Pop15rate[HV_Male,HV_DP_H_No,t] := (HV_GetDPPop(p,t,HV_Male,HV_DP_A10_14)+
//                         HV_GetDPPop(p,t,HV_Male,HV_DP_A15_19))/10;
//
//  _Pop15rate[HV_Female,HV_DP_H_No,t] := (HV_GetDPPop(p,t,HV_Female,HV_DP_A10_14)+
//                           HV_GetDPPop(p,t,HV_Female,HV_DP_A15_19))/10;
//end;

procedure THVProj.Calc_NewEntrantsAIM(t:integer);
var
  s,h,d : byte;
begin
  for s := HV_Male to HV_Female do
  begin
    (* The survival probabilities correct for the transition 14 to 15 year old HIV- *)
    _Pop15rate[s, HV_DP_H_No, t] := HV_GetDPAIDSAging15(p, t, s, HV_DP_H_No, HV_DP_D_HIVNeg);//*(1-_NonAIDSDeathRate[s,t-1]-0*_RateofAging[s,t]);

    for h := HV_CD4_GT500 to HV_CD4_LT50 do
    begin
      //HIV+, No ART
        _Pop15rate[s, h, t] := HV_GetDPAIDSAging15(p, t, s, h, HV_DP_D_HIVPos);//*(1-_NonAIDSDeathRate[s,t-1]-_mu[h+1,s]-0*_RateofHIVAging[s,h+1,t]);

       //HIV+, On ART
       _Pop15rate[s, h+10, t] := 0;
       for d := HV_DP_D_ARTlt6m to HV_DP_D_ARTgt12m do
        _Pop15rate[s, h+10, t] := _Pop15rate[s, h+10, t] + HV_GetDPAIDSAging15(p, t, s, h, d);//*(1-_NonAIDSDeathRate[s,t-1]-_alpha[h+8,s]-0*_RateofARTAging[s,h+8,t]);
    end;{h loop}
  end;
end;



procedure THVProj.Calc_MigrRate(t:integer);
var
  s,r,h,v : byte;
  PopSize   : array[HV_Male..HV_Female, HV_Neg..HV_Pos] of double;
begin

  (* Calculate migration rate denominators *)
  for s := HV_Male to HV_Female do
  begin
    PopSize[s,HV_Neg] := 0.0;
    PopSize[s,HV_Pos] := 0.0;
    for r := HV_None to HV_MSMIDU do
      for v := HV_RN_UnV to HV_RN_NoProt do
        PopSize[s,HV_Neg] := PopSize[s,HV_Neg] + _Adults[s,r,HV_Negative,v,t];
    PopSize[s,HV_Pos] := _Populations[s,t] - PopSize[s,HV_Neg];
  end;

  for s := HV_Male to HV_Female do
    for h := HV_Neg to HV_Pos do
    begin
      _MigrRate[t,s,h] := 0.0;
      if (PopSize[s,h] > 0) then
        _MigrRate[t,s,h] := _MigrNum[t,s,h] / PopSize[s,h];
    end;

  for s := HV_Male to HV_Female do
    for r := HV_None to HV_MSMIDU do
      for v := HV_RN_UnV to HV_RN_NoProt do
      begin
        h := HV_Negative;
        _Migrant[s,r,h,v] := _Adults[s,r,h,v,t] * _MigrRate[t,s,HV_Neg];
        if _Migrant[s,r,h,v] + _Adults[s,r,h,v,t] < 0 then
          _Migrant[s,r,h,v] := -_Adults[s,r,h,v,t];

        for h := HV_Primary to HV_CD4_LT50_ART do
        begin
          _Migrant[s,r,h,v] := _Adults[s,r,h,v,t] * _MigrRate[t,s,HV_Pos];
          if _Migrant[s,r,h,v] + _Adults[s,r,h,v,t] < 0 then
            _Migrant[s,r,h,v] := -_Adults[s,r,h,v,t];
        end;
      end;
end;

procedure THVProj.Calc_MigrNum(t:integer);
var
  s, a       : byte;
  migr, prev : double;
begin
  (* Calculate numbers of net migrants by HIV status from DemProj *)
  for s := HV_Male to HV_Female do
  begin
    _MigrNum[t,s,HV_Neg] := 0.0;
    _MigrNum[t,s,HV_Pos] := 0.0;
    for a := HV_DP_A15_19 to HV_DP_A45_49 do
    begin
      (* We use last year's HIV prevalence to approximate net migration of
         PLHIV because this time t prevalence is not known yet *)
      migr := (HV_GetDPMigrRate(p, t-1, s, a) + HV_GetDPMigrRate(p, t, s, a)) / 2;
      if HV_GetDPPop(p, t-1, s, a) > 0.0 then
        prev := HV_GetDPHIV(p, t-1, s, a) / HV_GetDPPop(p, t-1, s, a)
      else
        prev := 0.0;
      _MigrNum[t,s,HV_Neg] := _MigrNum[t,s,HV_Neg] + migr * (1.0 - prev);
      _MigrNum[t,s,HV_Pos] := _MigrNum[t,s,HV_Pos] + migr * prev;
    end;
  end;
end;

procedure THVProj.Calc_ProgressNoRiskHIVNeg(t,_ts:integer);
var
  s,v     : byte;
  value   : double;
begin
  for s := HV_Male to HV_Female do
  begin
    if ((max(_AgeFirstSex[s,t]-15,1)) -
            (_NonAIDSDeathRate[s,t]+
             _RateofAging[s,t]) *
             _Adults[s,HV_None,HV_Negative,HV_RN_UnV,t] +
             _Migrant[s,HV_None,HV_Negative,HV_RN_UnV] -
             _NewVaccinations[s,HV_None,t] *
             _Adults[s,HV_None,HV_Negative,HV_RN_UnV,t] = 0) or
            (_Adults[s,HV_None,HV_AllHIV,HV_RN_UnV,t] +
            (_Adults[s,HV_None,HV_Negative,HV_RN_Take,t] +
             _Adults[s,HV_None,HV_Negative,HV_RN_Partial,t] +
             _Adults[s,HV_None,HV_Negative,HV_RN_NoProt,t]) = 0) then
      value := 0
    else
      value := _Adults[s,HV_None,HV_Negative,HV_RN_UnV,t] +
              (1/HV_TimeStep)*(_Pop15rate[s,HV_DP_H_No,t] - _Adults[s,HV_None,HV_Negative,HV_RN_UnV,t]/
              (max(_AgeFirstSex[s,t]-15,1)) -
              (_NonAIDSDeathRate[s,t]+
               _RateofAging[s,t]) *
               _Adults[s,HV_None,HV_Negative,HV_RN_UnV,t]+
               _Migrant[s,HV_None,HV_Negative,HV_RN_UnV] -
               _NewVaccinations[s,HV_None,t] *
               _Adults[s,HV_None,HV_Negative,HV_RN_UnV,t] /
               _Adults[s,HV_None,HV_AllHIV,HV_RN_UnV,t] +
              (_Adults[s,HV_None,HV_Negative,HV_RN_Take,t] +
               _Adults[s,HV_None,HV_Negative,HV_RN_Partial,t] +
               _Adults[s,HV_None,HV_Negative,HV_RN_NoProt,t]) /
               _VaccineEffectiveness[HV_RN_Duration]);
    if value < 0 then value := 0;
    _Adults[s,HV_None,HV_Negative,HV_RN_UnV,t] := value;
  end; {s loop}

  (* At this point _Adults[s,HV_None,HV_Negative,HV_RN_UnV,t] is updated, so use _AdultsTS, which represents
     _Adults before the update, in certain references to _Adults below *)
  for s := HV_Male to HV_Female do
  begin
    for v := HV_RN_Take to HV_RN_NoProt do
    begin
      if (_AdultsTS[s,HV_None,HV_AllHIV,HV_RN_UnV,_ts] *
               _Veffect[_TNType,v] -
               _Adults[s,HV_None,HV_Negative,v,t]) = 0 then
        value := 0
      else if ((max(_AgeFirstSex[s,t] - 15,1)) -
              (_NonAIDSDeathRate[s,t]+
               _RateofAging[s,t]) *
               _Adults[s,HV_None,HV_Negative,v,t] -
               _Adults[s,HV_None,HV_Negative,v,t]) = 0 then
        value := 0
      else if (_VaccineEffectiveness[HV_RN_Duration] ) = 0 then
        value := 0
      else
      value := _Adults[s,HV_None,HV_Negative,v,t] +
                (1/HV_TimeStep)*(_NewVaccinations[s,HV_None,t] *
               _AdultsTS[s,HV_None,HV_Negative,HV_RN_UnV,_ts] / {use _AdultsTS}
               _AdultsTS[s,HV_None,HV_AllHIV,HV_RN_UnV,_ts] *   {use _AdultsTS}
               _Veffect[_TNType,v] -
               _Adults[s,HV_None,HV_Negative,v,t] /
              (max(_AgeFirstSex[s,t] - 15,1)) -
              (_NonAIDSDeathRate[s,t]+
               _RateofAging[s,t]) *
               _Adults[s,HV_None,HV_Negative,v,t] +
               _Migrant[s,HV_None,HV_Negative,v] -
               _Adults[s,HV_None,HV_Negative,v,t] /
               _VaccineEffectiveness[HV_RN_Duration]);

      if value < 0 then value := 0;
      _Adults[s,HV_None,HV_Negative,v,t] := value;
    end; {v loop}
  end; {s loop}
end;

procedure THVProj.Calc_ProgressNoRiskHIVPos(t:integer);
var
  s,h       : byte;
  value     : double;
begin
  for s := HV_Male to HV_Female do
  begin
    //not on ART
    for h := HV_CD4_GT500 to HV_CD4_LT50 do
    begin
      value := _Adults[s,HV_None,h,HV_RN_UnV,t] +
                (1/HV_TimeStep)*(_Pop15rate[s,h,t] -
                _Adults[s,HV_None,h,HV_RN_UnV,t]/max(_AgeFirstSex[s,t]-15,1) -
                _StageExits[s,HV_None,h,HV_RN_UnV]+
                _StageProgressors[s,HV_None,h-1,HV_RN_UnV]+
                _Migrant[s,HV_None,h,HV_RN_UnV]);

      if value < 0 then
        value := 0;
      _Adults[s,HV_None,h,HV_RN_UnV,t] := value;
    end;{h loop}

    //on ART
    for h := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
    begin

      value := _Adults[s,HV_None,h,HV_RN_UnV,t] +
                (1/HV_TimeStep)*(_Pop15rate[s,h,t] -
                _Adults[s,HV_None,h,HV_RN_UnV,t]/max(_AgeFirstSex[s,t]-15,1) -
                _StageExits[s,HV_None,h,HV_RN_UnV]+
                _StageProgressors[s,HV_None,h-1,HV_RN_UnV]+
                _Migrant[s,HV_None,h,HV_RN_UnV]);

      if value < 0 then
        value := 0;
      _Adults[s,HV_None,h,HV_RN_UnV,t] := value;
    end;{h loop}
  end; {s loop}
end;

procedure THVProj.Calc_ProgressHIVneg(t,_ts:integer);
var
  s,r,v,nr,r1,
  r2,r3,r4        : byte;
  value,value2    : double;
begin
  (*PopActiveHIV-*)
  for s := HV_Male to HV_Female do
  begin
    value2 := _Adults[s,HV_LRH,HV_AllHIV,HV_RN_UnV,t];
    if value2 <> 0 then
       value2 := _Adults[s,HV_LRH,HV_Negative,HV_RN_UnV,t] /
                 _Adults[s,HV_LRH,HV_AllHIV,HV_RN_UnV,t];

    if s = HV_Male then
    begin
      r1 := HV_MRH;
    //        r2 := HV_HRH_M;
    //        r3 := HV_IDU_M;
      r4 := HV_LRH;
    end
    else
    begin
      r1 := HV_MRH_F3;
    //        r2 := HV_HRH_F3;
    //        r3 := HV_IDU_F3;
      r4 := HV_LRH_F3;
    end;

    value := _Adults[s,HV_LRH,HV_Negative,HV_RN_UnV,t]+
             (1/HV_TimeStep)*(_AdultsTS[s,HV_None,HV_Negative,HV_RN_UnV,_ts] / //use _AdultsTS
             (max(_AgeFirstSex[s,t] - 15,1)) *
            _RiskGroupPercent[r4] +
             _Adults[s,HV_MRH,HV_Negative,HV_RN_UnV,t] *
             _BehaviorChangeRate[r1] -
             //
             {HV_HRH->HV_MRH, HV_IDU_HV_HRH, as is done for behavior transitions in HIV+ categories}
             //
             {_Adults[s,HV_HRH,HV_Negative,HV_RN_UnV,t] /
             _Behavior[r2,HV_AvgDur] +
             _Adults[s,HV_IDU,HV_Negative,HV_RN_UnV,t] /
             _Behavior[r3,HV_AvgDur]}
            (_NonAIDSDeathRate[s,t] +
             _RateofAging[s,t]) *
             _Adults[s,HV_LRH,HV_Negative,HV_RN_UnV,t] +
             _Migrant[s,HV_LRH,HV_Negative,HV_RN_UnV]  -
             _NewVaccinations[s,HV_LRH,t] *
             value2 +
            (_Adults[s,HV_LRH,HV_Negative,HV_RN_Take,t] +
             _Adults[s,HV_LRH,HV_Negative,HV_RN_Partial,t] +
             _Adults[s,HV_LRH,HV_Negative,HV_RN_NoProt,t]) /
             _VaccineEffectiveness[HV_RN_Duration]);
    _Adults[s,HV_LRH,HV_Negative,HV_RN_UnV,t] := value;
  end; {s loop}

  for s := HV_Male to HV_Female do
  begin

    if s = HV_Male then
    begin
  //        r1 := HV_MRH_M;
      r2 := HV_HRH;
      r3 := HV_IDU;
  //    r4 := HV_LRH_M;
    end
    else
    begin
  //        r1 := HV_MRH_F3;
      r2 := HV_HRH_F3;
      r3 := HV_IDU_F3;
  //    r4 := HV_LRH_F3;
    end;

    for r := HV_MRH to HV_MSMIDU do
    begin
      if not((s = HV_Female) and (r >= HV_MSM)) then
      begin
        if s = HV_Male then nr := r
        else nr := r + HV_Total;
        value2 := _Adults[s,r,HV_AllHIV,HV_RN_UnV,t];
        if value2 <> 0 then
          value2 :=  _Adults[s,r,HV_Negative,HV_RN_UnV,t] /
                    _Adults[s,r,HV_AllHIV,HV_RN_UnV,t];

        if HV_GetRNVaccineTargeting(p) = 1 then
          value2 := 1;

        if nr <= HV_IDU_F3 then
        begin

          value := _Adults[s,r,HV_Negative,HV_RN_UnV,t] +
                     (1/HV_TimeStep)*(_AdultsTS[s,HV_None,HV_Negative,HV_RN_UnV,_ts] /  //use  _AdultsTS
                    (max(_AgeFirstSex[s,t] - 15,1)) *
                     _RiskGroupPercent[nr] -
                     (_NonAIDSDeathRate[s,t] +
                      _RateofAging[s,t]+
                     _BehaviorChangeRate[nr]) *
                     _Adults[s,r,HV_Negative,HV_RN_UnV,t] +
                     _Migrant[s,r,HV_Negative,HV_RN_UnV]  -
                     _NewVaccinations[s,r,t] *
                     value2 +
                    (_Adults[s,r,HV_Negative,HV_RN_Take,t] +
                     _Adults[s,r,HV_Negative,HV_RN_Partial,t] +
                     _Adults[s,r,HV_Negative,HV_RN_NoProt,t]) /
                     _VaccineEffectiveness[HV_RN_Duration]);

          if r = HV_MRH then
            value := _Adults[s,r,HV_Negative,HV_RN_UnV,t] +
                     (1/HV_TimeStep)*(_AdultsTS[s,HV_None,HV_Negative,HV_RN_UnV,_ts] /  //use  _AdultsTS
                    (max(_AgeFirstSex[s,t] - 15,1)) *
                     _RiskGroupPercent[nr] +
                    _Adults[s,HV_HRH,HV_Negative,HV_RN_UnV,t] *
                    _BehaviorChangeRate[r2] -
                     (_NonAIDSDeathRate[s,t] +
                      _RateofAging[s,t]+
                      _BehaviorChangeRate[nr]) *
                     _Adults[s,r,HV_Negative,HV_RN_UnV,t] +
                     _Migrant[s,r,HV_Negative,HV_RN_UnV]  -
                     _NewVaccinations[s,r,t] *
                     value2 +
                    (_Adults[s,r,HV_Negative,HV_RN_Take,t] +
                     _Adults[s,r,HV_Negative,HV_RN_Partial,t] +
                     _Adults[s,r,HV_Negative,HV_RN_NoProt,t]) /
                     _VaccineEffectiveness[HV_RN_Duration]);

            if r = HV_HRH then
              value := _Adults[s,r,HV_Negative,HV_RN_UnV,t] +
                         (1/HV_TimeStep)*(_AdultsTS[s,HV_None,HV_Negative,HV_RN_UnV,_ts] /  //use  _AdultsTS
                        (max(_AgeFirstSex[s,t] - 15,1)) *
                         _RiskGroupPercent[nr]  +
                        _Adults[s,HV_IDU,HV_Negative,HV_RN_UnV,t] *
                        _BehaviorChangeRate[r3] -
                         ( _NonAIDSDeathRate[s,t] +
                          _RateofAging[s,t]+
                         _BehaviorChangeRate[nr]) *
                         _Adults[s,r,HV_Negative,HV_RN_UnV,t] +
                         _Migrant[s,r,HV_Negative,HV_RN_UnV]  -
                         _NewVaccinations[s,r,t] *
                         value2 +
                        (_Adults[s,r,HV_Negative,HV_RN_Take,t] +
                         _Adults[s,r,HV_Negative,HV_RN_Partial,t] +
                         _Adults[s,r,HV_Negative,HV_RN_NoProt,t]) /
                         _VaccineEffectiveness[HV_RN_Duration]);

            if r = HV_IDU then
              value := _Adults[s,r,HV_Negative,HV_RN_UnV,t] +
                       (1/HV_TimeStep)*(_AdultsTS[s,HV_None,HV_Negative,HV_RN_UnV,_ts] /  //use  _AdultsTS
                      (max(_AgeFirstSex[s,t] - 15,1)) *
                       _RiskGroupPercent[nr] -
                       ( _NonAIDSDeathRate[s,t] +
                        _RateofAging[s,t]+
                       _BehaviorChangeRate[nr]) *
                       _Adults[s,r,HV_Negative,HV_RN_UnV,t] +
                       _Migrant[s,r,HV_Negative,HV_RN_UnV]  -
                       _NewVaccinations[s,r,t] *
                       value2 +
                      (_Adults[s,r,HV_Negative,HV_RN_Take,t] +
                       _Adults[s,r,HV_Negative,HV_RN_Partial,t] +
                       _Adults[s,r,HV_Negative,HV_RN_NoProt,t]) /
                       _VaccineEffectiveness[HV_RN_Duration]);


          _Adults[s,r,HV_Negative,HV_RN_UnV,t] := value;
        end; {nr < HV_HRH_F_SP}
      end; {not((s = HV_Female) and (r >= HV_MSM)}
    end; {r loop}
  end; {s loop}

  (* At this point _Adults[S,HV_None,HV_Negative,V,t] and _Adults[s,r,HV_Negative,HV_RN_UnV,t]
      is updated, so use _AdultsTS, which represents _Adults before the update, in certain references
      to _Adults below *)
  for s := HV_Male to HV_Female do
  begin
    for v := HV_RN_Take to HV_RN_NoProt do
    begin
      value2 := _Adults[s,HV_LRH,HV_AllHIV,HV_RN_UnV,t];
      if value2 <> 0 then
        value2 := _AdultsTS[s,HV_LRH,HV_Negative,HV_RN_UnV,_ts] /
               _AdultsTS[s,HV_LRH,HV_AllHIV,HV_RN_UnV,_ts]; {use AdultsTS}

      if HV_GetRNVaccineTargeting(p) = 1 then value2 := 1;

      if s = HV_Male then
      begin
        r1 := HV_MRH;
  //          r2 := HV_HRH_M;
  //          r3 := HV_IDU_M;
        r4 := HV_LRH;
      end
      else
      begin
        r1 := HV_MRH_F3;
  //          r2 := HV_HRH_F3;
  //          r3 := HV_IDU_F3;
        r4 := HV_LRH_F3;
      end;

      value := _Adults[s,HV_LRH,HV_Negative,v,t] +
               (1/HV_TimeStep)*(_AdultsTS[s,HV_None,HV_Negative,v,_ts]/  // use _AdultsTS
               (max(_AgeFirstSex[s,t] - 15,1)) *
               _RiskGroupPercent[r4] +
               _NewVaccinations[s,HV_LRH,t] *
               value2 *
               _Veffect[_TNType,v] +
               _Adults[s,HV_MRH,HV_Negative,v,t] *
               _BehaviorChangeRate[r1] -
               //
               {HV_HRH->HV_MRH, HV_IDU_HV_HRH, as is done for behavior transitions in HIV+ categories} {+
               //
               _Adults[s,HV_HRH,HV_Negative,v,t] /
               _Behavior[r2,HV_AvgDur] +
               _Adults[s,HV_IDU,HV_Negative,v,t] /
               _Behavior[r3,HV_AvgDur]}
              (_NonAIDSDeathRate[s,t] +
               _RateofAging[s,t]) *
               _Adults[s,HV_LRH,HV_Negative,v,t] +
               _Migrant[s,HV_LRH,HV_Negative,v]  -
               _AdultsTS[s,HV_LRH,HV_Negative,v,_ts] / {use AdultsTS}
               _VaccineEffectiveness[HV_RN_Duration]);

      _Adults[s,HV_LRH,HV_Negative,v,t] := value;
    end;{v loop}
  end;{s loop}

  for s := HV_Male to HV_Female do
  begin
    if s = HV_Male then
    begin
  //  r1 := HV_MRH_M;
      r2 := HV_HRH;
      r3 := HV_IDU;
  //  r4 := HV_LRH_M;
    end
    else
    begin
  //  r1 := HV_MRH_F3;
      r2 := HV_HRH_F3;
      r3 := HV_IDU_F3;
  //  r4 := HV_LRH_F3;
    end;

    for r := HV_MRH to HV_MSMIDU do
    begin

      for v := HV_RN_Take to HV_RN_NoProt do
      begin

        if not((s = HV_Female) and (r >= HV_MSM)) then
        begin
          if s = HV_Male then
            nr := r
          else
            nr := r + HV_Total;

          value2 := _Adults[s,r,HV_AllHIV,HV_RN_UnV,t];
          if value2 <> 0 then
            value2 := _Adults[s,r,HV_Negative,HV_RN_UnV,t] /
                 _Adults[s,r,HV_AllHIV,HV_RN_UnV,t];

          if HV_GetRNVaccineTargeting(p) = 1 then value2 := 1;

          if nr <= HV_IDU_F3 then
          begin

             value := _Adults[s,r,HV_Negative,v,t] +
                       (1/HV_TimeStep)*(_AdultsTS[s,HV_None,HV_Negative,v,_ts] / //use _AdultsTS
                       (max(_AgeFirstSex[s,t] - 15,1)) *
                       _RiskGroupPercent[nr] +
                       _NewVaccinations[s,r,t] *
                       value2 *
                       _Veffect[_TNType,v] -
                       (_NonAIDSDeathRate[s,t] +
                       _RateofAging[s,t]+
                       _BehaviorChangeRate[nr]) *
                       _Adults[s,r,HV_Negative,v,t] +
                       _Migrant[s,r,HV_Negative,v]  -
                       _Adults[s,r,HV_Negative,v,t] /
                       _VaccineEffectiveness[HV_RN_Duration]);


            if r = HV_MRH then
                value := _Adults[s,r,HV_Negative,v,t] +
                       (1/HV_TimeStep)*(_AdultsTS[s,HV_None,HV_Negative,v,_ts] / //use _AdultsTS
                       (max(_AgeFirstSex[s,t] - 15,1)) *
                       _RiskGroupPercent[nr] +
                       _Adults[s,HV_HRH,HV_Negative,v,t] *
                       _BehaviorChangeRate[r2] +
                       _NewVaccinations[s,r,t] *
                       value2 *
                       _Veffect[_TNType,v] -
                       (_NonAIDSDeathRate[s,t] +
                       _RateofAging[s,t]+
                       _BehaviorChangeRate[nr]) *
                       _Adults[s,r,HV_Negative,v,t] +
                       _Migrant[s,r,HV_Negative,v]  -
                       _Adults[s,r,HV_Negative,v,t] /
                       _VaccineEffectiveness[HV_RN_Duration]);

              if r = HV_HRH then
                value := _Adults[s,r,HV_Negative,v,t] +
                       (1/HV_TimeStep)*(_AdultsTS[s,HV_None,HV_Negative,v,_ts] / // use _AdultsTS
                       (max(_AgeFirstSex[s,t] - 15,1)) *
                       _RiskGroupPercent[nr] +
                       _Adults[s,HV_IDU,HV_Negative,v,t] *
                       _BehaviorChangeRate[r3] +
                       _NewVaccinations[s,r,t] *
                       value2 *
                       _Veffect[_TNType,v] -
                       (_NonAIDSDeathRate[s,t] +
                       _RateofAging[s,t]+
                       _BehaviorChangeRate[nr]) *
                       _Adults[s,r,HV_Negative,v,t] +
                       _Migrant[s,r,HV_Negative,v]  -
                       _Adults[s,r,HV_Negative,v,t] /
                       _VaccineEffectiveness[HV_RN_Duration]);

                  if r = HV_IDU then
                  value := _Adults[s,r,HV_Negative,v,t] +
                         (1/HV_TimeStep)*(_AdultsTS[s,HV_None,HV_Negative,v,_ts] / //_AdultsTS
                         (max(_AgeFirstSex[s,t] - 15,1)) *
                         _RiskGroupPercent[nr] +
                         _NewVaccinations[s,r,t] *
                         value2 *
                         _Veffect[_TNType,v] -
                         (_NonAIDSDeathRate[s,t] +
                         _RateofAging[s,t]+
                         _BehaviorChangeRate[nr]) *
                         _Adults[s,r,HV_Negative,v,t] +
                         _Migrant[s,r,HV_Negative,v]  -
                         _Adults[s,r,HV_Negative,v,t] /
                         _VaccineEffectiveness[HV_RN_Duration]);

            _Adults[s,r,HV_Negative,v,t] := value;

          end;{if nr <= HV_IDU_F3 then}
        end;{if not((s = HV_Female) and (r >= HV_MSM)) then}
      end;{v loop}
    end;{r loop}
  end;{s loop}
end;

procedure THVProj.Calc_ProgressHIVPrimary(t:integer);
var
  s,r,v,r1,r2,r3,VE   : byte;
  value,value2        : double;
begin
  {* Primary Infection *}
  for s := HV_Male to HV_Female do
  begin
    value2 := _Adults[s,HV_LRH,HV_AllHIV,HV_RN_UnV,t];
    if value2 <> 0 then
      value2 := _Adults[s,HV_LRH,HV_Primary,HV_RN_UnV,t] /
                _Adults[s,HV_LRH,HV_AllHIV,HV_RN_UnV,t];

    if s = HV_Male then
    begin
      r1 := HV_MRH;
      r2 := HV_HRH;
      r3 := HV_IDU;
    end
    else
    begin
      r1 := HV_MRH_F3;
      r2 := HV_HRH_F3;
      r3 := HV_IDU_F3;
    end;

    value := _Adults[s,HV_LRH,HV_Primary,HV_RN_UnV,t] +
             (1/HV_TimeStep)*(- _StageExits[s,HV_LRH,HV_Primary,HV_RN_UnV] +
              _Migrant[s,HV_LRH,HV_Primary,HV_RN_UnV] +
              _BehaviorChangeRate[r1] * _Adults[s,HV_MRH,HV_Primary,HV_RN_UnV,t] -
             _NewVaccinations[s,HV_LRH,t] * value2 +
             (_Adults[s,HV_LRH,HV_Primary,HV_RN_Take,t] +
             _Adults[s,HV_LRH,HV_Primary,HV_RN_Partial,t] +
             _Adults[s,HV_LRH,HV_Primary,HV_RN_NoProt,t]) /
             _VaccineEffectiveness[HV_RN_Duration]);

    _Adults[s,HV_LRH,HV_Primary,HV_RN_UnV,t] := value;

    for r := HV_MRH to HV_IDU do
    begin
      value2 := _Adults[s,r,HV_AllHIV,HV_RN_UnV,t];
      if value2 <> 0 then
        value2 := _Adults[s,r,HV_Primary,HV_RN_UnV,t] /
                   _Adults[s,r,HV_AllHIV,HV_RN_UnV,t];

      if r = HV_MRH then
        value := _Adults[s,r,HV_Primary,HV_RN_UnV,t] +
               (1/HV_TimeStep)*(- _StageExits[s,r,HV_Primary,HV_RN_UnV] +
                _Migrant[s,r,HV_Primary,HV_RN_UnV] +
                _BehaviorChangeRate[r2] * _Adults[s,HV_HRH,HV_Primary,HV_RN_UnV,t] -
               _NewVaccinations[s,r,t] *value2 +
               (_Adults[s,r,HV_Primary,HV_RN_Take,t] +
               _Adults[s,r,HV_Primary,HV_RN_Partial,t] +
               _Adults[s,r,HV_Primary,HV_RN_NoProt,t]) /
               _VaccineEffectiveness[HV_RN_Duration]);
      if r = HV_HRH then
        value := _Adults[s,r,HV_Primary,HV_RN_UnV,t] +
               (1/HV_TimeStep)*(- _StageExits[s,r,HV_Primary,HV_RN_UnV] +
                _Migrant[s,r,HV_Primary,HV_RN_UnV] +
                _BehaviorChangeRate[r3] * _Adults[s,HV_IDU,HV_Primary,HV_RN_UnV,t] -
                _NewVaccinations[s,r,t] * value2 +
               (_Adults[s,r,HV_Primary,HV_RN_Take,t] +
               _Adults[s,r,HV_Primary,HV_RN_Partial,t] +
               _Adults[s,r,HV_Primary,HV_RN_NoProt,t]) /
               _VaccineEffectiveness[HV_RN_Duration]);
      if r = HV_IDU then
        value := _Adults[s,r,HV_Primary,HV_RN_UnV,t] +
               (1/HV_TimeStep)*(- _StageExits[s,r,HV_Primary,HV_RN_UnV] +
                _Migrant[s,r,HV_Primary,HV_RN_UnV] -
               _NewVaccinations[s,r,t] * value2 +
               (_Adults[s,r,HV_Primary,HV_RN_Take,t] +
               _Adults[s,r,HV_Primary,HV_RN_Partial,t] +
               _Adults[s,r,HV_Primary,HV_RN_NoProt,t]) /
               _VaccineEffectiveness[HV_RN_Duration]);

      _Adults[s,r,HV_Primary,HV_RN_UnV,t] := value;
    end;{r loop}

    (* This check is needed to ensure the calcs are not repeated when s=HV_Female *)
    if s = HV_Male then
    begin
      for r := HV_MSM to HV_MSMIDU do
      begin
        value2 := _Adults[s,r,HV_AllHIV,HV_RN_UnV,t];
        if value2 <> 0 then
          value2 := _Adults[HV_Male,r,HV_Primary,HV_RN_UnV,t] /
                 _Adults[s,r,HV_AllHIV,HV_RN_UnV,t];

        value := _Adults[HV_Male,r,HV_Primary,HV_RN_UnV,t] +
                 (1/HV_TimeStep)*(- _StageExits[HV_Male,r,HV_Primary,HV_RN_UnV] +
                  _Migrant[HV_Male,r,HV_Primary,HV_RN_UnV] -
                 _NewVaccinations[HV_Male,r,t] * value2 +
                 (_Adults[HV_Male,r,HV_Primary,HV_RN_Take,t] +
                 _Adults[HV_Male,r,HV_Primary,HV_RN_Partial,t] +
                 _Adults[HV_Male,r,HV_Primary,HV_RN_NoProt,t]) /
                 _VaccineEffectiveness[HV_RN_Duration]);

        _Adults[HV_Male,r,HV_Primary,HV_RN_UnV,t] := value;
      end;{r loop}
    end;{if s = HV_Male, for MSM}

    for v := HV_RN_Take to HV_RN_NoProt do
    begin
      value2 := _Adults[s,HV_LRH,HV_AllHIV,HV_RN_UnV,t];
      if value2 <> 0 then
        value2 := _Adults[s,HV_LRH,HV_Primary,HV_RN_UnV,t] /
               _Adults[s,HV_LRH,HV_AllHIV,HV_RN_UnV,t];

      if HV_GetRNVaccineTargeting(p) = 1 then
        value2 := 0;

      if v = HV_RN_NoProt then
        VE := 1
      else
        VE := 0;

      value := _Adults[s,HV_LRH,HV_Primary,v,t] +
               (1/HV_TimeStep)*(- _StageExits[s,HV_LRH,HV_Primary,v] +
                _Migrant[s,HV_LRH,HV_Primary,v] +
               _NewVaccinations[s,HV_LRH,t] * value2 * VE +
                _BehaviorChangeRate[r1] * _Adults[s,HV_MRH,HV_Primary,v,t] -
               _Adults[s,HV_LRH,HV_Primary,v,t] / _VaccineEffectiveness[HV_RN_Duration]);

      _Adults[s,HV_LRH,HV_Primary,v,t] := value;

      for r := HV_MRH to HV_IDU do
      begin
        value2 := _Adults[s,r,HV_AllHIV,HV_RN_UnV,t];
        if value2 <> 0 then
          value2 := _Adults[s,r,HV_Primary,HV_RN_UnV,t] /
                 _Adults[s,r,HV_AllHIV,HV_RN_UnV,t];

        if v = HV_RN_NoProt then
          VE := 1
        else
          VE := 0;

        if r = HV_MRH then
          value := _Adults[s,r,HV_Primary,v,t] +
                 (1/HV_TimeStep)*(- _StageExits[s,r,HV_Primary,v] +
                  _Migrant[s,r,HV_Primary,v] +
                 _BehaviorChangeRate[r2] * _Adults[s,HV_HRH,HV_Primary,v,t] +
                _NewVaccinations[s,r,t] * value2 * VE -
                _Adults[s,r,HV_Primary,v,t] / _VaccineEffectiveness[HV_RN_Duration]);
        if r = HV_HRH then
          value := _Adults[s,r,HV_Primary,v,t] +
                 (1/HV_TimeStep)*(- _StageExits[s,r,HV_Primary,v] +
                  _Migrant[s,r,HV_Primary,v] +
                   _BehaviorChangeRate[r3] * _Adults[s,HV_IDU,HV_Primary,v,t] +
                  _NewVaccinations[s,r,t] * value2 * VE -
                  _Adults[s,r,HV_Primary,v,t] / _VaccineEffectiveness[HV_RN_Duration]);
        if r = HV_IDU then
          value := _Adults[s,r,HV_Primary,v,t] +
                (1/HV_TimeStep)*(- _StageExits[s,r,HV_Primary,v] +
                 _Migrant[s,r,HV_Primary,v] +
                _NewVaccinations[s,r,t] * value2 * VE -
                _Adults[s,r,HV_Primary,v,t] / _VaccineEffectiveness[HV_RN_Duration]);

        _Adults[s,r,HV_Primary,v,t] := value;
      end; {r loop}

     (* This check is needed to ensure the calcs are not repeated when s=HV_Female *)
     if s = HV_Male then
     begin
      for r := HV_MSM to HV_MSMIDU do
      begin
        value2 := _Adults[HV_Male,r,HV_AllHIV,HV_RN_UnV,t];
        if value2 <> 0 then
          value2 := _Adults[HV_Male,r,HV_Primary,HV_RN_UnV,t] /
                 _Adults[HV_Male,r,HV_AllHIV,HV_RN_UnV,t];

        value := _Adults[HV_Male,r,HV_Primary,v,t] +
                 (1/HV_TimeStep)*(- _StageExits[HV_Male,r,HV_Primary,v] +
                  _Migrant[HV_Male,r,HV_Primary,v] +
                 _NewVaccinations[HV_Male,r,t] * value2 * VE -
                 _Adults[HV_Male,r,HV_Primary,v,t] / _VaccineEffectiveness[HV_RN_Duration]);

        _Adults[HV_Male,r,HV_Primary,v,t] := value;
        end; {r loop}
      end;  {if s = HV_Male, for MSM}
    end;  {v loop}
  end; {s loop}

  for s := HV_Male to HV_Female do
    for r := HV_None to HV_MSMIDU do
      for v := HV_RN_UnV to HV_RN_NoProt do
        _Adults[s,r,HV_Primary,v,t]:=max(0,_Adults[s,r,HV_Primary,v,t]);
end;

procedure THVProj.Calc_ProgressHIVpos(t,_ts:integer);
var
  s,r,h,v,r1,r2,r3,r4,VE  : byte;
  value,value2            : double;
begin


   (* Asymptomatic and Symptomatic HIV *)
   (* Loop over CD4 stages *)
  for h := HV_CD4_GT500 to HV_CD4_LT50 do
  begin
    for s := HV_Male to HV_Female do
    begin
      value2 := _Adults[s,HV_LRH,HV_AllHIV,HV_RN_UnV,t];
      if value2 <> 0 then
        value2 := _Adults[s,HV_LRH,h,HV_RN_UnV,t] /
                  _Adults[s,HV_LRH,HV_AllHIV,HV_RN_UnV,t];

      if HV_GetRNVaccineTargeting(p) = 1 then value2 := 0;

      if s = HV_Male then
      begin
        r1 := HV_MRH;
        r2 := HV_HRH;
        r3 := HV_IDU;
        r4 := HV_LRH;
      end
      else
      begin
        r1 := HV_MRH_F3;
        r2 := HV_HRH_F3;
        r3 := HV_IDU_F3;
        r4 := HV_LRH_F3;
      end;

      value := _Adults[s,HV_LRH,h,HV_RN_UnV,t] +
               (1/HV_TimeStep)*( _StageProgressors[s,HV_LRH,h-1,HV_RN_UnV] -
                _StageExits[s,HV_LRH,h,HV_RN_UnV] +
                _Migrant[s,HV_LRH,h,HV_RN_UnV] +
                _BehaviorChangeRate[r1] * _Adults[s,HV_MRH,h,HV_RN_UnV,t] -
               _NewVaccinations[s,HV_LRH,t] * value2 +
               (_Adults[s,HV_LRH,h,HV_RN_Take,t] +
               _Adults[s,HV_LRH,h,HV_RN_Partial,t] +
               _Adults[s,HV_LRH,h,HV_RN_NoProt,t]) /
               _VaccineEffectiveness[HV_RN_Duration]+
               _AdultsTS[s,HV_None,h,HV_RN_UnV,_ts]/(max(_AgeFirstSex[s,t] - 15,1)) *
              _RiskGroupPercent[r4]);

      _Adults[s,HV_LRH,h,HV_RN_UnV,t] := value;

      for r := HV_MRH to HV_IDU do
      begin
        value2 := _Adults[s,r,HV_AllHIV,HV_RN_UnV,t];
        if value2 <> 0 then
          value2 := _Adults[s,r,h,HV_RN_UnV,t] /
                _Adults[s,r,HV_AllHIV,HV_RN_UnV,t];

        if HV_GetRNVaccineTargeting(p) = 1 then
          value2 := 0;

        if r = HV_MRH then
          value := _Adults[s,r,h,HV_RN_UnV,t] +
                  (1/HV_TimeStep)*(_StageProgressors[s,r,h-1,HV_RN_UnV] -
                  _StageExits[s,r,h,HV_RN_UnV] +
                  _Migrant[s,r,h,HV_RN_UnV] +
                   _BehaviorChangeRate[r2] * _Adults[s,HV_HRH,h,HV_RN_UnV,t] -
                  _NewVaccinations[s,r,t] *  value2 +
                  (_Adults[s,r,h,HV_RN_Take,t] +
                  _Adults[s,r,h,HV_RN_Partial,t] +
                  _Adults[s,r,h,HV_RN_NoProt,t]) /
                  _VaccineEffectiveness[HV_RN_Duration]+
                  _AdultsTS[s,HV_None,h,HV_RN_UnV,_ts]/(max(_AgeFirstSex[s,t] - 15,1)) *
                  _RiskGroupPercent[r1]);
        if r = HV_HRH then
          value := _Adults[s,r,h,HV_RN_UnV,t] +
                  (1/HV_TimeStep)*(_StageProgressors[s,r,h-1,HV_RN_UnV] -
                  _StageExits[s,r,h,HV_RN_UnV] +
                  _Migrant[s,r,h,HV_RN_UnV] +
                   _BehaviorChangeRate[r3] * _Adults[s,HV_IDU,h,HV_RN_UnV,t] -
                  _NewVaccinations[s,r,t] *  value2 +
                  (_Adults[s,r,h,HV_RN_Take,t] +
                  _Adults[s,r,h,HV_RN_Partial,t] +
                  _Adults[s,r,h,HV_RN_NoProt,t]) /
                  _VaccineEffectiveness[HV_RN_Duration]++
                  _AdultsTS[s,HV_None,h,HV_RN_UnV,_ts]/(max(_AgeFirstSex[s,t] - 15,1)) *
                  _RiskGroupPercent[r2]);
        if r = HV_IDU then
          value := _Adults[s,r,h,HV_RN_UnV,t] +
                  (1/HV_TimeStep)*(_StageProgressors[s,r,h-1,HV_RN_UnV] -
                   _StageExits[s,r,h,HV_RN_UnV] +
                   _Migrant[s,r,h,HV_RN_UnV] -
                  _NewVaccinations[s,r,t] *  value2 +
                  (_Adults[s,r,h,HV_RN_Take,t] +
                  _Adults[s,r,h,HV_RN_Partial,t] +
                  _Adults[s,r,h,HV_RN_NoProt,t]) /
                  _VaccineEffectiveness[HV_RN_Duration]+
                  _AdultsTS[s,HV_None,h,HV_RN_UnV,_ts]/(max(_AgeFirstSex[s,t] - 15,1)) *
                  _RiskGroupPercent[r3]);

        _Adults[s,r,h,HV_RN_UnV,t] := value;
      end;

      if(s = HV_Male) then
      begin
        for r := HV_MSM to HV_MSMIDU do
        begin
          value2 := _Adults[HV_Male,r,HV_AllHIV,HV_RN_UnV,t];
          if value2 <> 0 then
            value2 := _Adults[HV_Male,r,h,HV_RN_UnV,t] /
                   _Adults[HV_Male,r,HV_AllHIV,HV_RN_UnV,t];

         if HV_GetRNVaccineTargeting(p) = 1 then
            value2 := 0;

          value := _Adults[HV_Male,r,h,HV_RN_UnV,t] +
                   (1/HV_TimeStep)*(_StageProgressors[HV_Male,r,h-1,HV_RN_UnV] -
                   _StageExits[HV_Male,r,h,HV_RN_UnV] +
                   _Migrant[HV_Male,r,h,HV_RN_UnV] -
                   _NewVaccinations[HV_Male,r,t] * value2 +
                   (_Adults[HV_Male,r,h,HV_RN_Take,t] +
                   _Adults[HV_Male,r,h,HV_RN_Partial,t] +
                   _Adults[HV_Male,r,h,HV_RN_NoProt,t]) /
                   _VaccineEffectiveness[HV_RN_Duration]+
                    _AdultsTS[s,HV_None,h,HV_RN_UnV,_ts]/(max(_AgeFirstSex[s,t] - 15,1)) *
                    _RiskGroupPercent[r]);

          _Adults[HV_Male,r,h,HV_RN_UnV,t] := value;
        end;
      end;{if s=HV_Male}

      for v := HV_RN_Take to HV_RN_NoProt do
      begin
        value2 := _Adults[s,HV_LRH,HV_AllHIV,HV_RN_UnV,t];
        if value2 <> 0 then
          value2 := _Adults[s,HV_LRH,h,HV_RN_UnV,t] /
                 _Adults[s,HV_LRH,HV_AllHIV,HV_RN_UnV,t];

        if HV_GetRNVaccineTargeting(p) = 1 then
          value2 := 0;

        if v = HV_RN_NoProt then
          VE := 1
        else
          VE := 0;

        value := _Adults[s,HV_LRH,h,v,t] +
                 (1/HV_TimeStep)*(_StageProgressors[s,HV_LRH,h-1,v] -
                  _StageExits[s,HV_LRH,h,v] +
                  _Migrant[s,HV_LRH,h,v] +
                  _NewVaccinations[s,HV_LRH,t] *  value2 * VE +
                 _BehaviorChangeRate[r1] * _Adults[s,HV_MRH,h,v,t] +
  //                 1 / _Behavior[r2,HV_AvgDur] *
  //                 _Adults[s,HV_HRH,h,v,t] +
  //                 1 / _Behavior[r3,HV_AvgDur] *
  //                 _Adults[s,HV_IDU,h,v,t] -
                 - _Adults[s,HV_LRH,h,v,t] /_VaccineEffectiveness[HV_RN_Duration]);

        _Adults[s,HV_LRH,h,v,t] := value;

        for r := HV_MRH to HV_IDU do
        begin
          value2 := _Adults[s,r,HV_AllHIV,HV_RN_UnV,t];
          if value2 <> 0 then
            value2 := _Adults[s,r,h,HV_RN_UnV,t] /
                   _Adults[s,r,HV_AllHIV,HV_RN_UnV,t];

          if HV_GetRNVaccineTargeting(p) = 1 then
            value2 := 0;

          if v = HV_RN_NoProt then
            VE := 1
          else
            VE := 0;

          value := _Adults[s,r,h,v,t] +
                   1/HV_TimeStep*(_StageProgressors[s,r,h-1,v] -
                    _StageExits[s,r,h,v] +
                    _Migrant[s,r,h,v] +
                   _NewVaccinations[s,r,t] * value2 * VE -
                   _Adults[s,r,h,v,t] /_VaccineEffectiveness[HV_RN_Duration]);

          _Adults[s,r,h,v,t] := value;
        end;

        if (s = HV_Male) then
        begin
          for r := HV_MSM to HV_MSMIDU do
          begin
            value2 := _Adults[HV_Male,r,HV_AllHIV,HV_RN_UnV,t];
            if value2 <> 0 then
              value2 := _Adults[HV_Male,r,h,HV_RN_UnV,t] /
                     _Adults[HV_Male,r,HV_AllHIV,HV_RN_UnV,t];

            if HV_GetRNVaccineTargeting(p) = 1 then
              value2 := 0;

            if v = HV_RN_NoProt then
              VE := 1
            else
              VE := 0;

            value := _Adults[HV_Male,r,h,v,t]+
                     (1/HV_TimeStep)*(_StageProgressors[HV_Male,r,h-1,v] -
                      _StageExits[HV_Male,r,h,v] +
                      _Migrant[HV_Male,r,h,v] +
                     _NewVaccinations[HV_Male,r,t] *value2 * VE  -
                     _Adults[HV_Male,r,h,v,t] /_VaccineEffectiveness[HV_RN_Duration]);

            _Adults[HV_Male,r,h,v,t] := value;
          end; {r loop}
        end; {if s=HV_Male}
      end; {v loop}
    end;{ s loop}
  end;{h loop : outer CD4 loop}

  for s := HV_Male to HV_Female do
    for r := HV_None to HV_MSMIDU do
      for h := HV_CD4_GT500 to HV_CD4_LT50 do
        for v := HV_RN_UnV to HV_RN_NoProt do
          _Adults[s,r,h,v,t]:=max(0,_Adults[s,r,h,v,t]);


   (* Calculate progression below current art eligibility threshold *)
  for s := HV_Male to HV_Female do
    for h := HV_CD4_GT500 to HV_CD4_LT50 do
    begin
      if ( (_CD4_lower_limit[h] < _CD4ThreshHoldAdults) and (h=_h_cov) and (_h_cov<>HV_Primary) )  then
      begin
        for r := HV_None to HV_MSMIDU do
          for v := HV_RN_UnV to HV_RN_NoProt do //note offset on h for _NewlyEligibleART
            _NewlyEligibleART[s,h+10,t]:= _NewlyEligibleART[s,h+10,t]+(1/HV_TimeStep)*_StageProgressors[s,r,h-1,v];
      end;{if}
    end;{s,h loop}
end;

procedure THVProj.Calc_ProgressART(t,_ts:integer);
var
  s,r,h,v,r1,r2,r3,r4     : byte;
  value                   : double;
  bAddEntrants            : integer;
begin
  (* Loop over ART CD4 stages *)
  for h := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
  begin
    (* _NewInfectionsT[HV_RN_Take]=0, so will never receive ART *)
    for v := HV_RN_UnV to HV_RN_NoProt do
    begin
      if(v=HV_RN_UnV) then
        bAddEntrants:=1
      else
        bAddEntrants:=0;

      for s := HV_Male to HV_Female do
      begin
        (* [CDP] note offsett on h for ExitrateART
           Note behavior change for ART receiver *)
        if s = HV_Male then
        begin
          r1 := HV_MRH;
          r2 := HV_HRH;
          r3 := HV_IDU;
          r4 := HV_LRH;
        end
        else
        begin
          r1 := HV_MRH_F3;
          r2 := HV_HRH_F3;
          r3 := HV_IDU_F3;
          r4 := HV_LRH_F3;
        end;

        value := _Adults[s,HV_LRH,h,v,t] +
                 (1/HV_TimeStep)*(- _StageExits[s,HV_LRH,h,v] + _Migrant[s,HV_LRH,h,v]+
                  _BehaviorChangeRate[r1] * _Adults[s,HV_MRH,h,v,t]+
                 bAddEntrants*_AdultsTS[s,HV_None,h,HV_RN_UnV,_ts]/(max(_AgeFirstSex[s,t] - 15,1)) *
                 _RiskGroupPercent[r4]);

        _Adults[s,HV_LRH,h,v,t] := value;

        for r := HV_MRH to HV_IDU do
        begin
          if r = HV_MRH then
            value := _Adults[s,r,h,v,t] +
                    (1/HV_TimeStep)*(- _StageExits[s,r,h,v] + _Migrant[s,r,h,v]+
                     _BehaviorChangeRate[r2] * _Adults[s,HV_HRH,h,v,t]+
                     bAddEntrants*_AdultsTS[s,HV_None,h,HV_RN_UnV,_ts]/(max(_AgeFirstSex[s,t] - 15,1)) *
                    _RiskGroupPercent[r1]);
          if r = HV_HRH then
            value := _Adults[s,r,h,v,t] +
                    (1/HV_TimeStep)*(-_StageExits[s,r,h,v] + _Migrant[s,r,h,v] +
                     _BehaviorChangeRate[r3] * _Adults[s,HV_IDU,h,v,t]+
                     bAddEntrants*_AdultsTS[s,HV_None,h,HV_RN_UnV,_ts]/(max(_AgeFirstSex[s,t] - 15,1)) *
                    _RiskGroupPercent[r2]);
          if r = HV_IDU then
            value := _Adults[s,r,h,v,t] +
                    (1/HV_TimeStep)*(- _StageExits[s,r,h,v] + _Migrant[s,r,h,v]+
                     bAddEntrants*_AdultsTS[s,HV_None,h,HV_RN_UnV,_ts]/(max(_AgeFirstSex[s,t] - 15,1)) *
                    _RiskGroupPercent[r3]);

           _Adults[s,r,h,v,t] := value;
        end;

        if(s=HV_Male) then
        begin
          for r := HV_MSM to HV_MSMIDU do
          begin
            value := _Adults[HV_Male,r,h,v,t] +
                     (1/HV_TimeStep)*(- _StageExits[HV_Male,r,h,v] + _Migrant[HV_Male,r,h,v] +
                      bAddEntrants*_AdultsTS[s,HV_None,h,HV_RN_UnV,_ts]/(max(_AgeFirstSex[s,t] - 15,1)) *
                      _RiskGroupPercent[r]);

            _Adults[HV_Male,r,h,v,t] := value;
          end;{r loop}
        end;//{if s=HV_Male}
      end;{s loop}
    end;{v loop}
  end;{h loop : outer CD4 loop}

  for s := HV_Male to HV_Female do
    for r := HV_None to HV_MSMIDU do
      for h := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
        for v := HV_RN_UnV to HV_RN_NoProt do
          _Adults[s,r,h,v,t]:=max(0,_Adults[s,r,h,v,t]);
end;

procedure THVProj.Init_GetFirstYearART;
var
  t : integer;
  h : byte;
begin
  _ARTYrIdx := GetGBCalcYearIdx(p);
  (* Number or percent used *)
  if( (HV_GetDPARTCoverageSelection(p) = HV_DP_NumOrPercent)) then
  begin
    for t := GetGBCalcYearIdx(p) to HV.GetFinalIndex do
        if( ((HV_GetDPHAARTBySex(p,t,DP_Male) > 0) or (HV_GetDPHAARTBySex(p,t,DP_Female) > 0))and (_ARTYrIdx = 1)) then
            _ARTYrIdx:=t;
  end//CD4 coverage used
  else if(HV_GetDPARTCoverageSelection(p)= HV_DP_CD4Percent + 1) then
  begin
    for t := GetGBCalcYearIdx(p) to HV.GetFinalIndex do
       for h := HV_CD4_GT500 to HV_CD4_LT50 do
           if((HV_GetDPCD4Coverage(p,t,h-1)>0) and (_ARTYrIdx=1)) then
              _ARTYrIdx:=t;
  end;
end;

procedure THVProj.Calc_DistributeNewRecruits(t:integer);
var
  s,r,nr        : byte;
  value,value2  : double;
begin
  for s := HV_Male to HV_Female do
  begin
    value2 := 0;
    if _Adults[s,HV_AllRisk,HV_AllHIV,HV_RN_AllVacc,t] <> 0 then
    begin
      for r := HV_LRH to HV_MSMIDU do
      begin
        if s = HV_Male then
          nr := r
        else
          nr := r + HV_Total;
        if nr <= HV_IDU_F3 then
        begin
        (* allow negative adjustment to reduce relative risk group size *)
        if (_Behavior[nr,HV_PercPop]/100 -
        (_Adults[s,r,HV_AllHIV,HV_RN_AllVacc,t]/
         _Adults[s,HV_AllRisk,HV_AllHIV,HV_RN_AllVacc,t])) <> 0 then
          begin
            if (_Adults[s,HV_AllRisk,HV_AllHIV,HV_RN_AllVacc,t] = 0) or
               (_Behavior[nr,HV_AvgDur] = 0) or
               (_Adults[s,HV_None,HV_Negative,HV_RN_UnV,t] / {use AdultsTS}
                 (max(_AgeFirstSex[s,t] - 15,1)) = 0)
            then value := 0
            else
              value := (_Behavior[nr,HV_PercPop]/100) +
                       ( ((_Behavior[nr,HV_PercPop]/100) -
                       (1-_BehaviorChangeRate[nr])*(_Adults[s,r,HV_AllHIV,HV_RN_AllVacc,t]/
                       (_Adults[s,HV_AllRisk,HV_AllHIV,HV_RN_AllVacc,t])) ) *
                       (_IncRecruitment[s,r]/100)+
                       (0*(_Behavior[nr,HV_PercPop]/100)/
                         _Behavior[nr,HV_AvgDur]))*
                       ( (_Adults[s,HV_AllRisk,HV_AllHIV,HV_RN_AllVacc,t])/
                       (_Adults[s,HV_None,HV_Negative,HV_RN_UnV,t] / {use AdultsTS}
                       (max(_AgeFirstSex[s,t] - 15,1))));

            (* no movement back to not sexually active group *)
            _RiskGroupPercent[nr] := max(value,0);
          end;
          value2 := value2 + _RiskGroupPercent[nr];
        end; {nr <= HV_HRH_F_S}
      end; {r loop}

      {* Adjust for MSM risk groups not being used *}
      if s = HV_Male then
      begin
        if (HV.GetNumMSMRiskGroups = 1) then {using four groups}
        begin
          value2 := value2 - _RiskGroupPercent[HV_MSM];
          _RiskGroupPercent[HV_MSM] := 0
        end
        else {using 1 group}
        begin
          for r := HV_MSMLR to HV_MSMIDU do
          begin
            value2 := value2 - _RiskGroupPercent[r];
            _RiskGroupPercent[r] := 0;
          end;
        end;
      end; {if s=HV_Male}

      {* Rescale overall adjustment *}
      if s = HV_Male then
      begin
        for r := HV_LRH to HV_MSMIDU do
        begin
        value:=0;
        if(value2>0) then
                value := _RiskGroupPercent[r]/value2;
          _RiskGroupPercent[r] := value;
        end; {r loop}
      end
      else
      begin
        for nr := HV_LRH_F3 to HV_IDU_F3 do
        begin
          value := 0;
          if(value2>0) then
                  value := _RiskGroupPercent[nr]/value2;
          _RiskGroupPercent[nr] := value;
        end;
      end; {if s=HV_Male}
    end; {adult <> 0}
  end; {s loop}
end;

procedure THVProj.Calc_HIVMortTotals(t:integer);
var
  s,r,v     : byte;
  value     : double;
begin
  (* Sum AIDS Mortality array over vaccination state *)
  for s := HV_Male to HV_Female do
  begin
    for r := HV_None to HV_MSMIDU do
    begin
        value := 0;
        for v := HV_RN_UnV to HV_RN_NoProt do
          value := value + _AIDSDeaths[s,r,v,t];

        _AIDSDeaths[s,r,HV_RN_AllVacc,t] := value;
    end;{r loop}
  end;{s loop}

 (* Sum AIDS Mortality ART array over Risk group *)
  for s := HV_Male to HV_Female do
  begin
    value:=0;
    for r := HV_None to HV_MSMIDU do
        value := value + _AIDSDeathsART[s,r,t];

    _AIDSDeathsART[s,HV_AllRisk,t] := value;
  end;{s loop}

  for r := HV_None to HV_MSMIDU do
      _AIDSDeathsART[HV_Bothsexes,r,t]:=_AIDSDeathsART[HV_Male,r,t]+_AIDSDeathsART[HV_Female,r,t];

  _AIDSDeathsART[HV_Bothsexes,HV_AllRisk,t]:=_AIDSDeathsART[HV_Male,HV_AllRisk,t]+_AIDSDeathsART[HV_Female,HV_AllRisk,t];
end;

procedure THVProj.Calc_HIVPopTotals(t:integer);
var
  s,r,h     : byte;
  value     : double;
begin
  (* _TotalAdultsHIV *)
  for s := HV_Male to HV_Female do
  begin
    value := 0;
    for r := HV_None to HV_MSMIDU do
    begin
      for h := HV_Primary to HV_CD4_LT50_ART do
      begin
        if not((s = HV_Female) and (r >= HV_MSM)) then
        begin
          value := value + _Adults[s,r,h,HV_RN_UnV,t] +
                   _Adults[s,r,h,HV_RN_Take,t] +
                   _Adults[s,r,h,HV_RN_Partial,t] +
                   _Adults[s,r,h,HV_RN_NoProt,t];
        end;
      end;
    end;
    _TotalAdultsHIV[s,t] := value;//[CDP] not on ART
  end;

  (* _TotalAdultsHIV on ART *)
  for s := HV_Male to HV_Female do
  begin
    value := 0;
    for r := HV_None to HV_MSMIDU do
    begin
      for h := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
      begin
        if not((s = HV_Female) and (r >= HV_MSM)) then
        begin
          value := value + _Adults[s,r,h,HV_RN_UnV,t] +
                   _Adults[s,r,h,HV_RN_Take,t] +
                   _Adults[s,r,h,HV_RN_Partial,t] +
                   _Adults[s,r,h,HV_RN_NoProt,t];
        end;
      end;
    end;
    _TotalAdultsART[s,t] := value;
  end;
end;

procedure THVProj.Calc_TotalPops(t,t1:integer);
var
  s,r,h,v       : byte;
  value,value2  : double;
begin
  for s := HV_Male to HV_Female do
  begin
    value2 := 0;
    for r := HV_None to HV_MSMIDU do
    begin
      value := 0;
      for h := HV_Negative to HV_CD4_LT50_ART do
        for v := HV_RN_UnV to HV_RN_NoProt do
          value := value + _Adults[s,r,h,v,t];

       value2 := value2 + value;
      _TotalPops[s,r,t1,t] := value;
    end;{r loop}
    _TotalPops[s,HV_AllRisk,t1,t] := value2;
  end;{s loop}

  for r := HV_None to HV_MSMIDU do
    _TotalPops[HV_Bothsexes,r,t1,t] := _TotalPops[HV_Male,r,t1,t] + _TotalPops[HV_Female,r,t1,t];

 _TotalPops[HV_Bothsexes,HV_AllRisk,t1,t] := _TotalPops[HV_Male,HV_AllRisk,t1,t] + _TotalPops[HV_Female,HV_AllRisk,t1,t];
end;

procedure THVProj.Calc_TotalEntrants(t,_ts:integer);
var
  s,r,h,v,nr       :  byte;
  active,migrants,
  value,value2     :  double;

begin
  for s := HV_Male to HV_Female do
  begin
    value:=(1/HV_TimeStep) * (_Pop15rate[s,HV_DP_H_No,t] + _Migrant[s,HV_None,HV_Negative,HV_RN_UnV]);

    (* not at risk,hiv+ *)
    for h := HV_CD4_GT500 to HV_CD4_LT50 do
      value := value+(1/HV_TimeStep) * (_Pop15rate[s,h,t] + _Migrant[s,HV_None,h,HV_RN_UnV]);

    (* not at risk,hiv+ *)
    for h := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
      value:=value+(1/HV_TimeStep) * (_Pop15rate[s,h,t] + _Migrant[s,HV_None,h,HV_RN_UnV]);

    _TotalEntrants[s,HV_None,1,t] := _TotalEntrants[s,HV_None,1,t]+value;
    value2 := value;

    for r := HV_LRH to HV_MSMIDU do
    begin
      if not((s = HV_Female) and (r >= HV_MSM)) then
      begin
        (* at risk HIV- *)
        migrants:=0;
        for v := HV_RN_UnV to HV_RN_NoProt do
          migrants := migrants + _Migrant[s,r,HV_Negative,v];

        if s = HV_Male then
        nr := r
        else
        nr := r + HV_Total;

        active:=0;
        for v := HV_RN_UnV to HV_RN_NoProt do
        active := active + _AdultsTS[s,HV_None,HV_Negative,v,_ts]/(max(_AgeFirstSex[s,t] - 15,1)) *
          _RiskGroupPercent[nr];

        for h := HV_CD4_GT500 to HV_CD4_LT50_ART do
        begin
          for v := HV_RN_UnV to HV_RN_NoProt do
            migrants := migrants + _Migrant[s,r,h,v];

          if s = HV_Male then
            nr := r
          else
            nr := r + HV_Total;

          for v := HV_RN_UnV to HV_RN_NoProt do
            active := active + _AdultsTS[s,HV_None,h,v,_ts]/(max(_AgeFirstSex[s,t] - 15,1)) *
              _RiskGroupPercent[nr];
        end;

        value := (1/HV_TimeStep) * (active + migrants);

        _TotalEntrants[s,r,1,t] := _TotalEntrants[s,r,1,t] + value;
        value2 := value2 + value;

        (* HIV+ not on ART *)
        for h := HV_Primary to HV_CD4_LT50 do
        begin
           migrants:=0;
           for v := HV_RN_UnV to HV_RN_NoProt do
              migrants := migrants + _Migrant[s,r,h,v];
           value := (1/HV_TimeStep) * (migrants);
        end;
        _TotalEntrants[s,r,1,t] := _TotalEntrants[s,r,1,t] + value;
        value2 := value2 + value;

        (* HIV+ on ART *)
        for h := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
        begin
           migrants:=0;
           for v := HV_RN_UnV to HV_RN_NoProt do
              migrants := migrants + _Migrant[s,r,h,v];
           value := (1/HV_TimeStep) * (migrants);
        end;

        _TotalEntrants[s,r,1,t] := _TotalEntrants[s,r,1,t] + value;
        value2 := value2 + value;

      end;{if s}
    end;{r loop}

    _TotalEntrants[s,HV_AllRisk,1,t] := _TotalEntrants[s,HV_AllRisk,1,t] + value2;
  end;{s loop}

  for r := HV_None to HV_MSMIDU do
   _TotalEntrants[HV_Bothsexes,r,1,t] := _TotalEntrants[HV_Male,r,1,t] + _TotalEntrants[HV_Female,r,1,t];

  _TotalEntrants[HV_Bothsexes,HV_AllRisk,1,t] := _TotalEntrants[HV_Male,HV_AllRisk,1,t] + _TotalEntrants[HV_Female,HV_AllRisk,1,t];

  (* total new migrants *)
  for s := HV_Male to HV_Female do
    for r := HV_None to HV_MSMIDU do
      for h := HV_Negative to HV_CD4_LT50 do
        for v := HV_RN_UnV to HV_RN_NoProt do
           if not((s = HV_Female) and (r >= HV_MSM)) then
              _TotalMigrants[t] := _TotalMigrants[t] + (1/HV_TimeStep) * _Migrant[s,r,h,v];

  for s := HV_Male to HV_Female do
  begin
    (* total new 15 year olds *)
    _TotalNewAdults[t] := _TotalNewAdults[t] + (1/HV_TimeStep) * _Pop15rate[s, HV_DP_H_No, t];

    for h := HV_CD4_GT500 to HV_CD4_LT50 do
    begin
         _TotalNewAdults[t] := _TotalNewAdults[t]+(1/HV_TimeStep) * _Pop15rate[s, h, t];
         _TotalNewAdults[t] := _TotalNewAdults[t]+(1/HV_TimeStep) * _Pop15rate[s, h+10, t];
    end;{h loop}
  end;
end;

procedure THVProj.Calc_TotalExits(t,_ts:integer);
var
  s,r,h,v         : byte;
  value,value2    : double;
begin
  for s := HV_Male to HV_Female do
  begin
  value2:=0;

  (* not at risk *)
  value:=0;
  for v := HV_RN_UnV to HV_RN_NoProt do
    value:=value+(1/HV_TimeStep)*(_NonAIDSDeathRate[s,t]+_RateofAging[s,t]+1/(max(_AgeFirstSex[s,t] - 15,1)))*_AdultsTS[s,HV_None,HV_Negative,v,_ts];

  (* not at risk,hiv+ *)
  for h := HV_CD4_GT500 to HV_CD4_LT50 do
    value:=value+(1/HV_TimeStep)*(_NonAIDSDeathRate[s,t]+_mu[h,s]*_ScaleMuART[h,s]+_RateofHIVAging[s,h,t]+1/(max(_AgeFirstSex[s,t] - 15,1)))*_AdultsTS[s,HV_None,h,HV_RN_UnV,_ts];

   (* not at risk on ART *)
  for h := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
    value:=value+(1/HV_TimeStep)*(_NonAIDSDeathRate[s,t]+_alpha[h,s]+_RateofARTAging[s,h,t]+1/(max(_AgeFirstSex[s,t] - 15,1)))*_AdultsTS[s,HV_None,h,HV_RN_UnV,_ts];

  _TotalExits[s,HV_None,1,t]:=_TotalExits[s,HV_None,1,t]+value;
  value2:=value2+value;

  for r := HV_LRH to HV_MSMIDU do
  begin
     if not((s = HV_Female) and (r >= HV_MSM)) then
      begin
       (* at risk HIV- *)
        value:=0;
        for v := HV_RN_UnV to HV_RN_NoProt do
          value:=value+(1/HV_TimeStep)*(_NonAIDSDeathRate[s,t]+_RateofAging[s,t])*_AdultsTS[s,r,HV_Negative,v,_ts];
        _TotalExits[s,r,1,t]:=_TotalExits[s,r,1,t]+value;

        value2:=value2+value;
      end; {if}
  end;{r loop}

  for r := HV_LRH to HV_MSMIDU do
  begin
     if not((s = HV_Female) and (r >= HV_MSM)) then
      begin
        (* HIV+ not on ART *)
        value:=0;
        for h := HV_Primary to HV_CD4_LT50 do
        for v := HV_RN_UnV to HV_RN_NoProt do
           value:=value+(1/HV_TimeStep)*(_NonAIDSDeathRate[s,t]+_mu[h,s]*_ScaleMuART[h,s]+_RateofHIVAging[s,h,t])*_AdultsTS[s,r,h,v,_ts];
        _TotalExits[s,r,1,t]:=_TotalExits[s,r,1,t]+value;
        value2:=value2+value;

        (* HIV+ on ART *)
        value:=0;
        for h := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
          for v := HV_RN_UnV to HV_RN_NoProt do
            value:=value+(1/HV_TimeStep)*(_NonAIDSDeathRate[s,t]+_alpha[h,s]+_RateofARTAging[s,h,t])*_AdultsTS[s,r,h,v,_ts];
        _TotalExits[s,r,1,t]:=_TotalExits[s,r,1,t]+value;
        value2:=value2+value;
      end;{if s}
    end;{r loop}
      _TotalExits[s,HV_AllRisk,1,t]:=_TotalExits[s,HV_AllRisk,1,t]+value2;
  end;{s loop}

  for r := HV_None to HV_MSMIDU do
    _TotalExits[HV_Bothsexes,r,1,t]:=_TotalExits[HV_Male,r,1,t]+_TotalExits[HV_Female,r,1,t];

  _TotalExits[HV_Bothsexes,HV_AllRisk,1,t]:=_TotalExits[HV_Male,HV_AllRisk,1,t]+_TotalExits[HV_Female,HV_AllRisk,1,t];
end;


procedure THVProj.PreInnerLoopCalcs(t: integer);
var
  s,r,h,v : byte;
begin

   (* Update behavioral parameters through impact matrix *)
   if(t > HV.GetGoalsBaseYearIdx) then
      Calc_BehavioralCoverage(t);

  (* Set CD4-based ART eligibility threshold *)
  _CD4ThreshHoldAdults := HV_GetDPCD4ThreshHoldAdults(p,t);
  Calc_SetARTEligibility(t);

  (* CD4 transition parameters *)
  Calc_CD4Trans(t);

  (* calculate non-AIDS mortality and aging rates *)
  Calc_NonAIDSMortAIM(t);

  (* Calculate new 15 year olds *)
  Calc_NewEntrantsAIM(t);

  (* calc migration rate and indicate its done *)
  Calc_MigrNum(t);

  (* Initialize all groups to zero change *)
  for s := HV_Bothsexes to HV_Female do
    for r := HV_AllRisk to HV_MSMIDU do
      for h := HV_Negative to HV_AllHIV do
          for v := HV_RN_AllVacc to HV_RN_NoProt do
              _Adults[s,r,h,v,t]:=_Adults[s,r,h,v,t-1];

  (* SumPopGroups here when t=2, else at the end of t loop *)
  if t = 2 then
    Calc_SumPopGroups(t);

  Calc_AdditionalARTCov(t);   //Calculate additional coverage for risk group
end;


function THVProj.GetRiskGrpARTCov(s,r,h,v,t:byte): double;
 begin
  Result := _RiskGrpARTCov[s,r,h,v,t];
 end;

 function THVProj.GetDP_RiskGrpARTPerc(s,h,t:byte): double;
 begin
  Result := _DP_RiskGrpARTPerc[s,h,t];
 end;

procedure THVProj.CalcNumOnART(t, t1: integer; var numOnART : GB_TMaleFemale);
var
  s, h,  ts: Integer;
begin
  ts := HV_TimeStep*t + t1;

  //HV_BothSexes..HV_Female,HV_AllRisk..HV_MSMIDU,HV_Negative..HV_AllHIV,HV_RN_AllVacc..HV_RN_NoProt
  for s := GB_Male to GB_Female do
  begin
    numOnART[s] := 0;
    for h := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
      numOnART[s] := numOnART[s] + _AdultsTS[s, HV_AllRisk, h, HV_RN_AllVacc, ts];
  end;
end;

procedure  THVProj.Calc_AdditionalARTCov(t: integer);
var
  s,r,h,v,
  found,i,
  specPop              : byte;
  TotalWomen,
  TotalNumART,
  TotalNumPop,
  TimeBtwPreg,
  dummy                : double;
  Filename             : string;
begin
  for s := HV_Male to HV_Female do
    for r := HV_None to HV_MSMIDU do
      _AddlRiskGrpARTCov[s,r,t] := 0;

  for r :=  DP_EligTreatPregnantWomen to DP_EligTreatPopsMax do
    if not(GetDPData(p).getPopsEligTreat(r).Eligible) or (GetDPData(p).getPopsEligTreat(r).Year > (GetGBCalcYear(p) + t - GetGBCalcYearIdx(p))) then
       GetDPData(p).setCovPopsEligTreat(t,r,0);

  //FSW
  _AddlRiskGrpARTCov[HV_Female, HV_HRH, t] := max(_AddlRiskGrpARTCov[HV_Female, HV_HRH, t], GetDPData(p).getCovPopsEligTreat(t,DP_EligTreatSexWorkers));
  //MSM
  _AddlRiskGrpARTCov[HV_Male, HV_MSM, t] := max(_AddlRiskGrpARTCov[HV_Male, HV_MSM, t], GetDPData(p).getCovPopsEligTreat(t,DP_EligTreatMSM));
  _AddlRiskGrpARTCov[HV_Male, HV_MSMIDU, t] := max(_AddlRiskGrpARTCov[HV_Male, HV_MSMIDU, t], GetDPData(p).getCovPopsEligTreat(t,DP_EligTreatMSM));
  //IDU
  _AddlRiskGrpARTCov[HV_Male, HV_MSMIDU, t] := max(_AddlRiskGrpARTCov[HV_Male, HV_MSMIDU, t], GetDPData(p).getCovPopsEligTreat(t,DP_EligTreatIDU));
  _AddlRiskGrpARTCov[HV_Female, HV_IDU, t] := max(_AddlRiskGrpARTCov[HV_Female, HV_IDU, t], GetDPData(p).getCovPopsEligTreat(t,DP_EligTreatIDU));
  _AddlRiskGrpARTCov[HV_Male, HV_IDU, t] := max(_AddlRiskGrpARTCov[HV_Male, HV_IDU, t], GetDPData(p).getCovPopsEligTreat(t,DP_EligTreatIDU));

  TimeBtwPreg := 1;
//  TimeBtwPreg := HV_GetBirthInterval(p);

   //if pregnant women
  TotalWomen := 0;
  if (GetDPData(p).getCovPopsEligTreat(t,DP_EligTreatPregnantWomen) > 0) then
    for h := HV_CD4_GT500 to HV_CD4_LT50 do
      for r := HV_None to HV_MSMIDU do
        for v := HV_RN_UnV to HV_RN_NoProt do
          TotalWomen :=   TotalWomen + _Adults[HV_Female,r,h,v,t] + _Adults[HV_Female,r,h+10,v,t] ;

  for s := HV_male to HV_Female do
    for r := HV_None to HV_MSMIDU do
      for h := HV_CD4_GT500 to HV_CD4_LT50 do
        for v := HV_RN_AllVacc to HV_RN_NoProt do
        begin
         (*Coverage of all except TB, SDC, and other pop*)
            _RiskGrpARTCov[s,r,h,v,t] := max(0, (_Adults[s,r,h,v,t] + _Adults[s,r,h+10,v,t])
                                              * (100 - GetDPData(p).getPopsEligTreat(DP_EligTreatTB_HIV).PercentHIV
                                                 - GetDPData(p).getPopsEligTreat(DP_EligTreatDiscordantCouples).PercentHIV
                                                 - GetDPData(p).getPopsEligTreat(DP_EligTreatOtherPop).PercentHIV)/ 100
                                              * max(_AddlRiskGrpARTCov[s,r,t],HV_GetDPCD4Coverage(p,t,h))/ 100);
            (*Add coverage of  TB*)
            SpecPop := DP_EligTreatTB_HIV;
            _RiskGrpARTCov[s,r,h,v,t] := _RiskGrpARTCov[s,r,h,v,t]
                                        + max(0, (_Adults[s,r,h,v,t] + _Adults[s,r,h+10,v,t])
                                              * (GetDPData(p).getPopsEligTreat(SpecPop).PercentHIV/ 100)
                                              * max(GetDPData(p).getCovPopsEligTreat(t,SpecPop)/100, max(_AddlRiskGrpARTCov[s,r,t],HV_GetDPCD4Coverage(p,t,h))/100 ));
            (*Add coverage of  SDC*)
            SpecPop := DP_EligTreatDiscordantCouples;
            _RiskGrpARTCov[s,r,h,v,t] := _RiskGrpARTCov[s,r,h,v,t]
                                         + max(0, (_Adults[s,r,h,v,t] + _Adults[s,r,h+10,v,t])
                                              * (GetDPData(p).getPopsEligTreat(SpecPop).PercentHIV/ 100)
                                              * max(GetDPData(p).getCovPopsEligTreat(t,SpecPop)/100, max(_AddlRiskGrpARTCov[s,r,t],HV_GetDPCD4Coverage(p,t,h))/100 ));

            (*Add coverage of 'Other population'*)
            SpecPop := DP_EligTreatOtherPop;
            _RiskGrpARTCov[s,r,h,v,t] := _RiskGrpARTCov[s,r,h,v,t]
                                         + max(0, (_Adults[s,r,h,v,t] + _Adults[s,r,h+10,v,t])
                                              * (GetDPData(p).getPopsEligTreat(SpecPop).PercentHIV/ 100)
                                              * max(GetDPData(p).getCovPopsEligTreat(t,SpecPop)/100, max(_AddlRiskGrpARTCov[s,r,t],HV_GetDPCD4Coverage(p,t,h))/100 ));

            (*Add coverage of PW*)
            (*Coverage so far*)
            dummy := 0;
            if (_Adults[HV_Female,r,h,v,t] + _Adults[HV_Female,r,h+10,v,t]) > 0 then
              dummy :=  _RiskGrpARTCov[HV_female,r,h,v,t]/ (_Adults[HV_Female,r,h,v,t] + _Adults[HV_Female,r,h+10,v,t]);
            if (GetDPData(p).getCovPopsEligTreat(t,DP_EligTreatPregnantWomen) > 0) and (s = HV_Female) then
              (* Note: GetDPNeedPMTCT   uses t-1 instead of t because at this point this indicator has not yet been estimated*)
              _RiskGrpARTCov[HV_female,r,h,v,t] := _RiskGrpARTCov[HV_female,r,h,v,t]
                                                 + max(0, (GetDPData(p).getNeedPMTCT(t-1) * (_Adults[HV_Female,r,h,v,t] + _Adults[HV_Female,r,h+10,v,t] )/ TotalWomen)
                                                         * TimeBtwPreg * max(GetDPData(p).getCovPopsEligTreat(t,DP_EligTreatPregnantWomen)/100 - dummy ,0));
          end;

  //For use in AIM
  for s := HV_male to HV_Female do
    for h := HV_CD4_GT500 to HV_CD4_LT50 do
    begin
      TotalNumART := 0;
      TotalNumPop := 0;
      for r := HV_None to HV_MSMIDU do
        for v := HV_RN_AllVacc to HV_RN_NoProt do
        begin
          TotalNumART := TotalNumART + _RiskGrpARTCov[s,r,h,v,t];
          TotalNumPop := TotalNumPop + _Adults[s,r,h,v,t] + _Adults[s,r,h+10,v,t];
        end;
      if TotalNumPop > 0 then
        _DP_RiskGrpARTPerc[s,h,t] :=  TotalNumART/ TotalNumPop
      else
        _DP_RiskGrpARTPerc[s,h,t] := 0;
    end;

  //For use in GOALS
  for s := HV_male to HV_Female do
    for h := HV_CD4_GT500 to HV_CD4_LT50 do
      for r := HV_None to HV_MSMIDU do
        for v := HV_RN_AllVacc to HV_RN_NoProt do
          if (_Adults[s,r,h,v,t] + _Adults[s,r,h+10,v,t]) > 0 then
            _HV_RiskGrpARTPerc[s,r,h,v,t] :=  _RiskGrpARTCov[s,r,h,v,t]/ (_Adults[s,r,h,v,t] + _Adults[s,r,h+10,v,t] )
          else
            _HV_RiskGrpARTPerc[s,r,h,v,t] := 0;

//       for h := HV_CD4_GT500 to HV_CD4_350_500 do
//       if (_RiskGrpARTPerc[HV_male,h,t] > 0) or (_RiskGrpARTPerc[HV_female,h,t] > 0) then
//       HV_SetDPCD4Coverage(p,t,h, max(_RiskGrpARTPerc[HV_male,h,t], _RiskGrpARTPerc[HV_Female,h,t]))  ;
end;

//function THVProj.HV_GetBirthInterval(p : byte) : double;
//var
//  CountryCode,
//  row                  : integer;
//  sheet                : GB_TFile;
//  found                : boolean;
//const
//  CountryCodeCol = 0;
//  BirthIntervalCol = 6;
//begin
//  (* set the found to false, when the country is located, found will be set to true *)
//  found := false;
//
//  (* Get the code of the country used in the projection *)
//  CountryCode := GetGBCountryRegionCode(p);
//
//  (* Create a sheet to read in the prevention technologies file *)
//  sheet := GB_TFile.Create;
//
//  try
//    sheet.ExtractSheet(GB_GetModDataFQName(GB_RN), HV_RN_PreventionTechnologies);
//
//    row := 2;
//    (* Read in rows until we find out which row matches the Country Code. *)
//    while ((not found) and (Trim(sheet.cells[CountryCodeCol, row]) <> '')) do
//    begin
//      if StrToInt(Trim(sheet.cells[CountryCodeCol, row])) = CountryCode then
//        found := true
//      else
//        Inc(row);
//    end;
//
//    (* if a match was found, read in the Epidemic type, PrEP coverage,
//       proportion discordant couples, and proportion first births *)
//    if found then
//      Result := sheet.floats[BirthIntervalCol,row]
//    else
//      Result := 1;
//
//  finally
//    sheet.Free;
//  end;
//end;

procedure THVProj.InnerLoopCalcs(t,t1 : integer);
var
  s,r,h,v : byte;
begin
  (* Keep track of the inner time index *)
  _ts:=_ts+1;

  (* Need a copy of the _Adults array before any updates to it *)
  for s := HV_Bothsexes to HV_Female do
    for r := HV_AllRisk to HV_MSMIDU do
      for h := HV_Negative to HV_AllHIV do
        for v := HV_RN_AllVacc to HV_RN_NoProt do
          _AdultsTS[s,r,h,v,_ts]:=_Adults[s,r,h,v,t];

  (* Sum to get total pops, needed in migration rate, vaccination and incidence calcs *)
  Calc_PopulationTotals(t);

  (* Calc Migration rate, turn to number for s,r,h,v *)
   Calc_MigrRate(t);

  (* Aging Rate *)
  if(t1=1) then
  begin
    Calc_HIVNegAgingAIM(t,_ts);
    if bHAsHIV or GetGBInCalcStateMode(p) then
      Calc_HIVPosAgingAIM(t,_ts);
  end;

  (* Get CD4 parameters from AIM, used in exit rates *)
  //if ((t1 = 1) and (bHAsHIV)) then
   //  Calc_CD4TransAIM(t,t1,_ts);

  (* calc HIV exit and progression rate *)
  Calc_HIVExitRates(t);
  (* apply rates and get exits and progressors by s,r,h,v *)
  Calc_HIVExitsProgressors(t);

  (* Get ART coverage from AM *)
 //Calc_ARTCoverage(t);//falls away

  {***Vaccinations***}

  (* Total Vaccinations and Un-Vaccinated *)
  Calc_VaccinationTotals(t);
  (* New vaccinations to maintain coverage *)
  if t1 = 1 then
    Calc_NewVaccinations(t);

  {***Not at risk population***}

  (* Update not at risk population *)
  (* HIV-, not at risk, to themselves or others *)
  Calc_ProgressNoRiskHIVNeg(t,_ts);

  (* HIV+, not at risk others *)
  Calc_ProgressNoRiskHIVPos(t);


  {***At risk population***}

  (* Get HIV prevalence *)
  Calc_Prevalence(t,_ts);
  (* Calc new infections if HIV already introduced *)
  if bHasHIV or GetGBInCalcStateMode(p) then
    bHasNewInfections := Calc_NewInfections(t);

  (* At this point AIM could progress for one t1 time step using NewInfections/HV_TimeStep
     Goals can adjust progression rates based on aging, progression and mortlaity in AIM  *)

  (* Calculate distribution of new recruits *)
  if t >= 2 then
    Calc_DistributeNewRecruits(t);

  (* Progress HIV- at risk population and distribute new recruits *)
  Calc_ProgressHIVneg(t,_ts);

  {***HIV+: Primary, Not on ART and on ART***}

  (* Progress HIV+:Primary, Not on ART *)
  if bHasHIV or GetGBInCalcStateMode(p) then
  begin
    Calc_ProgressHIVPrimary(t);
    Calc_ProgressHIVpos(t,_ts);
  end;

  (* Progress HIV+ on ART, ART enrollment,after first year of specified ART coverage/numbers *)
  if (bHasHIV or GetGBInCalcStateMode(p)) and (t >= _ARTYrIdx) then
  begin
    Calc_ProgressART(t,_ts);
    Calc_AllocateART(t);
  end;

  {***Add new infections and initial pulse***}
  _NewHIV[HV_Male]:=0;_NewHIV[HV_Female]:=0;
  (* Distribute New infections *)
  if bHasNewInfections then
    Calc_AddNewInfections(t,HV_TimeStep);

  (* Distribute initial Pulse *)
//  if GetGBFinalYear(p) >= (GetGBCalcYear(p) + _EpidemicStYr - GetGBCalcYearIdx(p)) then
  if (GetGBCalcYear(p) <= _EpidemicStYr) and (GetGBFinalYear(p) >= _EpidemicStYr) then
//    if ((t = _EpidemicStYr-1) and (t1 = HV_TimeStep)) then
    if ((t = _EpidemicStYr - GetGBCalcYear(p) + 1) and (t1 = HV_TimeStep)) then
    begin
        Calc_InitialPulse(t,_ts);
        (* if HIV not already added to IDU through sexual mixing *)
        if (_Adults[HV_Male,HV_IDU,HV_Primary,HV_RN_UnV,t]=0) then
          Calc_InitialPulseIDU(t,_ts);
    end;

  (* New infection accumulate in _NewHIV, in: Calc_AddNewInfections, Calc_InitialPulse or Calc_InitialPulseIDU*)
  if( ((_NewHIV[HV_Male]>0) or (_NewHIV[HV_Female]>0)) ) then
      Calc_AddNewInfectionsAIM(t);

  (* Check if HIV+ individuals were introduced*)
  bHasHIV:=Calc_HasHIV(t);

  (* SumPopGroups *)
  Calc_SumPopGroups(t);

  (* Track indicators at time scale t1 *)
  Calc_TrackAIDSDeaths(t,_ts);
  Calc_TrackNewInfections(t);

  (* Stockflow Info *)
  if (bDoStockFlow = True) then
  begin
    (* Keep track of total exits to occur in time t. Accumulate over t1. includes aging rates. *)
    Calc_TotalExits(t,_ts);
    (* Keep track of new entrants during time t. Accumulate over t1. includes migrants. *)
    Calc_TotalEntrants(t,_ts);
  end;
end;

procedure THVProj.PostInnerLoopCalcs(t: integer);
begin
  (* record total pop after prior to t1 loop for 'stockflow' *)
  if (bDoStockFlow = True) then
    Calc_TotalPops(t,2);

  (* Update _AIDSDeaths and _AIDSDeathsART and sum totals *)
  Calc_HIVMortTotals(t);

  (* Update _AdultsHIV and _AdultsART, and sum totals *)
  Calc_HIVPopTotals(t);

  Calc_IncSexRatio(t);

end;

procedure THVProj.PostOuterLoopCalcs;
begin
  Calc_Incidence;
  Calc_PercentPop;
  Calc_ARTCoverageDisplay;

  if bDoStockFlow = True then
    Calc_WriteStockFlow;
end;

procedure THVProj.Calc_WriteStockFlow;
var
  value,
  value1 : double;
  s,r,t  : byte;
  col,
  row    : integer;
  sheet  : GB_TFile;
begin
  sheet := GB_TFile.Create;
  try
    row := 0;
    col := 0;
    for t := 2 to HV.GetFinalIndex do
    begin
      sheet.Ints[col,row] := t;
      inc(col);
      sheet.Floats[col,row] := _TotalNewAdults[t];
      inc(col);
      sheet.Floats[col,row] := _TotalMigrants[t];
      inc(col);
      sheet.Floats[col,row] := -1;
      inc(col);

      for s := HV_Bothsexes to HV_Bothsexes do
        for r := HV_AllRisk to HV_MSMIDU do
        begin
          sheet.Floats[col,row] := _TotalPops[s,r,1,t];
          inc(col);

          value1 := _TotalPops[s,r,2,t];
          sheet.Floats[col,row] := value1;
          inc(col);

          value := _TotalPops[s,r,1,t]+
                   _TotalEntrants[s,r,1,t]-
                   _TotalExits[s,r,1,t];
          sheet.Floats[col,row] := value;
          inc(col);

          if(value > 0) then
             sheet.Floats[col,row] := (value - value1) / value * 100
          else
            sheet.Floats[col,row] := 0;
          inc(col);

          sheet.Floats[col,row] := -1;
          inc(col);
        end;{s,r loop}
      inc(row);
    end;
    sheet.saveToCSV(GetGBAbrvProjFName(GetGBActiveProj) + '_StockFlow' + GB_CSV_DOTEXT);
  finally
    sheet.free;
  end;
end;

procedure THVProj.Init_IDUBehavior;
var
  t: integer;
  value: Double;
begin
  for t := GetGBCalcYearIdx(p) to GetGBFinalYearIdx(p) do
  begin
    (* IDU: condom use *)
    value := HV.GetCondomPercent(HV_MRH,t);
    HV.SetCondomPercent(HV_IDU,t,value);
    (* IDU: number of partners *)
    value := HV.GetNumPart(HV_MRH,t);
    HV.SetNumPart(HV_IDU,t,value);
    value := HV.GetNumPart(HV_MRH_F3,t);
    HV.SetNumPart(HV_IDU_F3,t,value);
    (* IDU: sex acts *)
    value := HV.GetSexActs(HV_MRH,t);
    HV.SetSexActs(HV_IDU,t,value);
    value := HV.GetSexActs(HV_MRH_F3,t);
    HV.SetSexActs(HV_IDU_F3,t,value);

    (* MSMIDU: condom use *)
    value := HV.GetCondomPercent(HV_MSMMR,t);
    HV.SetCondomPercent(HV_MSMIDU,t,value);
    (* MSMIDU: number of partners *)
    value := HV.GetNumPart(HV_MSMMR,t);
    HV.SetNumPart(HV_MSMIDU,t,value);
    (* MSMIDU: sex acts *)
    value := HV.GetSexActs(HV_MSMMR,t);
    HV.SetSexActs(HV_MSMIDU,t,value);
 end;
end;

procedure THVProj.Calc_BehavioralCoverage(t:integer);
var
  i,r,s,
  CurYrIdx    : Integer;
  Value       : Double;
begin
  (* Calc adjusted coverage *)
  CurYrIdx := HV.GetGoalsBaseYearIdx;
  for i := 1 to HV_RN_IDUDrugSub do
    _AdjCoverage[i, CurYrIdx] := HV_GetRNCoverage(p, CurYrIdx, i)/100;

  //for t := GetHVGoalsCurYearIdx(p) + 1 to GetHVFinalIndex(p) do
  //begin
    for i := 1 to HV_RN_IDUDrugSub do
    begin
      case i of
        HV_RN_SecondaryTeachers : begin

         (* Adjust coverage of school-based intervention for proportion of adults
             that are in school *)
          _AdjCoverage[i, t] := (HV_GetRNCoverage(p, t, i)/100 -
              _AdjCoverage[i, CurYrIdx]) *
              ((HV_GetRNPopulationSizes(p,t,HV_RN_POP_SecondarySchoolMale) +
                HV_GetRNPopulationSizes(p,t,HV_RN_POP_SecondarySchoolFemale)) / 2) / 100 *
              HV_DP_Get_Pop_Age(p, t-1, HV_DP_BothSexes,DP_A15_19, DP_A15_19) /
              HV_DP_Get_Pop_Age(p, t-1, HV_DP_BothSexes,DP_A15_19, DP_A45_49)
        end;
        HV_RN_OutOfSchoolYouth : begin
         (* Adjust coverage of youth-based interventions for proportion of adults
            that are youth *)
          _AdjCoverage[i, t] := (HV_GetRNCoverage(p, t, i)/100 -
              _AdjCoverage[i, CurYrIdx]) * HV_DP_Get_Pop_Age(p, t-1, HV_DP_BothSexes,DP_A15_19, DP_A15_19) /
              HV_DP_Get_Pop_Age(p, t-1, HV_DP_BothSexes,DP_A15_19, DP_A45_49)
        end;
        HV_RN_CashTransfers : begin
          (* Adjust coverage of workplace interventions for proportion of adults
             that are in the formal sector labor force *)
          _AdjCoverage[i, t] := (HV_GetRNCoverage(p, t, i)/100 -
              _AdjCoverage[i, CurYrIdx]) *
//              (HV_GetRNPopulationSizes(p,GetHVGoalsBaseYearIdx(p),HV_RN_POP_FormalSectorEmployees)/
//               HV_DP_Get_Pop_Age(p, GetHVGoalsBaseYearIdx(p), HV_DP_BothSexes,DP_A15_19, DP_A45_49));
              HV_DP_Get_Pop_Age(p, t-1, HV_DP_Female,DP_A15_19, DP_A20_24) /
              HV_DP_Get_Pop_Age(p, t-1, HV_DP_BothSexes,DP_A15_19, DP_A45_49);

        end;
        else
          _AdjCoverage[i, t] := HV_GetRNCoverage(p, t, i)/100 -
              _AdjCoverage[i, CurYrIdx];
      end;
    end;
  //end;{t loop}
  (* Calc adjusted behaviors *)

  (* Condom use *)
  //for t := CurYrIdx + 1 to GetHVFinalIndex(p) do
  //begin
   (* Low risk *)
    value := 1;
    for i := 1 to HV_RN_IDUDrugSub do
      if not (i in HV_RN_NotInCEAorIM) then
        value := value * (1 + _AdjCoverage[i, t] * HV.GetImpactMatrix(HV_Data, i, HV_RGReducCondomNonUse_Low)/100);
    value := 1 - (1-HV.GetCondomPercent(HV_LRH,CurYrIdx)) * value;
    value := max(value,0);
    value := min(value,1);
    HV.SetCondomPercent(HV_LRH,t, Value);

    (* Med Risk *)
    value := 1;
    for i := 1 to HV_RN_IDUDrugSub do
      if not (i in HV_RN_NotInCEAorIM) then
        value := value * (1 + _AdjCoverage[i, t] * HV.GetImpactMatrix(HV_Data, i, HV_RGReducCondomNonUse_Med)/100);
    value := 1 - (1-HV.GetCondomPercent(HV_MRH,CurYrIdx)) * value;
    value := max(value,0);
    value := min(value,1);
    HV.SetCondomPercent(HV_MRH,t, Value);

    (* High Risk *)
    value := 1;
    for i := 1 to HV_RN_IDUDrugSub do
      if not (i in HV_RN_NotInCEAorIM) then
        if (i <> HV_RN_MSMOutreach) and (i <> HV_RN_MSMLub) then
          value := value * (1 + _AdjCoverage[i, t] * HV.GetImpactMatrix(HV_Data, i, HV_RGReducCondomNonUse_High)/100);
    value := 1 - (1-HV.GetCondomPercent(HV_HRH,CurYrIdx)) * value;
    value := max(value,0);
    value := min(value,1);
    HV.SetCondomPercent(HV_HRH,t, Value);

    (* IDU *)
    value := 1;
    for i := 1 to HV_RN_IDUDrugSub do
      if not (i in HV_RN_NotInCEAorIM) then
        if (i <> HV_RN_MSMOutreach) and (i <> HV_RN_MSMLub) then
          value := value * (1 + _AdjCoverage[i, t] * HV.GetImpactMatrix(HV_Data, i, HV_RGReducCondomNonUse_IDU)/100);
    value := 1 - (1-HV.GetCondomPercent(HV_IDU,CurYrIdx)) * value;
    value := max(value,0);
    value := min(value,1);
    HV.SetCondomPercent(HV_IDU,t, Value);

    (* MSM Risk *)
    value := 1;
    for i := 1 to HV_RN_IDUDrugSub do
      if not (i in HV_RN_NotInCEAorIM) then
        value := value * (1 + _AdjCoverage[i, t] * HV.GetImpactMatrix(HV_Data, i, HV_RGReducCondomNonUse_MSM)/100);
    value := 1 - (1-HV.GetCondomPercent(HV_MSM,CurYrIdx)) * value;
    value := max(value,0);
    value := min(value,1);
    HV.SetCondomPercent(HV_MSM,t, Value);

    if HV.GetNumMSMRiskGroups > 0 then
    begin
      (* Low risk *)
      value := 1;
      for i := 1 to HV_RN_IDUDrugSub do
        if not (i in HV_RN_NotInCEAorIM) then
          value := value * (1 + _AdjCoverage[i, t] * HV.GetImpactMatrix(HV_Data, i, HV_RGReducCondomNonUse_MSM_Low)/100);
      value := 1 - (1-HV.GetCondomPercent(HV_MSMLR,CurYrIdx)) * value;
      value := max(value,0);
      value := min(value,1);
      HV.SetCondomPercent(HV_MSMLR,t, Value);

      (* Med Risk *)
      value := 1;
      for i := 1 to HV_RN_IDUDrugSub do
        if not (i in HV_RN_NotInCEAorIM) then
          value := value * (1 + _AdjCoverage[i, t] * HV.GetImpactMatrix(HV_Data, i, HV_RGReducCondomNonUse_MSM_Med)/100);
      value := 1 - (1-HV.GetCondomPercent(HV_MSMMR,CurYrIdx)) * value;
      value := max(value,0);
      value := min(value,1);
      HV.SetCondomPercent(HV_MSMMR,t, Value);

      (* High Risk *)
      value := 1;
      for i := 1 to HV_RN_IDUDrugSub do
        if not (i in HV_RN_NotInCEAorIM) then
          value := value * (1 + _AdjCoverage[i, t] * HV.GetImpactMatrix(HV_Data, i, HV_RGReducCondomNonUse_MSM_High)/100);
      value := 1 - (1-HV.GetCondomPercent(HV_MSMHR,CurYrIdx)) * value;
      value := max(value,0);
      value := min(value,1);
      HV.SetCondomPercent(HV_MSMHR,t, Value);

      (* MSM: IDU *)
      value := 1;
      for i := 1 to HV_RN_IDUDrugSub do
        if not (i in HV_RN_NotInCEAorIM) then
          value := value * (1 + _AdjCoverage[i, t] * HV.GetImpactMatrix(HV_Data, i, HV_RGReducCondomNonUse_MSM_IDU)/100);
      value := 1 - (1-HV.GetCondomPercent(HV_MSMIDU,CurYrIdx)) * value;
      value := max(value,0);
      value := min(value,1);
      HV.SetCondomPercent(HV_MSMIDU,t, Value);
    end;
  //end;{t loop}

  (* Number of partners - male and Female *)
  //for t := CurYrIdx + 1 to GetHVFinalIndex(p) do
  //begin
    for s := HV_Male to HV_Female do
    begin
     (* Low risk *)
      value := 1;
      for i := 1 to HV_RN_IDUDrugSub do
        if not (i in HV_RN_NotInCEAorIM) then
          value := value * (1 + _AdjCoverage[i, t] * HV.GetImpactMatrix(HV_Data, i, HV_RGReducInNumPartners_Low)/ 100);
      if s = HV_Male then
      begin
        value :=  HV.GetNumPart(HV_LRH,CurYrIdx) * value;
        value := max(value,0);
        HV.SetNumPart(HV_LRH,t,value);
      end
      else
      begin
        value :=  HV.GetNumPart(HV_LRH_F3,CurYrIdx) * value;
        value := max(value,0);
        HV.SetNumPart(HV_LRH_F3,t,value);
      end;

      (* Med risk *)
      value := 1;
      for i := 1 to HV_RN_IDUDrugSub do
        if not (i in HV_RN_NotInCEAorIM) then
          value := value * (1 + _AdjCoverage[i, t] * HV.GetImpactMatrix(HV_Data, i, HV_RGReducInNumPartners_Med)/ 100);
      if s = HV_Male then
      begin
        value :=  HV.GetNumPart(HV_MRH,CurYrIdx) * value;
        value := max(value,0);
        HV.SetNumPart(HV_MRH,t,value);
      end
      else
      begin
        value :=  HV.GetNumPart(HV_MRH_F3,CurYrIdx) * value;
        value := max(value,0);
        HV.SetNumPart(HV_MRH_F3,t,value);
      end;

      (* High risk *)
      value := 1;
      for i := 1 to HV_RN_IDUDrugSub do
        if not (i in HV_RN_NotInCEAorIM) then
          if (i <> HV_RN_MSMOutreach) and (i <> HV_RN_MSMLub) then
            value := value * (1 + _AdjCoverage[i, t] * HV.GetImpactMatrix(HV_Data, i, HV_RGReducInNumPartners_High) /100);
      if s = HV_Male then
      begin
        value :=  HV.GetNumPart(HV_HRH,CurYrIdx) * value;
        value := max(value,0);
        HV.SetNumPart(HV_HRH,t,value);
      end
      else
      begin
        value :=  HV.GetNumPart(HV_HRH_F3,CurYrIdx) * value;
        value := max(value,0);
        HV.SetNumPart(HV_HRH_F3,t,value);
      end;

      (* IDU risk *)
      value := 1;
      for i := 1 to HV_RN_IDUDrugSub do
        if not (i in HV_RN_NotInCEAorIM) then
          value := value * (1 + _AdjCoverage[i, t] * HV.GetImpactMatrix(HV_Data, i, HV_RGReducInNumPartners_IDU)/ 100);
      if s = HV_Male then
      begin
        value :=  HV.GetNumPart(HV_IDU,CurYrIdx) * value;
        value := max(value,0);
        HV.SetNumPart(HV_IDU,t,value);
      end
      else
      begin
        value :=  HV.GetNumPart(HV_IDU_F3,CurYrIdx) * value;
        value := max(value,0);
        HV.SetNumPart(HV_IDU_F3,t,value);
      end;
    end;

   (* MSM risk *)
    value := 1;
    for i := 1 to HV_RN_IDUDrugSub do
      if not (i in HV_RN_NotInCEAorIM) then
        value := value * (1 + _AdjCoverage[i, t] * HV.GetImpactMatrix(HV_Data, i, HV_RGReducInNumPartners_MSM) / 100);
    value :=  HV.GetNumPart(HV_MSM,CurYrIdx) * value;
    value := max(value,0);
    HV.SetNumPart(HV_MSM,t,value);

    if HV.GetNumMSMRiskGroups > 0 then
    begin
      (* Low risk *)
      value := 1;
      for i := 1 to HV_RN_IDUDrugSub do
        if not (i in HV_RN_NotInCEAorIM) then
          value := value * (1 + _AdjCoverage[i, t] * HV.GetImpactMatrix(HV_Data, i, HV_RGReducInNumPartners_MSM_Low)/ 100);
      value :=  HV.GetNumPart(HV_MSMLR,CurYrIdx) * value;
      value := max(value,0);
      HV.SetNumPart(HV_MSMLR,t,value);

      (* Med risk *)
      value := 1;
      for i := 1 to HV_RN_IDUDrugSub do
        if not (i in HV_RN_NotInCEAorIM) then
          value := value * (1 + _AdjCoverage[i, t] * HV.GetImpactMatrix(HV_Data, i, HV_RGReducInNumPartners_MSM_Med)/ 100);
      value :=  HV.GetNumPart(HV_MSMMR,CurYrIdx) * value;
      value := max(value,0);
      HV.SetNumPart(HV_MSMMR,t,value);

      (* High risk *)
      value := 1;
      for i := 1 to HV_RN_IDUDrugSub do
        if not (i in HV_RN_NotInCEAorIM) then
          value := value * (1 + _AdjCoverage[i, t] * HV.GetImpactMatrix(HV_Data, i, HV_RGReducInNumPartners_MSM_High) /100);
      value :=  HV.GetNumPart(HV_MSMHR,CurYrIdx) * value;
      value := max(value,0);
      HV.SetNumPart(HV_MSMHR,t,value);

      (* MSM_IDU *)
      value := 1;
      for i := 1 to HV_RN_IDUDrugSub do
        if not (i in HV_RN_NotInCEAorIM) then
        value := value * (1 + _AdjCoverage[i, t] * HV.GetImpactMatrix(HV_Data, i, HV_RGReducInNumPartners_MSM_IDU) /100);
      value :=  HV.GetNumPart(HV_MSMIDU,CurYrIdx) * value;
      value := max(value,0);
      HV.SetNumPart(HV_MSMIDU,t,value);
    end;
  //end;{t loop}

  (* Age at first Sex *)
  //for t := CurYrIdx + 1 to GetHVFinalIndex(p) do
  //begin
    for s := HV_Male to HV_Female do
    begin
      value := 0;
      //value := AdjCoverage[t, HV_IntCommMobilization]/100 * GetHVImpactMatrix(p, HV_Data, HV_IntCommMobilization, HV_RGIncreaseAgeFirstSex)/100;
      for i := 2 to HV_RN_IDUDrugSub do
      if not (i in HV_RN_NotInCEAorIM) then
      begin
          if (i = HV_RN_SecondaryTeachers) or (i = HV_RN_OutOfSchoolYouth) then
          begin
            value := value + _AdjCoverage[i, t] / (HV_DP_Get_Pop_Age(p, t-1, HV_DP_BothSexes,
                  DP_A15_19, DP_A15_19) / HV_DP_Get_Pop_Age(p, t-1, HV_DP_BothSexes,
                  DP_A15_19, DP_A45_49)) * HV.GetImpactMatrix(HV_Data, i, HV_RGIncreaseAgeFirstSex)
      end
      else
        value := value + _AdjCoverage[i, t] * HV.GetImpactMatrix(HV_Data, i, HV_RGIncreaseAgeFirstSex);
      end;
      value :=  HV.GetAgeFirstSex(s,CurYrIdx) + value;
      value := max(value,0);
      HV.SetAgeFirstSex(s,t,value)
    end;
 // end;{t loop}

  (* IDU *)
  //for t := CurYrIdx + 1 to GetHVFinalIndex(p) do
  //begin
    value := 1;
    for i := 1 to HV_RN_IDUDrugSub do
      if not (i in HV_RN_NotInCEAorIM) then
        value := value * (1 + _AdjCoverage[i, t] * HV.GetImpactMatrix(HV_Data, i, HV_RGReducNeedleSharing_IDU) / 100)
                       * (1 + _AdjCoverage[i, t] * HV.GetImpactMatrix(HV_Data, i, HV_RGReducUnsafeInjectBehavior) / 100);
    value := HV.GetPerIDUSharing(CurYrIdx) * Value;
    value := max(value,0);
    HV.SetPerIDUSharing(t,value);
  //end;{t loop}

  (* Transfer local must be done here to affect behavioral change.
     Otherwise, behavior change editor must be opened twice *)

  (* condom use *)
  for r := HV_None to HV_MSMHR do
  begin
    _CondomPercent[r,t] := HV.GetCondomPercent(r,t);
    (* Transfer condom-use values to InterpolatedVal condom use editor
       in order to visualize changes that may have resulted from the impact matrix *)
    HV.SetCondomUseInterpolatedVal(r, t, HV.GetCondomPercent(r,t));
  end;

  (* partners and acts *)
  for r := HV_AllRisk to HV_IDU_F3 do
  begin
    _NumPart[r,t] := HV.GetNumPart(r,t);
    _SexActs[r,t] := HV.GetSexActs(r,t);
  end;

  (* needle sharing *)
  _PerIDUSharing[t] := HV.GetPerIDUSharing(t)/100;

  (* age of first sex *)
  for s := HV_BothSexes to HV_Female do
    _AgeFirstSex[s,t] := HV.GetAgeFirstSex(s,t);

end;

procedure THVProj.Calc_ARTCoverageDisplay;
var
  s,r,t,v,h       : byte;
  TotalNeed,
  TotalReceiving  : HV_Sex;
  CD4LowerLimit   : array[HV_CD4_GT500..HV_CD4_LT50] of double;
  NeedART,
  ReceivingART    : array[HV_Male..HV_Female,HV_None..HV_MSMIDU] of HV_TDoubleDynYearArray;
  I: Integer;
  j: Integer;
begin
  for I := HV_Male to HV_Female do
    for j := HV_None to HV_MSMIDU do
    begin
      SetLength(NeedART[i, j], GetGBFinalYearIdx(p) + 1);
      SetLength(ReceivingART[i, j], GetGBFinalYearIdx(p) + 1);
    end;

  (* Initialize arrays *)
  CD4LowerLimit[HV_CD4_GT500] := 500;
  CD4LowerLimit[HV_CD4_350_500] := 350;
  CD4LowerLimit[HV_CD4_250_349] := 250;
  CD4LowerLimit[HV_CD4_200_249] := 200;
  CD4LowerLimit[HV_CD4_100_199] := 100;
  CD4LowerLimit[HV_CD4_50_99] := 50;
  CD4LowerLimit[HV_CD4_LT50] := 0;

  for s := HV_Male to HV_Female do
    for r := HV_None to HV_MSMIDU do
      for t := HV_Min_Year to GetGBFinalYearIdx(p) do
      begin
        NeedART[s,r,t] := 0.0;
        ReceivingART[s,r,t] := 0.0;
      end;

  for t := GetGBCalcYearIdx(p) to HV.GetFinalIndex do
  begin
    (* Determine number needing ART and number receiving ART *)
    for s := HV_Male to HV_Female do
      for r := HV_None to HV_MSMIDU do
        for v := HV_RN_UnV to HV_RN_NoProt do
          for h := HV_CD4_GT500 to HV_CD4_LT50 do
          begin
           // if CD4LowerLimit[h] < HV_GetDPCD4ThreshHoldAdults(p,t)  then
            //  NeedART[s,r,t]:= NeedART[s,r,t]+ GetHVAdults(p,s,r,h,v,t);
            //ReceivingART[s,r,t]:= ReceivingART[s,r,t] + GetHVAdults(p,s,r,h+10,v,t);
           NeedART[s,r,t]:= NeedART[s,r,t]+ _Adults[s,r,h,v,t] ;//_RiskGrpARTCov[s,r,h,v,t]  ;
           ReceivingART[s,r,t]:= ReceivingART[s,r,t] + _Adults[s,r,h+10,v,t] ;
          end;
    (* Determine ART coverage by risk group *)
    for s := HV_Male to HV_Female do
      for r := HV_None to HV_MSMIDU do
      begin

        if (NeedART[s,r,t]+ReceivingART[s,r,t]>0 )  then
          HV.SetARTCoverageByRG(s,r,t, (ReceivingART[s,r,t] / (NeedART[s,r,t]+ReceivingART[s,r,t]))*100)
        else
          HV.SetARTCoverageByRG(s,r,t,0.0);
      end;
    (* Calculate ART Coverage for both sexes by risk group *)
    for r := HV_None to HV_MSMIDU do
    begin
        if ((NeedART[HV_Male,r,t] + NeedART[HV_Female,r,t]) + (ReceivingART[HV_Male,r,t] + ReceivingART[HV_Female,r,t])) > 0 then
           HV.SetARTCoverageByRG(HV_BothSexes,r,t,((ReceivingART[HV_Male,r,t] + ReceivingART[HV_Female,r,t])/
                                              ((NeedART[HV_Male,r,t] + NeedART[HV_Female,r,t]) +
                                                (ReceivingART[HV_Male,r,t] + ReceivingART[HV_Female,r,t]))*100))
         
       else
          HV.SetARTCoverageByRG(HV_BothSexes,r,t,0.0);
    end;

    (* Determine total number needing and receiving ART *)
    for s := HV_BothSexes to HV_Female do
    begin
      TotalNeed[s] := 0.0;
      TotalReceiving[s] := 0.0;
      for r := HV_None to HV_MSMIDU do
      begin
        if s = HV_BothSexes then
        begin
          TotalNeed[s] := TotalNeed[s] + NeedART[HV_Male,r,t] + NeedART[HV_Female,r,t];
          TotalReceiving[s] := TotalReceiving[s] + ReceivingART[HV_Male,r,t] + ReceivingART[HV_Female,r,t];
        end
        else
        begin
          TotalNeed[s] := TotalNeed[s] + NeedART[s,r,t];
          TotalReceiving[s] := TotalReceiving[s] + ReceivingART[s,r,t];
        end;
      end;
    end;
    (* Calculate Total ART Coverage *)
    for s := HV_BothSexes to HV_Female do
      if (TotalNeed[s] + TotalReceiving[s] > 0) then
        HV.SetTotalARTCoverage(s,t,(TotalReceiving[s] / (TotalNeed[s] + TotalReceiving[s]))*100)
      else
        HV.SetTotalARTCoverage(s,t,0.0);
  end;// t loop
end;

procedure THVProj.Init_AdjCoverage;
var
  i, t : byte;
begin
  for i := 1 to HV_RN_MaxInterventions do
    for t := HV_MIN_YEAR to GetGBFinalYearIdx(p) do
      _AdjCoverage[i, t] := 0;
end;

end.
