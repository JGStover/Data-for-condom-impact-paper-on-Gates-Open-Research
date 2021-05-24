unit DPProj;

{ Contains the calculations to project the population and related arrays for
  a single year projection }

interface

uses
  Forms, sysUtils, dialogs, controls, math, GBData, DPConst, DPDefs, DPData,
  GBSessionUnit, HAProj, HVProj, HAConst;

type
TDPProj = Class (TObject)
const
  Max_CD_categories = 8;
  D14Yrs            = 14;




  Retention_Rate_CD4_count    : array[1..Max_CD_categories] of double  = (0.986, 0.799, 0.761, 0.757, 0.666, 0.644, 0.460, 0.50);
  Retention_Rate_CD4_percent  : array[1..Max_CD_categories] of double  = (0.906, 0.759, 0.787, 0.795, 0.785, 0.756, 0.813, 0.700);
  Age_0_CD4_count             : array[1..Max_CD_categories] of double  = (5, 3, 2, 6, 6, 17, 16, 45);
  Age_0_CD4_percent           : array[1..Max_CD_categories] of double  = (1, 3, 5, 9, 10, 12, 26, 34);
  CD_count_categories         : array[1..Max_CD_categories] of integer = (10000, 2000, 1500, 1000, 750, 500, 350, 200);
  CD_percent_categories       : array[1..Max_CD_categories] of integer = (100, 40, 30, 25, 20, 15, 10, 5);

  Fert_age_discount : array[1..6] of single = (2.21, 1.05, 0.79, 0.75, 0.64, 0.58);
  Fert_CD4_discount : array[DP_CD4_GT500..DP_CD4_LT50] of single = (1.0, 0.91, 0.77, 0.57, 0.38, 0.29, 0.28);


  (* Distribution of HIV+ children at age 10 from HPPMCS
     <200 = 28%; 200-350 = 21%; 350-500 = 23%; >500 = 28%
     note that if there is a change to the values in A15_CD4_distribution array, they will also have to
     be changed for the ep4 file in EPPDataRW, labeled as HIVPOS15_CD4 *)
  A15_CD4_distribution        : array[DP_CD4_GT500..DP_CD4_LT50] of single = (0.056, 0.112, 0.112, 0.07, 0.14, 0.23, 0.28);
  A15_ART_CD4_distribution    : array[DP_CD4_250_349..DP_CD4_LT50] of single = (0.11, 0.23, 0.23, 0.14, 0.29);
private
  p,
  a10                          : byte;
  PTR1,PTR3,
  BFTR,
  pr,
  births1524,
  UnmetNeed,
  UnmetNeedMale,
  UnmetNeedFemale,
  FirstYearNewInfFromBF,
  SecondYearNewInfFromBF,
  TripleArtPrePreg,
  TripleArtDurPreg,
  sum1, sum2, sum3,
  births, HIVbirths, HIVbirthsToExistingInfections,
  CtxEff, CtxArtEff,
  CTXcoverage,
  TotalOnCTX,
  TotalNeedForCTX,
  NeedForCTX,
  temp, temp_all, AdjFactor,
  adults, adultMales,
  AdultFemales,
  DeathsAverted,
  PercentInProgram,
  NewInfFromBFLT6,
  NewInfFromBF6TO12,
  CumNewInfFromBF,
  PropNewInfFromBF,
  IncidenceAdjustmentFactor    : single;
  CD4_Threshold,
  PropLT200, Prop200to350,
  PropLT350,
  PropGE350,
  EligibleChildren,
  ExcessRatio,
  OptA_Transmission_Rate,
  OptB_Transmission_Rate,
//  ART_transmission_rate,
//  PropIncidentInfections,
  No_ART, On_ART,
  No_ART_14, On_ART_14,
  No_ART_50, On_ART_50,
  v1, v2, v3, v4, v5,
  v6, v7, IncidentInfectionsBF,
  alpha_mult                   : double;
  LastAge                      : byte;
  pop1                         : DP_TPopType;
  p5                           : FiveAgeArray;
  sa                           : SingleAgeArray;
  HAObj                        : THAProj;
  HVObj                        : THVProj;
  AidsDeathsByAge              : TDeathsBySingleAge;
  TotalOldAIDS,
  PTR2                         : array of single;
  durMatrix                    : array [DP_Male..DP_Female,0..(DP_MaxSingleAges - 14), DP_CD4_GT500 ..DP_CD4_LT50] of integer;      // duration of infection based on CD4 at death for orphan calcs
  durARTMatrix                 : array [DP_Male..DP_Female, DP_CD4_GT500 .. DP_CD4_LT50] of array of Integer; // duration of ART at death for orphan calcs
  mu                           : array [1 .. 4, DP_CD4_GT500 .. DP_CD4_LT50, DP_Male .. DP_Female] of double;   // duration of infection based on CD4 at death for orphan calcs
  alpha                        : array [1 .. 4, DP_CD4_GT500 .. DP_CD4_LT50, DP_Male .. DP_Female, DP_D_ARTlt6M .. DP_D_ARTgt12M] of double;
  PrevRatio                    : array [0 .. 80, DP_Male .. DP_Female] of single;
  Elig                         : array[0..14, DP_BothSexes ..DP_Female, DP_H_All .. DP_CD4_Per_LT5, DP_D_All .. DP_P_BF12] of double;
  NewHIVbyAge                  : array [0 .. 80] of single;
  NewHIV,
  NeededHIV, adultHIV          : array [DP_BothSexes .. DP_Female] of single;
  TreatPercent                 : array [DP_NoProphylaxis .. DP_ART_BF] of single;
  ARTThisYearBySex,
  ARTLastYearBySex             : array[DP_Male..DP_Female] of single;

  (* For CalcState *)
  Count_dist,
  Percent_dist                 : array[DP_A0..DP_A14, 1..DP_CD4_MaxCategories] of array of double;
  ARTTwoYearsAgoBySexCalcState : array[DP_Male..DP_Female] of double;
  GBSession                    : TGBSession;
  GB                           : TGBVariables;
  DP                           : TDPData;
public
  constructor Create(GBSession : TGBSession; Proj : byte);
  function ProjectPopulationBySingleYears : boolean;
  procedure CalcOrphans;
  function MaternalOrphanCalcs : Boolean;
  function PaternalOrphanCalcs : Boolean;
  function DualOrphanCalcs : Boolean;
  function TotalOrphanCalcs : Boolean;
  procedure CreateInfDur( t, s: Byte);
  function getAgeAtInfection (s, A, CD4 : byte): Byte;
  procedure setDurARTIni(s: Byte);
  procedure CalcBFtransmission(p, t, m1, m2: integer);
  procedure Do_AIDS_calculations(p, t: integer);
  procedure CalcNewAdultInfections(t : integer);
  procedure CalcHIVBirths(t : integer);
  procedure CalcNeedForCTX(t : integer);
  procedure CalcNumChildrenEligTreat(t : integer);
  function DP_Do_Calculations : boolean;
  procedure SumPopulationByHIVStatus(const t : integer);
  procedure Calc_BothSexes;
  procedure Calc_ExtraIndicators;
  destructor Destroy; override;
  (* CalcState*)
  procedure CaptureCalcStateData(pop1 : DP_TPopType);
  procedure RetrieveCalcStateData(pop1 : DP_TPopType);
end;

implementation

uses
  Windows, System.Diagnostics, GBConst, RetStrD, GBDefs, GBUtil, GBMathRtns, GBmdUtil,
  DPUtil, DPLink, DPProjUt, RNData, DpGbData, HAData, HAUtil, DPResultProj, DPDataUtil;

Constructor TDPProj.Create(GBSession : TGBSession; Proj : byte);
var
  a, c, s : byte;
begin
  p := Proj;
  self.GBSession := GBSession;
  GB := GBSession.GB;
  DP := GBSession.GetProjection(proj).DP;

  SetLength(TotalOldAIDS, GB.GetNumDynYears(p));
  SetLength(PTR2, GB.GetNumDynYears(p));


  for a := DP_A0 to DP_A14 do
    for c := 1 to DP_CD4_MaxCategories do
    begin
      SetLength(Count_dist[a,c], GB.GetNumDynYears(p));
      SetLength(Percent_dist[a,c], GB.GetNumDynYears(p));
    end;

  for s := DP_Male to DP_Female do
    for c := DP_CD4_GT500 to DP_CD4_LT50 do
      SetLength(durARTMatrix[s,c], GB.GetNumDynYears(p));
end;

(* Project the single year population values *)
function TDPProj.ProjectPopulationBySingleYears : boolean;
var
  PopAdj,
  sr, sr1, mr, dths,
  TempBirths,
  TempF,TempM,
  temp,
  RegTot,
  ratio,
  SumPLHIVOnART,
  SumSpec,
  SumEPP,
  AdjustFactor,
  TempDur,TempMDur,TempFDur : double;
  DeathsByAge               : TDeathsBySingleAge;
  CalcStateIdx,
  Idx ,r                    : integer;
  a,s,t,t1,c,d,d1,a5        : byte;
  migr                      : DP_TPopType2;
  CompOk                    : Boolean;
  Females                   : array [15..49] of double;
//  DataRetrieved             : boolean;
  value                     : double;
//  stopwatch                 : TStopWatch;
begin
//  stopwatch := TStopWatch.Create;
//  stopwatch.start;

  (* Nil the HVProj so that if it is not created properly, we can use Assigned
     to determine whether or not we should free the object in the destructor. *)
  HVObj := nil;
  HAObj := nil;

//  DataRetrieved := false;
  If GB.ModuleActive(p, GB_HV) then
    HVObj := THVProj.Create(p);

  if GB.ModuleActive(p, GB_HA) then
  begin
    HAObj := THAProj.Create(p);
    HA_PartnerChoiceCalc(p,1);
  end;

  result := true;

  for s := DP_Male to DP_Female do
    for a  := 0 to DP_MaxSingleAges do
    begin
      SetLength(AidsDeathsByAge[a,s], GB.GetNumDynYears(p));
      SetLength(DeathsByAge[a,s], GB.GetNumDynYears(p));
    end;

  for t := 1 to DP.GetFinalIndex do
    for s := DP_Male to DP_Female do
      for a  := 0 to DP_MaxSingleAges do
        DP.SetDeathsBySingleAge(t, a, s, 0);

  Try
    setlength(pop1,DP.GetFinalIndex+1,DP_Female+1,DP_MaxSingleAges+1,DP_H_MAX+1,DP_D_Max+1);
    setlength(migr,DP_MaxSingleAges+1,DP.GetFinalIndex+1,DP_Female+1,DP_H_MAX+1);

    { Initialize arrays to zero }
    for t := 1 to DP.GetFinalIndex do
    begin
      DP.SetAdultDeaths(t, 0);
      for s := DP_BothSexes to DP_Female do
      begin
        for a := DP_AllAges to DP_MAX_AGE do
        begin
          DP.SetBirths(t, a, s, 0);
          DP.SetDeaths(t, a, s, 0);
          DP.SetAIDSDeaths(t, a, s, 0);
        end;
        DP.SetChildARTCalc(DP_NumOnArt, s, t, 0);
      end;
    end;

    for t := 1 to DP.GetFinalIndex do
      for s := DP_Male to DP_Female do
        for a := DP_A0 to DP_A80 do
          DP.SetNewInfectionsBySingleAge(t, s, a, 0);

    (* initalize values to zero *)
    for t := 1 to DP.GetFinalIndex do
      for s := DP_BothSexes to DP_Female do
        for a := DP_A0_4 to DP_A80_up do
        begin
			(*EPW start*)
          DP.SetNewInfections(t,a,s,0);
          for d := DP_H_All to DP_H_Max do
            for t1 := 1 to 10 do
              DP.SetNewlyOnART(t,t1,s,a,d,0);

          if GB.ModuleActive(p, GB_TB) then
          for c := DP_H_All to DP_H_Max do
            for t1 := 1 to 10 do
              for d := DP_HIVNeg to DP_HIVARTgt12m do
              DP_TI_SetHIVbyAge(p,t,t1,s,a,c,d,0);

          for d := DP_HIVNeg to DP_HIVARTgt12m do
            DP.SetHIVAgingRate(t, s, a, d,0);
			(*EPW end*)
        end;

    (* initalize values to zero *)
    for t := 1 to DP.GetFinalIndex do
    begin
      DP.SetHIVPos15YearOlds(t, 0);
      DP.SetHIVPos50YearOlds(t, 0);
    end;

    for t := 0 to DP.GetFinalIndex do
      for s := DP_BothSexes to DP_Female do
        for a := 0 to DP_MaxSingleAges do
          for c := DP_H_All to DP_H_Max do
            for d := DP_D_All to DP_D_Max do
              pop1[t,s,a,c,d] := 0;

   { Initialize 5 year age group array to zero }
    for t := 1 to DP.GetFinalIndex do
    begin
      TotalOldAIDS[t] := 0;
      for s := DP_BothSexes to DP_Female do
        for a := DP_AllAges to DP_A80_Up do
        begin
          DP.SetHIV(t, s, a, 0);
        end;
    end;

    { Split base year population into single ages }
    if not GB.GetInCalcStateMode(p) then
    begin
      if DP.GetUseFYrSingleAge = false then
      begin
        for s := DP_Male to DP_Female do
        begin
          for a := 0 to DP_MaxSingleAges do
            sa[a] := 0;
          for a := DP_AllAges to DP_MAX_AGE do
            p5[a] := DP.GetPop(DP_MIN_YEAR, s, a);

          Calc_Single_Ages(GB, DP,p5, sa, s);

          for a := 0 to DP_MaxSingleAges do
          begin
            pop1[DP_MIN_YEAR, s, a, DP_H_No, DP_D_HIVNeg] := sa[a];
            pop1[DP_MIN_YEAR, s, a, DP_H_All,DP_D_All] := sa[a];
          end;
        end;
      end
      else
      begin
      (* Put Single Age Inputs in Population Array for First Year *)
        for s := DP_Male to DP_Female do
          for a := 0 to DP_MaxSingleAges do
          begin
            pop1[DP_MIN_YEAR, s, a, DP_H_No, DP_D_HIVNeg] := DP.GetBigPop(DP_MIN_YEAR, s, a);
            pop1[DP_MIN_YEAR, s, a, DP_H_All,DP_D_All] := DP.GetBigPop(DP_MIN_YEAR, s, a);
          end;
      end;
    end;
    (* If we are in CalcStateMode and we are not trying to generate the calc state,
       then retrieve the calcstate values that were saved *)
    if GB.GetInCalcStateMode(p) then
    begin
      RetrieveCalcStateData(pop1);
//      DataRetrieved := true;
    end;



    (* Calculate the number of immigrants for all ages by using the Net immigrants total for each sex from *)
    for t := GB.GetCalcYearIdx(p) to DP.GetFinalIndex do
      for a := 1 to DP_MAX_AGE do
        for s := DP_Male to DP_Female do
          DP.SetMigrRate(t, s, a, DP.GetMigrAgeDist(t, s, a) / 100 * DP.GetMigrRate(t, s, DP_AllAges));

     for t := GB.GetCalcYearIdx(p) to DP.GetFinalIndex do
       for a := DP_AllAges to DP_MAX_AGE do
         DP.SetMigrRate(t, DP_BothSexes, a, DP.GetMigrRate(t, DP_Male, a) + DP.GetMigrRate(t, DP_Female, a));

    (* Split migrants into single ages *)
    for t := GB.GetCalcYearIdx(p) to DP.GetFinalIndex do
    begin
      for s := DP_Male to DP_Female do
      begin
        for a := 0 to DP_MaxSingleAges do
          sa[a] := 0;

        for a := DP_AllAges to DP_Max_Age do
          p5[a] := DP.GetMigrRate(t,s,a);

        Calc_Single_Ages(GB,DP,p5,sa,s);

        for a := 0 to DP_MaxSingleAges do
          migr[a,t,s,DP_H_All] := sa[a];
      end;
    end;

   (* if creating CalcState so we are not in CalcState mode *)
    if not GB.GetInCalcStateMode(p) then
    begin
      t := 1;

      { Calculate births in base year }
      for a := DP_AllAges to DP_A45_49 do
        DP.SetBirths(t, a, DP_BothSexes, 0);
      for a := 15 to 49 do
      begin
        a5 := trunc(a/5) + 1;
        if DP.GetASFR (t,DP_AllAges) <> 0 then
        begin
          sr := Get_Survival_Rate(DP,t, DP_Female,a-1);
          Tempbirths := (pop1[t,DP_Female,a,DP_H_All,DP_D_All] + pop1[t,DP_Female,a-1,DP_H_All,DP_D_All] * sr)/2
                        * DP.GetTFR(t) * DP.GetASFR(t, a5) / DP.GetASFR(t, DP_AllAges) / 5
        end
        else
        begin
          TempBirths := 0;
          DP_AddLogEntry(GB, DP_DPAIM,'DPPROJ1','~3-10','DP.GetASFR (t,DP_AllAges)',
                             DP_EZ,DP.GetASFR (t,DP_AllAges), 't = ' + IntToStr(t));
        end;

        DP.SetBirths(t, DP_AllAges, DP_BothSexes,
                          DP.GetBirths(t, DP_AllAges,DP_BothSexes) + Tempbirths);
        DP.SetBirths(t, DP_AllAges, DP_Male,
                          DP.GetBirths(t, DP_AllAges,DP_Male)
                          + Tempbirths * GetDP_SexBirthRatioPercent(DP,t,DP_Male));
        DP.SetBirths(t, DP_AllAges, DP_Female,
                          DP.GetBirths(t, DP_AllAges,DP_Female)
                          + Tempbirths * GetDP_SexBirthRatioPercent(DP,t,DP_Female));
        DP.SetBirths(t, a5, DP_BothSexes,
                          DP.GetBirths(t, a5,DP_BothSexes) + TempBirths);
      end; {a loop }

      DP.SetPerinatalTransmission(t, DP.GetTransEffAssump(DP_Data, DP_NoProphylaxis, DP_Perinatal));

    { Calculate deaths in base year }
      for s := DP_Female downto DP_Male do
      begin
        for a := 1 to DP_MaxSingleAges - 1 do
        begin
          a5 := trunc(a/5) + 1;
          sr := Get_Survival_Rate(DP,t,s,a);
          mr := (migr[a-1,t,s,DP_H_All] + migr[a,t,s,DP_H_All])/2;
          dths := (pop1[t,s,a-1,DP_H_All,DP_D_All]) * (1 - sr)
                + mr*(1 - sr)/2;
          DeathsByAge[a,s,t] := dths;
          if a5 > DP_Max_Age then a5 := DP_Max_Age;
          DP.SetDeaths(t, a5, s, DP.GetDeaths(t, a5, s) + dths);
          DP.SetDeaths(t, DP_ALLAGES, s,
                 DP.GetDeaths(t, DP_ALLAGES, s) + dths);
          if (a >= 15) and (a <= 49) then
            DP.SetAdultDeaths(t, DP.GetAdultDeaths(t) + dths);
        end; { a loop }
        { Last age group }
        sr := Get_Survival_Rate(DP, t, s, 80);
        mr := (migr[79,t,s,DP_H_All] + migr[80,t,s,DP_H_All])/2;
        dths := (pop1[t,s,79,DP_H_All,DP_D_All] + pop1[t,s,80,DP_H_All,DP_D_All]) * (1 - sr)
              + mr * (1 - sr)/2;
        DeathsByAge[80,s,t] := dths;
        DP.SetDeaths(t, DP_MAX_AGE, s, DP.GetDeaths(t, DP_MAX_AGE, s) + dths);
        DP.SetDeaths(t, DP_ALLAGES, s,
                 DP.GetDeaths(t, DP_ALLAGES, s) + dths);
        { First age group }
        sr := Get_Survival_Rate(DP, t, s, 0);
       // mr := (migr[0,t,s,DP_H_All,DP_All])/2;
        dths := (DP.GetBirths(t, DP_AllAges,s)) * (1 - sr)
              + 0.5 * (1 - sr)/2;
        DeathsByAge[0,s,t] := dths;
        DP.SetDeaths(t, DP_A0_4, s, DP.GetDeaths(t, DP_A0_4, s) + dths);
        DP.SetDeaths(t, DP_ALLAGES, s,
                 DP.GetDeaths(t, DP_ALLAGES, s) + dths);
      end; { s loop }

      DP.SetDeaths(t, DP_ALLAGES, DP_BothSexes,
                 DP.GetDeaths(t, DP_ALLAGES, DP_Male)
                 + DP.GetDeaths(t, DP_ALLAGES, DP_Female));

      DP.SetAdjLE(t,DP_Male, DP.GetLE(t,DP_Male));
      DP.SetAdjLE(t,DP_Female, DP.GetLE(t,DP_Female));
      DP.SetAdjLE(t,DP_BothSexes, (DP.GetLE(t,DP_Male)
                                    * DP.GetPop(t,DP_Male,DP_AllAges)
                                    + DP.GetLE(t,DP_Female)
                                    * DP.GetPop(t,DP_Female,DP_AllAges))
                                    / DP.GetPop(t,DP_BothSexes,DP_AllAges));
	  end;
   { Project population and deaths for total region }
    for t := 2 to DP.GetFinalIndex do
    begin
      No_ART := 0;
      On_ART := 0;
      No_ART_14 := 0;
      On_ART_14 := 0;
      No_ART_50 := 0;
      On_ART_50 := 0;
      for s := DP_Female downto DP_Male do
      begin
        for a := 1 to DP_MaxSingleAges - 1 do
        begin
          a5 := trunc(a/5) + 1;
          sr := Get_Survival_Rate(DP,t,s,a);
          sr1 := sr;
          { Include excess IDU mortality }
          if (a>14) and (a<50) and (pop1[t-1,DP_Male,a-1,DP_H_All,DP_D_All] > 0) then
          begin
            if DP.GetEPPSexRatio(1,t) <> -1 then
            begin
              { Calculate proportion of HIV+ IDU for current sex }
              if s = DP_Male then ratio :=  1/(1+DP.GetEPPSexRatio(1,t))
                             else ratio :=  DP.GetEPPSexRatio(1,t)/(1+DP.GetEPPSexRatio(1,t));
              temp := sr;
              temp := 1 - temp;
              if (pop1[t-1,s,a-1,DP_H_All,DP_D_All] = 0) or (pop1[t-1,DP_Male,a-1,DP_H_All,DP_D_All]
                        + pop1[t-1,DP_Female,a-1,DP_H_All,DP_D_All] = 0) then
                temp := sr
              else
                temp := temp + max(0,(DP.GetIDUMortality(1)/100 - temp))  {max(0,(0.025 - temp))}
                      * DP.GetPropIDU_WB(0,t)
                      * ratio
                     /(pop1[t-1,s,a-1,DP_H_All,DP_D_All]
                      /(pop1[t-1,DP_Male,a-1,DP_H_All,DP_D_All]
                        + pop1[t-1,DP_Female,a-1,DP_H_All,DP_D_All]));
              sr1 := 1 - temp;
            end;
          end;
          mr := (migr[a-1,t,s,DP_H_All] + migr[a,t,s,DP_H_All])/2;

          if (a>14) and (a<50)
          then
            dths := (pop1[t-1,s,a-1,DP_H_No,DP_D_All] * (1 - sr)
                     + (pop1[t-1,s,a-1,DP_H_All,DP_D_All]
                        - pop1[t-1,s,a-1,DP_H_No,DP_D_All]) * (1 - sr1))
                    + mr * (1 - sr)/2
          else
            dths := pop1[t-1,s,a-1,DP_H_All,DP_D_All] * (1 - sr)
                    + mr * (1 - sr)/2;
          { Include COVID-19 deaths }
          if ((DP.GetFirstYear + t - 1) <= DP_COVID19FinalYr) and
             ((DP.GetFirstYear + t - 1) >= DP_COVID19FirstYr) then
          begin
           { Multiply total COVID deaths by percent by age
             Percent by age is for 5-year age group so divide by 5 }
            temp := DP.GetCOVID19DeathRate(t, s, 0) *
                       DP.GetCOVID19DeathAgeDist(DP_Data, t, s, a5) / 100 / 5;
            dths := dths + temp;
            if (pop1[t-1,s,a-1,DP_H_All,DP_D_All] > 0) then
              sr := max(0.0001, sr - temp / pop1[t-1,s,a-1,DP_H_All,DP_D_All])
          end;

          DeathsByAge[a,s,t] := dths;
          if a5 > DP_Max_Age then a5 := DP_Max_Age;
          DP.SetDeaths(t, a5, s, DP.GetDeaths(t, a5, s) + dths);
          if (a >= 15) and (a <= 49) then
            DP.SetAdultDeaths(t, DP.GetAdultDeaths(t) + dths);
          if pop1[t-1,s,a-1,DP_H_All,DP_D_All] = 0 then
          begin
            mr := 0;
            DP_AddLogEntry(GB, DP_DPAIM,'DPPROJ1','~3-20','pop1[t-1,s,a-1,DP_H_All,DP_D_All]',
                           DP_EZ,pop1[t-1,s,a-1,DP_H_All,DP_D_All],
                           'a = ' + IntToStr(a) + ',t = ' + IntToStr(t) + ',s = ' + IntToStr(s));
          end
          else
		        mr := ((migr[a-1,t,s,DP_H_All] + migr[a,t,s,DP_H_All]) / 2)
                       / pop1[t-1,s,a-1,DP_H_All,DP_D_All];

          { Progress those not infected }
          pop1[t,s,a,DP_H_No,DP_D_HIVNeg] := pop1[t-1,s,a-1,DP_H_No,DP_D_HIVNeg] * sr
                                 + pop1[t-1,s,a-1,DP_H_No,DP_D_HIVNeg]* mr * (1 + sr)/2;
          if pop1[t,s,a,DP_H_No,DP_D_HIVNeg] < 0 then pop1[t,s,a,DP_H_No,DP_D_HIVNeg] := 0;

          { Progress those who are infected }
          case a of
            0..4 : begin
                      for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
                      begin
                        { Not on ART }
                        for d := DP_P_Perinatal to DP_P_BF12 do
                        begin
                          pop1[t,s,a,c,d] := pop1[t-1,s,a-1,c,d] * sr
                                           + pop1[t-1,s,a-1,c,d] * mr * (1 + sr)/2;
                          if pop1[t,s,a,c,d] < 0 then pop1[t,s,a,c,d] := 0;
                        end;
                        { On ART }
                        for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
                        begin
                          pop1[t,s,a,c,d] := pop1[t-1,s,a-1,c,d] * sr
                                           + pop1[t-1,s,a-1,c,d] * mr * (1 + sr)/2;
                          if pop1[t,s,a,c,d] < 0 then pop1[t,s,a,c,d] := 0;
                        end;
                      end; {c loop}
                    end;
            5     : begin
                      for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
                      begin
                        { Not on ART }
                        for d := DP_P_Perinatal to DP_P_BF12 do
                        begin
                          pop1[t,s,a,c,d] := pop1[t-1,s,a-1,c,d] * sr
                                           + pop1[t-1,s,a-1,c,d] * mr * (1 + sr)/2;
                          if pop1[t,s,a,c,d] < 0 then pop1[t,s,a,c,d] := 0;
                        end;
                        { On ART }
                        for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
                        begin
                          pop1[t,s,a,c,d] := pop1[t-1,s,a-1,c,d] * sr
                                           + pop1[t-1,s,a-1,c,d] * mr * (1 + sr)/2;
                          if pop1[t,s,a,c,d] < 0 then pop1[t,s,a,c,d] := 0;
                        end;
                      end; {c loop}

                      { At age 5 convert from CD4 percent to CD4 count }
                      for d := DP_P_Perinatal to DP_D_ARTgt12m do
                      begin
                        v1 := pop1[t,s,a,DP_CD4_Per_LT5,  d];
                        v2 := pop1[t,s,a,DP_CD4_Per_5_10, d];
                        v3 := pop1[t,s,a,DP_CD4_Per_11_15,d];
                        v4 := pop1[t,s,a,DP_CD4_Per_16_20,d];
                        v5 := pop1[t,s,a,DP_CD4_Per_21_25,d];
                        v6 := pop1[t,s,a,DP_CD4_Per_26_30,d];
                        v7 := pop1[t,s,a,DP_CD4_Per_GT30, d];
                        for c := DP_CD4_Ped_GT1000 to DP_CD4_Ped_LT200 do
                          pop1[t,s,a,c,d] :=
                             v1 * DP.GetCD4Transition(DP_CD4_Per_LT5,   c) +
                             v2 * DP.GetCD4Transition(DP_CD4_Per_5_10,  c) +
                             v3 * DP.GetCD4Transition(DP_CD4_Per_11_15, c) +
                             v4 * DP.GetCD4Transition(DP_CD4_Per_16_20, c) +
                             v5 * DP.GetCD4Transition(DP_CD4_Per_21_25, c) +
                             v6 * DP.GetCD4Transition(DP_CD4_Per_26_30, c) +
                             v7 * DP.GetCD4Transition(DP_CD4_Per_GT30,  c);
                      end;
//                      { If age 5 convert from CD4 percent to CD4 count }
//                      if a = 5 then
//                      begin
//                        for d := DP_P_Perinatal to DP_P_BF12 do
//                        begin
//                          v1 := pop1[t,s,a,DP_CD4_Per_LT5,d];
//                          v2 := pop1[t,s,a,DP_CD4_Per_5_10,d];
//                          v3 := pop1[t,s,a,DP_CD4_Per_11_15,d];
//                          v4 := pop1[t,s,a,DP_CD4_Per_16_20,d];
//                          v5 := pop1[t,s,a,DP_CD4_Per_21_25,d];
//                          v6 := pop1[t,s,a,DP_CD4_Per_26_30,d];
//                          v7 := pop1[t,s,a,DP_CD4_Per_GT30,d];     d1 := DP_CD4_Ped_LT200;
//                          pop1[t,s,a,d1,d] :=
//                               v1 * DP.GetCD4Transition(DP_CD4_Per_LT5, d1) +
//                               v2 * DP.GetCD4Transition(DP_CD4_Per_5_10, d1) +
//                               v3 * DP.GetCD4Transition(DP_CD4_Per_11_15, d1) +
//                               v4 * DP.GetCD4Transition(DP_CD4_Per_16_20, d1) +
//                               v5 * DP.GetCD4Transition(DP_CD4_Per_21_25, d1) +
//                               v6 * DP.GetCD4Transition(DP_CD4_Per_26_30, d1) +
//                               v7 * DP.GetCD4Transition(DP_CD4_Per_GT30, d1);
//                          d1 := DP_CD4_Ped_200_349;
//                          pop1[t,s,a,d1,d] :=
//                               v1 * DP.GetCD4Transition(DP_CD4_Per_LT5, d1) +
//                               v2 * DP.GetCD4Transition(DP_CD4_Per_5_10, d1) +
//                               v3 * DP.GetCD4Transition(DP_CD4_Per_11_15, d1) +
//                               v4 * DP.GetCD4Transition(DP_CD4_Per_16_20, d1) +
//                               v5 * DP.GetCD4Transition(DP_CD4_Per_21_25, d1) +
//                               v6 * DP.GetCD4Transition(DP_CD4_Per_26_30, d1) +
//                               v7 * DP.GetCD4Transition(DP_CD4_Per_GT30, d1);
//                          d1 := DP_CD4_Ped_350_499;
//                          pop1[t,s,a,d1,d] :=
//                               v1 * DP.GetCD4Transition(DP_CD4_Per_LT5, d1) +
//                               v2 * DP.GetCD4Transition(DP_CD4_Per_5_10, d1) +
//                               v3 * DP.GetCD4Transition(DP_CD4_Per_11_15, d1) +
//                               v4 * DP.GetCD4Transition(DP_CD4_Per_16_20, d1) +
//                               v5 * DP.GetCD4Transition(DP_CD4_Per_21_25, d1) +
//                               v6 * DP.GetCD4Transition(DP_CD4_Per_26_30, d1) +
//                               v7 * DP.GetCD4Transition(DP_CD4_Per_GT30, d1);
//                          d1 := DP_CD4_Ped_500_749;
//                          pop1[t,s,a,d1,d] :=
//                               v1 * DP.GetCD4Transition(DP_CD4_Per_LT5, d1) +
//                               v2 * DP.GetCD4Transition(DP_CD4_Per_5_10, d1) +
//                               v3 * DP.GetCD4Transition(DP_CD4_Per_11_15, d1) +
//                               v4 * DP.GetCD4Transition(DP_CD4_Per_16_20, d1) +
//                               v5 * DP.GetCD4Transition(DP_CD4_Per_21_25, d1) +
//                               v6 * DP.GetCD4Transition(DP_CD4_Per_26_30, d1) +
//                               v7 * DP.GetCD4Transition(DP_CD4_Per_GT30, d1);
//                          d1 := DP_CD4_Ped_750_999;
//                          pop1[t,s,a,d1,d] :=
//                               v1 * DP.GetCD4Transition(DP_CD4_Per_LT5, d1) +
//                               v2 * DP.GetCD4Transition(DP_CD4_Per_5_10, d1) +
//                               v3 * DP.GetCD4Transition(DP_CD4_Per_11_15, d1) +
//                               v4 * DP.GetCD4Transition(DP_CD4_Per_16_20, d1) +
//                               v5 * DP.GetCD4Transition(DP_CD4_Per_21_25, d1) +
//                               v6 * DP.GetCD4Transition(DP_CD4_Per_26_30, d1) +
//                               v7 * DP.GetCD4Transition(DP_CD4_Per_GT30, d1);
//                          d1 := DP_CD4_Ped_GT1000;
//                          pop1[t,s,a,d1,d] :=
//                               v1 * DP.GetCD4Transition(DP_CD4_Per_LT5, d1) +
//                               v2 * DP.GetCD4Transition(DP_CD4_Per_5_10, d1) +
//                               v3 * DP.GetCD4Transition(DP_CD4_Per_11_15, d1) +
//                               v4 * DP.GetCD4Transition(DP_CD4_Per_16_20, d1) +
//                               v5 * DP.GetCD4Transition(DP_CD4_Per_21_25, d1) +
//                               v6 * DP.GetCD4Transition(DP_CD4_Per_26_30, d1) +
//                               v7 * DP.GetCD4Transition(DP_CD4_Per_GT30, d1);
//                        end;
//                        { Convert those on ART }
//                        for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
//                        begin
//                          v1 := pop1[t,s,a,DP_CD4_Per_LT5,d];
//                          v2 := pop1[t,s,a,DP_CD4_Per_5_10,d];
//                          v3 := pop1[t,s,a,DP_CD4_Per_11_15,d];
//                          v4 := pop1[t,s,a,DP_CD4_Per_16_20,d];
//                          v5 := pop1[t,s,a,DP_CD4_Per_21_25,d];
//                          v6 := pop1[t,s,a,DP_CD4_Per_26_30,d];
//                          v7 := pop1[t,s,a,DP_CD4_Per_GT30,d];
//                          pop1[t,s,a,DP_CD4_Ped_LT200,d] := v1 + v2 + 1/9 * v3;
//                          pop1[t,s,a,DP_CD4_Ped_200_349,d] := 8/9 * v3 + 3/13 * v4;
//                          pop1[t,s,a,DP_CD4_Ped_350_499,d] := 10/13 * v4 + 3/18 * v5;
//                          pop1[t,s,a,DP_CD4_Ped_500_749,d] := 15/18 * v5 + 6/15 * v6;
//                          pop1[t,s,a,DP_CD4_Ped_750_999,d] := 9/15 * v6 + 8/28 * v7;
//                          pop1[t,s,a,DP_CD4_Ped_GT1000,d] := 20/28 * v7;
//                          pop1[t,s,a,DP_CD4_Ped_Top,d] := 0;
//                        end;
//                      end; {a=5}
                    end;
            6..14 : begin
                      for c := DP_CD4_Ped_GT1000 to DP_CD4_Ped_LT200 do
                      begin
                        { Not on ART }
                        for d := DP_P_Perinatal to DP_P_BF12 do
                        begin
                          pop1[t,s,a,c,d] := pop1[t-1,s,a-1,c,d] * sr
                                           + pop1[t-1,s,a-1,c,d] * mr * (1 + sr)/2;
                          if pop1[t,s,a,c,d] < 0 then pop1[t,s,a,c,d] := 0;
                        end;
                        { On ART }
                        for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
                        begin
                          pop1[t,s,a,c,d] := pop1[t-1,s,a-1,c,d] * sr
                                           + pop1[t-1,s,a-1,c,d] * mr * (1 + sr)/2;
                          if pop1[t,s,a,c,d] < 0 then pop1[t,s,a,c,d] := 0;
                        end;
                      end; {c loop}
                    end;
            15    : begin
                      pop1[t,s,a,DP_CD4_GT500,DP_D_HIVPos] := 0;
                      pop1[t,s,a,DP_CD4_GT500,DP_D_ARTgt12M] := 0;
                      { Sum four highest child groups into highest adult group }
                      for c := DP_CD4_Ped_GT1000 to DP_CD4_Ped_500_749 do
                      begin
                        pop1[t,s,a,DP_CD4_GT500,DP_D_HIVPos] := pop1[t,s,a,DP_CD4_GT500,DP_D_HIVPos]
                                           + pop1[t-1,s,a-1,c,DP_P_Perinatal]
                                           + pop1[t-1,s,a-1,c,DP_P_BF0]
                                           + pop1[t-1,s,a-1,c,DP_P_BF7]
                                           + pop1[t-1,s,a-1,c,DP_P_BF12];
                        pop1[t,s,a,DP_CD4_GT500,DP_D_ARTgt12M] := pop1[t,s,a,DP_CD4_GT500,DP_D_ARTgt12M]
                                                  + pop1[t-1,s,a-1,c,DP_D_ARTlt6M]
                                                  + pop1[t-1,s,a-1,c,DP_D_ART6to12M]
                                                  + pop1[t-1,s,a-1,c,DP_D_ARTgt12M];
                      end;
                      { 350-500 group is the same for both adults and children }
                      c := DP_CD4_Ped_350_499;
                      pop1[t,s,a,DP_CD4_350_500,DP_D_HIVPos] := pop1[t-1,s,a-1,c,DP_P_Perinatal]
                                         + pop1[t-1,s,a-1,c,DP_P_BF0]
                                         + pop1[t-1,s,a-1,c,DP_P_BF7]
                                         + pop1[t-1,s,a-1,c,DP_P_BF12];
                      pop1[t,s,a,DP_CD4_350_500,DP_D_ARTgt12M] := pop1[t-1,s,a-1,c,DP_D_ARTlt6M]
                                                + pop1[t-1,s,a-1,c,DP_D_ART6to12M]
                                                + pop1[t-1,s,a-1,c,DP_D_ARTgt12M];
                      { Split 200-349 child group to 2 adult groups }
                      c := DP_CD4_Ped_200_349;
                      v1 := pop1[t-1,s,a-1,c,DP_P_Perinatal]
                            + pop1[t-1,s,a-1,c,DP_P_BF0]
                            + pop1[t-1,s,a-1,c,DP_P_BF7]
                            + pop1[t-1,s,a-1,c,DP_P_BF12];
                      pop1[t,s,a,DP_CD4_250_349,DP_D_HIVPos] := 2 / 3 * v1;
                      pop1[t,s,a,DP_CD4_200_249,DP_D_HIVPos] := 1 / 3 * v1;
                      v1 := pop1[t-1,s,a-1,c,DP_D_ARTlt6M]
                            + pop1[t-1,s,a-1,c,DP_D_ART6to12M]
                            + pop1[t-1,s,a-1,c,DP_D_ARTgt12M];
                      pop1[t,s,a,DP_CD4_250_349,DP_D_ARTgt12M] := 2 / 3 * v1;
                      pop1[t,s,a,DP_CD4_200_249,DP_D_ARTgt12M] := 1 / 3 * v1;
                     { Split <200 child group to 3 adult groups }
                      c := DP_CD4_Ped_LT200;
                      v1 := pop1[t-1,s,a-1,c,DP_P_Perinatal]
                            + pop1[t-1,s,a-1,c,DP_P_BF0]
                            + pop1[t-1,s,a-1,c,DP_P_BF7]
                            + pop1[t-1,s,a-1,c,DP_P_BF12];
                      pop1[t,s,a,DP_CD4_100_199,DP_D_HIVPos] := 0.35 * v1;
                      pop1[t,s,a,DP_CD4_50_99,DP_D_HIVPos] := 0.21 * v1;
                      pop1[t,s,a,DP_CD4_LT50,DP_D_HIVPos] := 0.44 * v1;
                      v1 := pop1[t-1,s,a-1,c,DP_D_ARTlt6M]
                            + pop1[t-1,s,a-1,c,DP_D_ART6to12M]
                            + pop1[t-1,s,a-1,c,DP_D_ARTgt12M];
                      pop1[t,s,a,DP_CD4_100_199,DP_D_ARTgt12M] := 0.35 * v1;
                      pop1[t,s,a,DP_CD4_50_99,DP_D_ARTgt12M] := 0.21 * v1;
                      pop1[t,s,a,DP_CD4_LT50,DP_D_ARTgt12M] := 0.44 * v1;

                      { Apply survival and migration to 15 year olds }
                      for c := DP_CD4_GT500 to DP_CD4_LT50 do
                      begin
                        pop1[t,s,a,c,DP_D_HIVPos] := pop1[t,s,a,c,DP_D_HIVPos] * (sr + mr * (1 + sr)/2);
                        No_ART := No_ART + pop1[t,s,a,c,DP_D_HIVPos];
                        pop1[t,s,a,c,DP_D_ARTgt12M] := pop1[t,s,a,c,DP_D_ARTgt12M]
                                           * (sr + mr * (1 + sr)/2);
                        On_ART := On_ART + pop1[t,s,a,c,DP_D_ARTgt12M];
                      end;

                      (* Save this for EPP - 9/25/14 *)
                      DP.SetHIVPos15YearOlds(t, No_ART);
                      DP.SetHIVPos15YearOldsART(t, On_ART);

                    end;
            16..80: begin
                      for c := DP_CD4_GT500 to DP_H_Max do
                      begin
                        pop1[t,s,a,c,DP_D_HIVPos] := pop1[t-1,s,a-1,c,DP_D_HIVPos] * sr1
                                                   + pop1[t-1,s,a-1,c,DP_D_HIVPos] * mr * (1 + sr1)/2;
                        if pop1[t,s,a,c,DP_D_HIVPos] < 0 then pop1[t,s,a,c,DP_D_HIVPos] := 0;
                        if a = 50 then No_ART_50 := No_ART_50 + pop1[t,s,a,c,DP_D_HIVPos];

                        for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
                        begin
                          pop1[t,s,a,c,d] := pop1[t-1,s,a-1,c,d] * sr1
                                           + pop1[t-1,s,a-1,c,d] * mr * (1 + sr1)/2;
                          if pop1[t,s,a,c,d] < 0 then pop1[t,s,a,c,d] := 0;
                          if a = 50 then On_ART_50 := On_ART_50 + pop1[t,s,a,c,d];
                        end;

                        if a = 50 then
                        begin
					  	            DP.SetHIVPos50YearOlds(t, No_ART_50);
                          DP.SetHIVPos50YearOldsART(t, On_ART_50);
                        end;
                      end;
                    end;
          end; { case a}
        end; {a loop}
        sr  := Get_Survival_Rate(DP, t, s, 80);
        sr1 := Get_Survival_Rate(DP, t, s, 81);
        mr := migr[80,t,s,DP_H_All] + migr[79,t,s,DP_H_All] / 2;
        if (pop1[t-1,s,79,DP_H_All,DP_D_All] + pop1[t-1,s,80,DP_H_All,DP_D_All]) <= 0
          then mr := 0
          else
            mr := mr / (pop1[t-1,s,79,DP_H_All,DP_D_All] + pop1[t-1,s,80,DP_H_All,DP_D_All]);
        dths := pop1[t-1,s,79,DP_H_All,DP_D_All] * (1 - sr)
              + pop1[t-1,s,79,DP_H_All,DP_D_All] * mr * (1 - sr) / 2
              + pop1[t-1,s,80,DP_H_All,DP_D_All] * (1 - sr1)
              + pop1[t-1,s,80,DP_H_All,DP_D_All] * mr * (1 - sr1) / 2;
       { Include COVID-19 deaths }
       if ((DP.GetFirstYear + t - 1) <= DP_COVID19FinalYr) and
          ((DP.GetFirstYear + t - 1) >= DP_COVID19FirstYr) then
       begin
       { Multiply total COVID deaths by percent to 80+ }
         temp := DP.GetCOVID19DeathRate(t, s, 0) *
                 DP.GetCOVID19DeathAgeDist(DP_Data, t, s, DP_A80_up) / 100 ;
         dths := dths + temp;
       end;

        DeathsByAge[80,s,t] := dths;
        DP.SetDeaths(t, DP_MAX_AGE, s, DP.GetDeaths(t, DP_Max_Age, s) + dths);

        for c := DP_H_No to DP_H_Max do
        begin
          if c = DP_H_No then d := DP_D_HIVNeg else d := DP_D_HIVPos;
          pop1[t,s,80,c,d] := pop1[t-1,s,79,c,d] * sr
                            + pop1[t-1,s,79,c,d] * mr * (1 + sr)/2
                            + pop1[t-1,s,80,c,d] * sr1
                            + pop1[t-1,s,80,c,d] * mr * (1 + sr1)/2;
          if pop1[t,s,80,c,d] < 0 then pop1[t,s,80,c,d] := 0;
          for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
          begin
            pop1[t,s,80,c,d] := pop1[t-1,s,79,c,d] * sr
                              + pop1[t-1,s,79,c,d] * mr * (1 + sr)/2
                              + pop1[t-1,s,80,c,d] * sr1
                              + pop1[t-1,s,80,c,d] * mr * (1 + sr1)/2;
            if pop1[t,s,80,c,d] < 0 then pop1[t,s,80,c,d] := 0;
          end;
        end;
      end; {s loop}

      // RLG: extra whitespace to preserve line numbering; safe to remove













      (* Sum Array to get all durations and status *)
      for a := 15 to 49 do
      begin
        Females[a] := 0;
        for c := DP_H_No to DP_H_Max do
          for d := DP_D_HIVNeg to DP_D_Max do
            Females[a] := Females[a]
                         + (pop1[t,DP_Female,a,c,d]
                            + pop1[t-1,DP_Female,a,c,d]) / 2;
      end; {a loop}

      { Set births for each fertile age }
      for s := DP_Female downto DP_Male do
      begin
        if s = DP_Female then
        begin
          for a := DP_AllAges to DP_A45_49 do
            DP.SetBirths(t, a, DP_BothSexes, 0);
          for a := 15 to 49 do
          begin
            a5 := trunc(a/5) + 1;
            { Calculate births }
            if DP.GetASFR (t,DP_AllAges) <> 0 then
              Tempbirths := DP.GetTFR(t) * Females[a] *
                            DP.GetASFR(t, a5) / DP.GetASFR(t, DP_AllAges) / 5
            else
            begin
              TempBirths := 0;
              DP_AddLogEntry(GB, DP_DPAIM,'DPPROJ1','~3-30','DP.GetASFR (t,DP_AllAges)',
                           DP_EZ,DP.GetASFR (t,DP_AllAges), 't = ' + IntToStr(t));
            end;

            DP.SetBirths(t, DP_AllAges, DP_BothSexes,
                        DP.GetBirths(t, DP_AllAges,DP_BothSexes) + Tempbirths);
            DP.SetBirths(t, a5, DP_BothSexes,
                        DP.GetBirths(t, a5,DP_BothSexes) + TempBirths);

          end; { a loop }
        end; { s = Female }
        temp := GetDP_SexBirthRatioPercent(DP,t,s);

        DP.SetBirths(t, DP_AllAges, s, DP.GetBirths(t, DP_AllAges,DP_BothSexes) *
                    temp);
        sr := Get_Survival_Rate(DP, t, s, 0);
        mr := migr[0,t,s,DP_H_All] / 2;
        if DP.GetBirths(t, DP_AllAges, s) <= 0
          then mr := 0
          else mr := mr / DP.GetBirths(t, DP_AllAges, s);
        dths := DP.GetBirths(t, DP_AllAges,s) * (1 - sr)
              + DP.GetBirths(t, DP_AllAges,s) * mr * (2 * (1 - sr))/3;
        DeathsByAge[0,s,t] := dths;
        DP.SetDeaths(t, DP_A0_4, s, DP.GetDeaths(t, DP_A0_4, s) + dths);
        pop1[t,s,0,DP_H_No,DP_D_HIVNeg] := DP.GetBirths(t, DP_AllAges,s) * sr
                               + DP.GetBirths(t, DP_AllAges,s) * mr * (1 + 2 * sr)/3;
      end; { s loop }

      (* Sum Array to get all durations and status *)
      SumPopulationByHIVStatus(t);

       (*Recording AIDS deaths by CD4 count for orphan calculations*)
      for s := DP_Male to DP_Female do
        for a := 15 to DP_MaxSingleAges do
          for c  := DP_CD4_GT500 to DP_CD4_LT50 do
          begin
            a10 := DP_cd4_15_24;
            case a of
               15 .. 24:
              a10 := DP_CD4_15_24;
              25 .. 34:
              a10 := DP_CD4_25_34;
              35 .. 44:
              a10 := DP_CD4_35_44;
              45 .. DP_MaxSingleAges:
              a10 := DP_CD4_45_54;
             end;
            mu[a10, c, s] := DP.GetAdultMortByCD4NoART(DP_Data, s, a10, c);
            d:= DP_D_HIVPos;

            DP.SetAIDSDeathsCD4NoART(t,a+1,s,c, d, mu[a10, c, s] * pop1[t, s, a, c, d])   ;
            DP.SetAIDSDeathsCD4NoART(t,DP_AllSingAges,s,c,d, DP.GetAIDSDeathsCD4NoART(t,DP_AllSingAges,s,c,d)
                                                             +DP.GetAIDSDeathsCD4NoART(t,a+1,s,c,d));

            alpha[a10, c, s, DP_D_ARTlt6m]   := DP.GetAdultMortByCD4WithART0to6(DP_Data, s, a10, c);
            alpha[a10, c, s, DP_D_ART6to12m] := DP.GetAdultMortByCD4WithART7to12(DP_Data, s, a10, c);
            alpha[a10, c, s, DP_D_ARTgt12m]  := DP.GetAdultMortByCD4WithARTGT12(DP_Data, s, a10, c);


            DP.SetAIDSDeathsCD4OnART(t,a+1,s,c,DP_D_All, alpha[a10,c,s,DP_D_ARTlt6m]   * pop1[t,s,a,c,DP_D_ARTlt6m]
                                                       + alpha[a10,c,s,DP_D_ARTgt12m]  * pop1[t,s,a,c,DP_D_ARTgt12m]
                                                       + alpha[a10,c,s,DP_D_ART6to12m] * pop1[t,s,a,c,DP_D_ART6to12m]);

            DP.SetAIDSDeathsCD4OnART(t,DP_AllSingAges,s,c, DP_D_All, DP.GetAIDSDeathsCD4OnART(t,DP_AllSingAges,s,c,DP_D_All)
                                                              + DP.GetAIDSDeathsCD4OnART(t,a,s,c,DP_D_All));
            for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
            begin
              if d = DP_D_ARTgt12m then
                d1 := DP_MortRates_GT12Mo
              else
                d1 := DP_MortRates_LT12Mo;
              DP.SetAIDSDeathsCD4OnART(t,a+1,s,c,d,  alpha[a10, c, s, d]
                * (DP.GetMortalityRates(DP_Data,d1,t) * DP.GetMortalityRatesMultiplier)
                * pop1[t, s, a, c, d]);
              DP.SetAIDSDeathsCD4OnART(t,DP_AllSingAges,s,c, d, DP.GetAIDSDeathsCD4OnART(t,DP_AllSingAges,s,c,d)
                                                                + DP.GetAIDSDeathsCD4OnART(t,a+1,s,c,d));
            end;
          end;   //end s,c,a

      for s := DP_Male to DP_Female do
        for c  := DP_CD4_GT500 to DP_CD4_LT50 do
        begin
          DP.SetAIDSDeathsCD4NoART(t,DP_A80,s,c,DP_D_HIVPos, DP.GetAIDSDeathsCD4NoART(t,DP_A79,s,c,DP_D_HIVPos));
          DP.SetAIDSDeathsCD4NoART(t,DP_AllSingAges,s,c,DP_D_HIVPos, DP.GetAIDSDeathsCD4NoART(t,DP_AllSingAges,s,c,DP_D_HIVPos)
                                                              + DP.GetAIDSDeathsCD4NoART(t,DP_A80,s,c,DP_D_HIVPos));
          for d in [DP_D_All, DP_D_ARTlt6m, DP_D_ART6to12m, DP_D_ARTgt12m] do
          begin
            DP.SetAIDSDeathsCD4OnART(t,DP_A80,s,c,d, DP.GetAIDSDeathsCD4OnART(t,DP_A79,s,c,d)) ;
            DP.SetAIDSDeathsCD4OnART(t,DP_AllSingAges,s,c,d, DP.GetAIDSDeathsCD4OnART(t,DP_AllSingAges,s,c,d)
                                                                +DP.GetAIDSDeathsCD4OnART(t,DP_A80,s,c,d));
          end;
        end;

      (*Recording HIV infected children for orphan calculations*)
      for a := 0 to 4 do
      begin
        temp := 0;
        for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
          for d := DP_P_Perinatal to DP_D_Max do
            temp := temp + pop1[t,DP_Female,a,c,d] + pop1[t,DP_Male,a,c,d];

        DP.SetHIVChildSingleAge(t,a+1, temp);
        DP.SetHIVChildSingleAge(t,DP_AllSingAges, DP.GetHIVChildSingleAge(t,DP_AllSingAges)
                                                + DP.GetHIVChildSingleAge(t,a+1));
      end;

      for a := 5 to 14 do
      begin
        temp := 0;
        for c := DP_CD4_Ped_GT1000 to DP_CD4_Ped_LT200 do
          for d := DP_P_Perinatal to DP_D_Max do
            temp := temp + pop1[t,DP_Female,a,c,d] + pop1[t,DP_Male,a,c,d];

        DP.SetHIVChildSingleAge(t,a+1, temp);
        DP.SetHIVChildSingleAge(t,DP_AllSingAges, DP.GetHIVChildSingleAge(t,DP_AllSingAges)
                                                + DP.GetHIVChildSingleAge(t,a+1));
      end;

      for a:= 15 to 18 do
      begin
        DP.SetHIVChildSingleAge(t,a+1, pop1[t,DP_Female,a,DP_H_All,DP_D_All] - pop1[t,DP_Female,a,DP_H_No,DP_D_All]
                                + pop1[t,DP_Male,a,DP_H_All,DP_D_All] - pop1[t,DP_Male,a,DP_H_No,DP_D_All]);
        DP.SetHIVChildSingleAge(t,DP_AllSingAges, DP.GetHIVChildSingleAge(t,DP_AllSingAges)
                                                + DP.GetHIVChildSingleAge(t,a+1));
      end;

      { Perform AIDS calculations if necessary }
      if (GetDP_IncludeAIDS(GB, DP)) or (p = GB.GetTempProj) then
        Do_AIDS_Calculations(p, t);

      (* Sum Array to get all durations and status *)
      SumPopulationByHIVStatus(t);

      (* If Regional population adjustment is being used *)
      if (DP.getRegionalAdjustPopCBState) and ((DP.getCustomPopStopRescalingYear -DP.GetFirstYear + 1) >= t) then
      begin
        for s := DP_Female downto DP_Male do
        begin
          DP.SetPopAdj(s,t,0);

          for a5 := DP_AllAges to DP_MAX_AGE do DP.SetPopAdjAmount(t,s,a5,0);

          for a := 0 to DP_MaxSingleAges do
          begin
           { Calculate ratio of population by age and sex from external file to Spectrum calculation for all HIV states and all durations}
           (* ratio := VALUE READ FROM FILE FOR POP(AGE=A, SEX=S, YEAR = T) / pop1[t,s,a, DP_H_All, DP_All] *)
           (* need to use a+1 in the RegionalAdjustPopData variable since the values go from AllSingleAges(0) to
              DP_A_80(81) *)
            a5 := trunc(a/5) + 1;

            DP.SetPopAdjAmount(t, s, a5, DP.GetPopAdjAmount(t,s,a5) +
                           DP.GetRegionalAdjustPopData(s,a+1,t) - pop1[t,s,a,DP_H_All, DP_D_All]);

            ratio := DP.GetRegionalAdjustPopData(s,a+1,t)/pop1[t,s,a,DP_H_All, DP_D_All];
           { Adjust all population values by ratio }
            for c := DP_H_No to DP_H_Max do
              for d := DP_D_HIVNeg to DP_D_Max do
              begin
                pop1[t,s,a,c,d] := pop1[t,s,a,c,d] * ratio;

                if (a >= 15) and (a <= 49) then
                  DP.SetPopAdj(s,t, DP.GetPopAdj(s,t) + (pop1[t,s,a,c,d] * (1 - ratio)));
              end;
          end;{a loop}
        end;

        (* Sum Array to get all durations and status *)
        SumPopulationByHIVStatus(t);
      end;

      (* if this is a subnational projection *)
      if (GB.GetSubnatRegionName(p) <> '') and (DP.GetPercentOfPop(t) > 0) and (not DP.GetRegionalAdjustPopCBState) then
      begin
        RegTot := 0;
        for a := 0 to DP_MaxSingleAges do
          for s := DP_Male to DP_Female do
            for c := DP_H_No to DP_H_Max do
              for d := DP_D_HIVNeg to DP_D_Max do
                RegTot := RegTot + pop1[t,s,a,c,d];

        (* if DP.GetCountryProjBigPop = 0 then that means that a country was not selected
           for the region, so we need to just keep the population the same as the regions *)
        if DP.GetCountryProjBigPop(t) = 0 then
          ratio := 1
        else
          ratio := (DP.GetCountryProjBigPop(t) * DP.GetPercentOfPop(t)/100) / RegTot;

        for a := 0 to DP_MaxSingleAges do
          for s := DP_Male to DP_Female do
            for c := DP_H_All to DP_H_Max do
              for d := DP_D_All to DP_D_Max do
                pop1[t,s,a,c,d] := pop1[t,s,a,c,d] * ratio;
      end;

      { Fill summary array indices }
      Cumulate_Pop_Totals(DP, t);

      pr := 0.85*(t/DP.GetFinalIndex);
      GB.SetProgressPanelPercent (Round (pr*100));

      { Calculate AIDS adjusted life expectancy }
      LEcalc(GB, DP,t, AidsDeathsByAge,pop1);

      (* If OHT is active, save the data in the AidsDeathsByAge and pop1 arrays
         as they are necessary for the calculation of OHT's Healthy life
         expectancy result. We can't do the Healthy life expectancy calculation
         here because it depends on results from LiST and NCD which have not
         been calculated yet at this point. *)
      if GB_OHT_On(GB, p) then
      begin
        for s := DP_Male to DP_Female do
          for a  := 0 to DP_MaxSingleAges do
            DP_SetUHAidsDeathsByAge(p, a, s, t, AidsDeathsByAge[a, s, t]);

        for s := DP_BothSexes to DP_Female do
          for a := 0 to DP_MaxSingleAges do
            for c := DP_H_All to DP_H_Max do
              for d := DP_D_All to DP_D_Max do
                DP_SetUHPop1(p, t, s, a, c, d, pop1[t,s,a,c,d]);
      end;

      (* Capture pop at one year prior to calc year and save to file.  This procedure
       is only called when the projection is calculating from the first year
       The pop is also captured when we are regenerating CalcState because the
       user changed the first year population. We do not what to capture CalcState for
       projections that begin prior to 1986 so we need to check the TempProj *)
      if (GB.GetTempProj > 1) and (t = DP.GetFinalIndex) then
      begin
        Idx := GBTempProjIdx(GB, GB_FirstYrIdxOfPrimProj);
        (* Need to set DPAIDSDeathsBySingleAge since it is needed by calcstate*)
        for s := DP_Male to DP_Female do
          for a := 0 to DP_MaxSingleAges do
          begin
            DP.SetDeathsBySingleAge(Idx, a, s, DeathsByAge[a,s,Idx]);
            DP.SetAIDSDeathsBySingleAge(Idx, a, s, AidsDeathsByAge[a,s,Idx]);
          end;
        GBSession.GetProjection(GB.GetPrimaryProj).DP.SetIncidenceAdjustmentFactor2012(DP.GetIncidenceAdjustmentFactor2012);
        CaptureCalcStateData(pop1);
//        if GB.ModuleActive(GB_HV) then
//          HV_CaptureCalcStateData(p);
      end;

      (* If we are in CalcState mode, we may (if the Insert users pop CB is checked in options)
         want to replace the first year population with the values contained in the First Year
         population editor prior to calculting the projection. We need to replace the first
         year population just in case it was changed by the user, because if it was changed
         a BackCast has occurred. When we BackCast, the values will be off just a bit especially
         if the values for each age were not adjusted by the same percent. We want the values
         displayed in the editor to be exactly the same as those entered by the user. *)
      if GB.GetInsertUsersPopBoolean and GB.GetInCalcStateMode(p) and (t = 2) then
      begin
        for a := 0 to DP_MaxSingleAges do
        begin
          for s := DP_Male to DP_Female do
          begin
            (* We need to know the adjustment made to the calculated population so we can
               apply the same adjustment to all HIV levels. We also need to make the same
               adjustment to the population at time t = 1 so that births and deaths can be
               calculated properly (births, for example, uses the current year and prior year
               divided by 2). The value in the pop1 array at this point is the calculated value
               for t = 2 and BigPop has the user entered value from the first year population editor *)
            if pop1[t,s,a,DP_H_All,DP_D_All] = 0 then
              PopAdj := 1
            else
              PopAdj := DP.GetBigPop(t,s,a) / pop1[t,s,a,DP_H_All,DP_D_All];

            (* Adjust all levels of HIV by the same adjustment that was done for the
               first year population *)
            for c := DP_H_All to DP_H_Max do
              for d := DP_D_All to DP_D_Max do
                pop1[t,s,a,c,d] := pop1[t,s,a,c,d] * PopAdj;

            (* set the user entered value from the first year pop editor into
               the pop1 array for t = 2 *)
            pop1[t,s,a,DP_H_All,DP_D_All] := DP.GetBigPop(t,s,a);
          end;

          (* set the values for bothsexes *)
          for c := DP_H_All to DP_H_Max do
            for d := DP_D_All to DP_D_Max do
              pop1[t,DP_BothSexes,a,c,d] := pop1[t,DP_Male,a,c,d] + pop1[t,DP_Female,a,c,d];
        end;
      end;

      (* Save the population at time t that was just calculated *)
      for a := 0 to DP_MaxSingleAges do
        for s := DP_BothSexes to DP_Female do
        begin
          if s = DP_BothSexes then
            DP.SetBigPop(t, DP_BothSexes, a,
                        pop1[t,DP_Male,a,DP_H_All,DP_D_All]+pop1[t,DP_Female,a,DP_H_All,DP_D_All])
          else
            DP.SetBigPop(t, s, a, pop1[t,s,a,DP_H_All,DP_D_All]);
        end; {s loop}

      { Calculate HIV-exposed but uninfected (HEU) children }
      { Calculated as sum over previous 14 years to capture all children
        aged 0-14 at time t
        New HEU are children born to HIV+ mothers minus those who became infected
        New HEU are survived to time t and cumulated over last 14 years }
      for s := DP_Male to DP_Female do
      begin
        v1 := 0;
        d := max(1,t-14);
        for t1 := d to t do
        begin
          temp := DP.GetNeedPMTCT(t1) * GetDP_SexBirthRatioPercent(DP,t1,s)
                  - DP.GetNewInfections(t1,DP_A0_4,s);
          for a:= 0 to t-t1 do
            temp := temp * Get_Survival_Rate(DP,t,s,a);
          v1 := v1 + temp;
        end;
        DP.SetHIVExpUninfChildren(s, t, v1);
      end;
      DP.SetHIVExpUninfChildren(DP_BothSexes, t,
                                    DP.GetHIVExpUninfChildren(DP_Male, t)
                                  + DP.GetHIVExpUninfChildren(DP_Female, t));

      { Calculate ART-EU children }
      { Calculated as children born to mother who received PMTCT
        adjusted for the percentage who became infected (PTR2) }
      for s := DP_Male to DP_Female do
      begin
        v1 := 0;
        d := max(1,t-14);
        for t1 := d to t do
        begin
          temp := DP.GetOnPMTCT(t1)
                  * GetDP_SexBirthRatioPercent(DP,t1,s)
                  * (1 - PTR2[t1]/100);
          for a:= 0 to t-t1 do
            temp := temp * Get_Survival_Rate(DP,t,s,a);
          v1 := v1 + temp;
        end;
        DP.SetARTExpUninfChildren(s, t, v1);
      end;
      DP.SetARTExpUninfChildren(DP_BothSexes, t,
                                    DP.GetARTExpUninfChildren(DP_Male, t)
                                  + DP.GetARTExpUninfChildren(DP_Female, t));

    end; { t loop }

//    (*Carel: ART info needed by TB module*)
//    for t := 1 to DP.GetFinalIndex do
//    begin
//      temp := 0;
//      temp_all := 0;
//      for c := DP_CD4_GT500 to DP_CD4_LT50 do
//        for a := 15 to DP_MaxSingleAges do
//        begin
//          temp := temp + pop1[t,DP_Male,a,c,DP_D_ARTgt12M] + pop1[t,DP_Female,a,c,DP_D_ARTgt12M];
//          for d := DP_D_ARTlt6M to DP_D_ARTgt12M do
//            temp_all := temp_all + pop1[t,DP_Male,a,c,d] + pop1[t,DP_Female,a,c,d];
//        end;{for c, a}
//
//        if(temp_all > 0 ) then
//          temp := temp / temp_all * 100
//        else
//          temp := 0;
//        DP.SetOnARTgt12M(t, temp);
//    end;{t loop}


    (* if Goals is active, run Goals PostOuterLoopCalcs;*)
    if (GB.ModuleActive(p, GB_HV)) then
    begin
      HVObj.PostOuterLoopCalcs;
      if not GB.GetInCalcStateMode(p) then
        (* at t=_EpidemicStYr-1, set new infections=initial pulse to zero *)
        if GB.GetFinalYear(p) >= (GB.GetCalcYear(p) + DP_GetHVEpidemicStYr(p) - GB.GetCalcYearIdx(p)) then
        begin
          t := DP_GetHVEpidemicStYr(p)-1;
          for s := DP_BothSexes to DP_Female do
             for a := DP_A0_4 to DP_A80_up do
                 DP.SetNewInfections(t,a,s,0);
      end;
    end;

    (* Save pop by single age *)
    for a := 0 to DP_MaxSingleAges do
      for t := 1 to DP.GetFinalIndex do
        for s := DP_BothSexes to DP_Female do
        begin
          if s = DP_BothSexes then
            DP.SetBigPop(t, DP_BothSexes, a,
                        pop1[t,DP_Male,a,DP_H_All,DP_D_All]+pop1[t,DP_Female,a,DP_H_All,DP_D_All])
          else
            DP.SetBigPop(t, s, a, pop1[t,s,a,DP_H_All,DP_D_All]);
        end; {s loop}

    (* Save PLHIV on ART by single age *)
    for t := 1 to DP.GetFinalIndex do
      for a := DP_A0 to DP_A80 do
      begin
        for s := DP_Male to DP_Female do
        begin
          SumPLHIVOnART := 0;
          if a <= DP_A4 then
          begin
            for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
              SumPLHIVOnART := SumPLHIVOnART + pop1[t,s,a-1,c,DP_D_ARTLt6m] +
                               pop1[t,s,a-1,c,DP_D_ART6to12m] + pop1[t,s,a-1,c,DP_D_ARTgt12m];
          end
          else if a <= DP_A14 then
          begin
            for c := DP_CD4_PED_GT1000 to DP_CD4_PED_LT200 do
              SumPLHIVOnART := SumPLHIVOnART + pop1[t,s,a-1,c,DP_D_ARTLt6m] +
                               pop1[t,s,a-1,c,DP_D_ART6to12m] + pop1[t,s,a-1,c,DP_D_ARTgt12m];
          end
          else
          begin
            for c := DP_CD4_GT500 to DP_CD4_LT50 do
              SumPLHIVOnART := SumPLHIVOnART + pop1[t,s,a-1,c,DP_D_ARTLt6m] +
                               pop1[t,s,a-1,c,DP_D_ART6to12m] + pop1[t,s,a-1,c,DP_D_ARTgt12m];
          end;

          DP.SetOnARTBySingleAge(t, s, a, SumPLHIVOnART);
        end; {s loop}
      end;

    (* Save HIV+ by single age for adults *)
    for t := 1 to DP.GetFinalIndex do
      for a := 0 to DP_MaxSingleAges do
      begin
        for s := DP_Male to DP_Female do
        begin
          Temp := 0;
          if a < 15 then
          begin
            for d := DP_P_Perinatal to DP_P_BF12 do
            begin
              if a < 5 then
                for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
                  Temp := Temp + pop1[t,s,a,c,d]
              else
                for c := DP_CD4_Ped_GT1000 to DP_CD4_Ped_LT200 do
                  Temp := Temp + pop1[t,s,a,c,d]
            end;
            for d := DP_D_ARTlt6M to DP_D_ARTgt12m do
            begin
              if a < 5 then
                for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
                  Temp := Temp + pop1[t,s,a,c,d]
              else
                for c := DP_CD4_Ped_GT1000 to DP_CD4_Ped_LT200 do
                  Temp := Temp + pop1[t,s,a,c,d]
            end;
          end
          else
          begin
            for c := DP_CD4_GT500 to DP_H_Max do
              Temp := Temp + pop1[t,s,a,c,DP_D_All];
          end;
          DP.SetHIVBySingleAge(t,a,s,Temp);
        end;
        DP.SetHIVBySingleAge(t, a, GB_BothSexes, DP.GetHIVBySingleAge(t, a, GB_Male) +
          DP.GetHIVBySingleAge(t, a, GB_Female));
      end;

     (* Save Deaths by single age, currently used for CS calculations *)
     for t := 1 to DP.GetFinalIndex do
        for s := DP_Male to DP_Female do
          for a := 0 to DP_MaxSingleAges do
          begin
            DP.SetDeathsBySingleAge(t, a, s, DeathsByAge[a,s,t]);
            DP.SetAIDSDeathsBySingleAge(t, a, s, AidsDeathsByAge[a,s,t]);
          end; {a loop}

    (* CD4 Distribution for 15+ *)
    for t := 1 to DP.GetFinalIndex do
    begin
      for c := DP_CD4_GT500 to DP_CD4_LT50 do
      begin
        tempM := 0;
        tempF := 0;
        for a := 15 to 80 do tempM := tempM + pop1[t,DP_Male,  a,c,DP_D_HIVPos];
        for a := 15 to 80 do tempF := tempF + pop1[t,DP_Female,a,c,DP_D_HIVPos];
        DP.SetCD4Distribution(t, DP_BothSexes, c, DP_NoTreat, tempM + tempF);
        DP.SetCD4Distribution(t, DP_Male,      c, DP_NoTreat, tempM);
        DP.SetCD4Distribution(t, DP_Female,    c, DP_NoTreat, tempF);
      end; {c loop}

      for c := DP_CD4_GT500 to DP_CD4_LT50 do
      begin
        tempM := 0;
        tempF := 0;
        for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
        begin
          tempMDur := 0;
          tempFDur := 0;
          for a := 15 to 80 do tempMDur := tempMDur + pop1[t,DP_Male,  a,c,d];
          for a := 15 to 80 do tempFDur := tempFDur + pop1[t,DP_Female,a,c,d];
          tempM := tempM + tempMDur;
          tempF := tempF + tempFDur;

          d1 := d - DP_D_ARTlt6m + 1; // Map d=6..8 to d1=1..3 for the TI module
          DP.SetCD4DistributionDur(t, DP_BothSexes, c, d1, tempMDur + tempFDur);
          DP.SetCD4DistributionDur(t, DP_Male,      c, d1, tempMDur);
          DP.SetCD4DistributionDur(t, DP_Female,    c, d1, tempFDur);
        end;

        DP.SetCD4Distribution(t, DP_BothSexes, c, DP_OnART, tempM + tempF);
        DP.SetCD4Distribution(t, DP_Male,      c, DP_OnART, tempM);
        DP.SetCD4Distribution(t, DP_Female,    c, DP_OnART, tempF);
      end; {c loop}
    end; {t loop}

    (* CD4 distribution for 15-49 *)
    for t := 1 to DP.GetFinalIndex do
    begin
      for c := DP_CD4_GT500 to DP_CD4_LT50 do
      begin
        tempM := 0;
        tempF := 0;
        for a := 15 to 49 do tempM := tempM + pop1[t,DP_Male,  a,c,DP_D_HIVPos];
        for a := 15 to 49 do tempF := tempF + pop1[t,DP_Female,a,c,DP_D_HIVPos];
        DP.SetCD4Distribution15_49(t, DP_BothSexes, c, DP_NoTreat, tempM + tempF);
        DP.SetCD4Distribution15_49(t, DP_Male,      c, DP_NoTreat, tempM);
        DP.SetCD4Distribution15_49(t, DP_Female,    c, DP_NoTreat, tempF);
      end; {c loop}

      for c := DP_CD4_GT500 to DP_CD4_LT50 do
      begin
        tempM := 0;
        tempF := 0;
        for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
        begin
          for a := 15 to 49 do tempM := tempM + pop1[t,DP_Male,  a,c,d];
          for a := 15 to 49 do tempF := tempF + pop1[t,DP_Female,a,c,d];
        end;
        DP.SetCD4Distribution15_49(t, DP_BothSexes, c, DP_OnART, tempM + tempF);
        DP.SetCD4Distribution15_49(t, DP_Male,      c, DP_OnART, tempM);
        DP.SetCD4Distribution15_49(t, DP_Female,    c, DP_OnART, tempF);
      end; {c loop}
    end; {t loop}
  finally
    for s := DP_Male to DP_Female do
    for a  := 0 to DP_MaxSingleAges do
    begin
      finalize(AidsDeathsByAge[a,s]);
      finalize(DeathsByAge[a,s]);
    end;
  end;

  (* CD4 percent distribution for those not on ART and aged 0-4 *)
  for t := 1 to DP.GetFinalIndex do
  begin
    for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
    begin
      for d := DP_P_Perinatal to DP_P_BF12 do
      begin
        tempM := 0;
        tempF := 0;
        for a := 0 to 4 do
        begin
          tempM := tempM + pop1[t,DP_Male,  a,c,d];
          tempF := tempF + pop1[t,DP_Female,a,c,d];
        end;
        DP.SetCD4DistributionNoARTChild(t,DP_BothSexes,DP_CD4_0t4,c,d,tempM + tempF);
        DP.SetCD4DistributionNoARTChild(t,DP_Male,     DP_CD4_0t4,c,d,tempM);
        DP.SetCD4DistributionNoARTChild(t,DP_Female,   DP_CD4_0t4,c,d,tempF);
      end;
    end; {c loop}
  end; {t loop}

  (* CD4 count distribution for those not on ART and aged 5-14 *)
  for t := 1 to DP.GetFinalIndex do
  begin
    for c := DP_CD4_Ped_GT1000 to DP_CD4_Ped_LT200 do
    begin
      for d := DP_P_Perinatal to DP_P_BF12 do
      begin
        tempM := 0;
        tempF := 0;
        for a := 5 to 14 do
        begin
          tempM := tempM + pop1[t,DP_Male,  a,c,d];
          tempF := tempF + pop1[t,DP_Female,a,c,d];
        end;
        DP.SetCD4DistributionNoARTChild(t,DP_BothSexes,DP_CD4_5t14,c,d,tempM + tempF);
        DP.SetCD4DistributionNoARTChild(t,DP_Male,     DP_CD4_5t14,c,d,tempM);
        DP.SetCD4DistributionNoARTChild(t,DP_Female,   DP_CD4_5t14,c,d,tempF);
      end;
    end; {c loop}
  end; {t loop}

  (* CD4 percent distribution for those on ART and aged 0-4 *)
  for t := 1 to DP.GetFinalIndex do
  begin
    for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
    begin
      for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
      begin
        tempM := 0;
        tempF := 0;
        for a := 0 to 4 do
        begin
          tempM := tempM + pop1[t,DP_Male,  a,c,d];
          tempF := tempF + pop1[t,DP_Female,a,c,d];
        end;
        DP.SetCD4DistributionWithARTChild(t,DP_BothSexes,DP_CD4_0t4,c,d,tempM + tempF);
        DP.SetCD4DistributionWithARTChild(t,DP_Male,     DP_CD4_0t4,c,d,tempM);
        DP.SetCD4DistributionWithARTChild(t,DP_Female,   DP_CD4_0t4,c,d,tempF);
      end;
    end; {c loop}
  end; {t loop}

  (* CD4 count distribution for those not on ART and aged 5-14 *)
  for t := 1 to DP.GetFinalIndex do
  begin
    for c := DP_CD4_Ped_GT1000 to DP_CD4_Ped_LT200 do
    begin
      for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
      begin
        tempM := 0;
        tempF := 0;
        for a := 5 to 14 do
        begin
          tempM := tempM + pop1[t,DP_Male,  a,c,d];
          tempF := tempF + pop1[t,DP_Female,a,c,d];
        end;
        DP.SetCD4DistributionWithARTChild(t,DP_BothSexes,DP_CD4_5t14,c,d,tempM + tempF);
        DP.SetCD4DistributionWithARTChild(t,DP_Male,     DP_CD4_5t14,c,d,tempM);
        DP.SetCD4DistributionWithARTChild(t,DP_Female,   DP_CD4_5t14,c,d,tempF);
      end;
    end; {c loop}
  end; {t loop}

  // RLG: extra whitespace to preserve line numbering; safe to remove













  for a := 0 to DP_MaxSingleAges do
    for t := 1 to DP.GetFinalIndex do
    begin
      value := GetDP_LECalcAge(GB, DP, t, a, pop1);                 //conversion, remove p
      DP.SetLESingleAge(a, t, value);
    end;

  { Calculate base year births and deaths }
{  Calc_Base_Year_Births_and_Deaths(p);}

  (* Calculate values for HIV pop by risk group and New HIV infections by risk group.  The values for Prev and Inc
     are taken from the SPT file for each risk group (region) *)
  for r := 1 to DP.GetNumOfEPPEpidemics do
    for t := 1 to DP.GetFinalIndex do
    begin
      DP.SetHIVPopByRiskGroup(r,t,DP.GetEPPPopData(r,t) * DP.GetEppPrevData(r,t)/100);
      if DP.GetEPPIncData(r,t) <= 0 then
        DP.SetNewHIVInfByRiskGroup(r,t,0.0)
      else
        DP.SetNewHIVInfByRiskGroup(r,t,(DP.GetEPPPopData(r,t) * (1 - (DP.GetEppPrevData(r,t)/100)))
                                                                        * (DP.GetEPPIncData(r,t)/100));
    end;

  (* Adjust the EPP values for HIV pop by risk group to match those from Spectrum *)
  for t := 1 to DP.GetFinalIndex do
  begin
    SumSpec := 0;
    for a := DP_A15_19 to DP_A45_49 do
      for s := DP_Male to DP_Female do
        SumSpec := SumSpec + DP.GetHIV(t, s, a);

    SumEPP := 0;
    for r := 1 to DP.GetNumOfEPPEpidemics do
      SumEPP := SumEPP + DP.GetHIVPopByRiskGroup(r,t);

    if SumEPP = 0 then
      AdjustFactor := 0
    else
      AdjustFactor := SumSpec/SumEPP;

    for r := 1 to DP.GetNumOfEPPEpidemics do
      DP.SetHIVPopByRiskGroup(r,t, AdjustFactor * DP.GetHIVPopByRiskGroup(r,t));
  end;

  (* Adjust the EPP values for new HIV infections by risk group to match those from Spectrum *)
  for t := 1 to DP.GetFinalIndex do
  begin
    SumSpec := 0;
    for a := DP_A15_19 to DP_A45_49 do
      for s := DP_Male to DP_Female do
        SumSpec := SumSpec + DP.GetNewInfections(t,a,s);

    SumEPP := 0;
    for r := 1 to DP.GetNumOfEPPEpidemics do
      SumEPP := SumEPP + DP.GetNewHIVInfByRiskGroup(r,t);

    if SumEPP = 0 then
      AdjustFactor := 0
    else
      AdjustFactor := SumSpec/SumEPP;

    for r := 1 to DP.GetNumOfEPPEpidemics do
      DP.SetNewHIVInfByRiskGroup(r,t, AdjustFactor * DP.GetNewHIVInfByRiskGroup(r,t));
  end;

  if GB.ModuleActive(p, GB_HA) and (not DP.GetIsFitting) then
  begin
    HA_CalculateCosts(p);
  end;

  if assigned(HVObj) then
    HVObj.Free;
  if assigned(HAObj) then
    HAObj.Free;
  GB.SetProgressPanelPercent (0);

   for s := DP_Male to DP_Female do                      //for .EP5 file
      for c := DP_CD4_Ped_GT1000 to DP_CD4_Ped_LT200 do
        for d := DP_P_Perinatal to DP_GT12monthsTreat do
          for t := GB.GetCalcYearIdx(p) to GB.GetFinalYearIdx(p) do
            DP.SetChAged14ByCD4Cat(t, s, c, d, pop1[t, s, 14, c, d]);

  if (not DP.GetIsFitting) then
  begin
    Calc_NewInfectionsLessImmigrants(DP);
    Calc_BothSexes;
    Calc_ExtraIndicators;
  end;

  if GB.GetDumpPop1Array then
    DumpPop1ToFile(GB, p, GetGBVariables.GetPrjPath(p), pop1);

//  stopwatch.stop;
//  OutputDebugString(PChar(Format('ProjectPopulationBySingleYears completed in %d ms', [stopwatch.ElapsedMilliseconds])));
end;

procedure TDPProj.CalcOrphans ;
var
  //f               : textfile;
  a, s, c         : Byte;
  //dummyMat        : array[1..10000] of integer;
  tempVal         : Double;
  //d, dd, value,
  t               : Integer;
  sum             : double;
begin
  for t := GB.GetCalcYearIdx(p)-DP_MaxChildAge to DP.GetFinalIndex do
    for s:= DP_Orphans_Maternal to DP_orphans_Dual do
      for c := DP_Orphans_Total to DP_Orphans_NON_AIDS do
        for a := DP_AllSingAges to DP_MaxChildAge do
        begin
          DP.SetOrphansTotal (t,s,c,0,a, 0.0);
          DP.SetOrphansTotal (t,s,c,1,a,0.0);
          DP.SetOrphansNew  (t,s,c,0,a, 0.0);
          DP.SetOrphansNew  (t,s,c,1,a, 0.0);
        end;

  MaternalOrphanCalcs();
  PaternalOrphanCalcs();

 (* ALl new orphan calculations have been estimated into 2 arrays DP.GetOrphansTotal and DP.GetOrphansNew,
 but the older orphan calculations had several array.
 So SET PREVIOUS ORPHAN CALCULATION ARRAY WITH THESE NEW ORPHAN CALCS*)

  (*.............TOTAL ORPHANS.............*)
  for t := GB.GetCalcYearIdx(p) to DP.GetFinalIndex do
  begin
    for a := DP_A0  to DP_MaxChildAge do
    begin
      DP.SetAgeLastBDay(t,DP_MaternalAIDS_ALB,a-1, (DP.GetOrphansTotal(t,DP_Orphans_Maternal,DP_Orphans_AIDS,0,a)
                                                +DP.GetOrphansTotal(t,DP_Orphans_Maternal,DP_Orphans_AIDS,1,a)
                                                +DP.GetOrphansTotal(t,DP_Orphans_Maternal,DP_Orphans_AIDS,0,max(DP_A0,a-1))
                                                +DP.GetOrphansTotal(t,DP_Orphans_Maternal,DP_Orphans_AIDS,1,max(DP_A0,a-1)))/2);

      DP.SetAgeLastBDay(t,DP_MaternalNonAIDS_ALB,a-1, (DP.GetOrphansTotal(t,DP_Orphans_Maternal,DP_Orphans_NON_AIDS,0,a)
                                                   + DP.GetOrphansTotal(t,DP_Orphans_Maternal,DP_Orphans_NON_AIDS,1,a)
                                                   +DP.GetOrphansTotal(t,DP_Orphans_Maternal,DP_Orphans_NON_AIDS,0,max(DP_A0,a-1))
                                                  + DP.GetOrphansTotal(t,DP_Orphans_Maternal,DP_Orphans_NON_AIDS,1,max(DP_A0,a-1)))/2);

      DP.SetAgeLastBDay(t,DP_PaternalAIDS_ALB,a-1, (DP.GetOrphansTotal(t,DP_Orphans_Paternal,DP_Orphans_AIDS,0,a)
                                                 + DP.GetOrphansTotal(t,DP_Orphans_Paternal,DP_Orphans_AIDS,1,a)
                                                 +DP.GetOrphansTotal(t,DP_Orphans_Paternal,DP_Orphans_AIDS,0,max(DP_A0,a-1))
                                                 + DP.GetOrphansTotal(t,DP_Orphans_Paternal,DP_Orphans_AIDS,1,max(DP_A0,a-1)))/2);

      DP.SetAgeLastBDay(t,DP_PaternalNonAIDS_ALB,a-1, (DP.GetOrphansTotal(t,DP_Orphans_Paternal,DP_Orphans_NON_AIDS,0,a)
                                                   + DP.GetOrphansTotal(t,DP_Orphans_Paternal,DP_Orphans_NON_AIDS,1,a)
                                                   +DP.GetOrphansTotal(t,DP_Orphans_Paternal,DP_Orphans_NON_AIDS,0,max(DP_A0,a-1))
                                                   + DP.GetOrphansTotal(t,DP_Orphans_Paternal,DP_Orphans_NON_AIDS,1,max(DP_A0,a-1)))/2);
    end;
  end;

  for t := GB.GetCalcYearIdx(p) to DP.GetFinalIndex do
  begin
    for a := DP_A0 to DP_MaxChildAge do
    begin
      DP.SetOrphans(DP_Maternal_AIDS,a-1, t, (DP.GetOrphansTotal(t,DP_Orphans_Maternal,DP_Orphans_AIDS,0,a)
                                                +DP.GetOrphansTotal(t,DP_Orphans_Maternal,DP_Orphans_AIDS,1,a)));

      DP.SetOrphans(DP_Maternal_NON_AIDS,a-1, t, (DP.GetOrphansTotal(t,DP_Orphans_Maternal,DP_Orphans_NON_AIDS,0,a)
                                                   + DP.GetOrphansTotal(t,DP_Orphans_Maternal,DP_Orphans_NON_AIDS,1,a)));

      DP.SetOrphans(DP_Paternal_AIDS ,a-1, t, (DP.GetOrphansTotal(t,DP_Orphans_Paternal,DP_Orphans_AIDS,0,a)
                                                 + DP.GetOrphansTotal(t,DP_Orphans_Paternal,DP_Orphans_AIDS,1,a)));

      DP.SetOrphans(DP_Paternal_NON_AIDS,a-1, t, (DP.GetOrphansTotal(t,DP_Orphans_Paternal,DP_Orphans_NON_AIDS,0,a)
                                                   + DP.GetOrphansTotal(t,DP_Orphans_Paternal,DP_Orphans_NON_AIDS,1,a)));

      DP.SetOrphans(DP_Maternal_Orphans,a-1,t,  DP.GetOrphansTotal(t,DP_Orphans_Maternal,DP_Orphans_Total,0,a)
                                             +DP.GetOrphansTotal(t,DP_Orphans_Maternal,DP_Orphans_Total,1,a));

      DP.SetOrphans(DP_Paternal_Orphans,a-1,t,  DP.GetOrphansTotal(t,DP_Orphans_Paternal,DP_Orphans_Total,0,a)
                                             +DP.GetOrphansTotal(t,DP_Orphans_Paternal,DP_Orphans_Total,1,a));
    end;
  end;

  for t := GB.GetCalcYearIdx(p) to DP.GetFinalIndex do
    for a := DP_A0 to DP_MaxChildAge do
    begin
      sum := 0;
      for c := DP_Maternal_AIDS to DP_Paternal_Orphans  do
      begin
        sum := sum + DP.GetOrphans(c,a-1,t);
        //DP.SetOrphans(DP_Total_Orphans,a-1,t, DP.GetOrphans(DP_Total_Orphans,a-1,t) + DP.GetOrphans(c,a-1,t));
      end;
      DP.SetOrphans(DP_Total_Orphans,a-1,t, sum);
    end;

  for t := GB.GetCalcYearIdx(p) to DP.GetFinalIndex do
  begin
    DP.SetTotalOrphans(DP_MaternalTotal,DP_AIDS,t, 0);
    DP.SetTotalOrphans(DP_PaternalTotal,DP_AIDS,t, 0) ;
    DP.SetTotalOrphans(DP_PaternalTotal,DP_Non_AIDS,t, 0);
    DP.SetTotalOrphans(DP_MaternalTotal,DP_Non_AIDS,t, 0);
    DP.SetTotalOrphans(DP_MaternalTotal,DP_Total_AIDS_NonAIDS,t,0);
    DP.SetTotalOrphans(DP_PaternalTotal,DP_Total_AIDS_NonAIDS,t,0);

    for a := DP_A0-1 to DP_MaxChildAge-1 do
    begin
      DP.SetTotalOrphans(DP_MaternalTotal,DP_AIDS,t, DP.GetTotalOrphans(DP_MaternalTotal,DP_AIDS,t)
                                                +  DP.GetAgeLastBDay(t,DP_MaternalAIDS_ALB,a))  ;
      DP.SetTotalOrphans(DP_PaternalTotal,DP_AIDS,t, DP.GetTotalOrphans(DP_PaternalTotal,DP_AIDS,t)
                                                +  DP.GetAgeLastBDay(t,DP_PaternalAIDS_ALB,a))  ;

      DP.SetTotalOrphans(DP_MaternalTotal,DP_Non_AIDS,t, DP.GetTotalOrphans(DP_MaternalTotal,DP_Non_AIDS,t)
                                                +  DP.GetAgeLastBDay(t,DP_MaternalNonAIDS_ALB,a))  ;
      DP.SetTotalOrphans(DP_PaternalTotal,DP_Non_AIDS,t, DP.GetTotalOrphans(DP_PaternalTotal,DP_Non_AIDS,t)
                                                +  DP.GetAgeLastBDay(t,DP_PaternalNonAIDS_ALB,a))  ;

      DP.SetTotalOrphans(DP_MaternalTotal,DP_Total_AIDS_NonAIDS,t, DP.GetTotalOrphans(DP_MaternalTotal,DP_Total_AIDS_NonAIDS,t)
                                                +  DP.GetAgeLastBDay(t,DP_MaternalNonAIDS_ALB,a)
                                                +  DP.GetAgeLastBDay(t,DP_MaternalAIDS_ALB,a))  ;
      DP.SetTotalOrphans(DP_PaternalTotal,DP_Total_AIDS_NonAIDS,t, DP.GetTotalOrphans(DP_PaternalTotal,DP_Total_AIDS_NonAIDS,t)
                                                +  DP.GetAgeLastBDay(t,DP_PaternalNonAIDS_ALB,a)
                                                +  DP.GetAgeLastBDay(t,DP_PaternalAIDS_ALB,a))  ;
    end;
  end;
 (*.............END TOTAL ORPHANS.............*)

 (*.............START NEW ORPHANS.............*)
  {for t := GB.GetCalcYearIdx(p) to DP.GetFinalIndex do
  begin
    for a := DP_A0 to DP_MaxChildAge do
    begin
      DP.SetNewOrphAgeLastBDay(t,DP_MaternalAIDS_ALB,a-1, (DP.GetOrphansNew(t,DP_Orphans_Maternal,DP_Orphans_AIDS,0,a)
                                                        +DP.GetOrphansNew(t,DP_Orphans_Maternal,DP_Orphans_AIDS,1,a)
                                                        +DP.GetOrphansNew(t,DP_Orphans_Maternal,DP_Orphans_AIDS,0,max(DP_A0,a-1))
                                                        +DP.GetOrphansNew(t,DP_Orphans_Maternal,DP_Orphans_AIDS,1,max(DP_A0,a-1)))/2);

      DP.SetNewOrphAgeLastBDay(t,DP_MaternalNonAIDS_ALB,a-1,( DP.GetOrphansNew(t,DP_Orphans_Maternal,DP_Orphans_NON_AIDS,0,a)
                                                         + DP.GetOrphansNew(t,DP_Orphans_Maternal,DP_Orphans_NON_AIDS,1,a)
                                                         +DP.GetOrphansNew(t,DP_Orphans_Maternal,DP_Orphans_NON_AIDS,0,max(DP_A0,a-1))
                                                         + DP.GetOrphansNew(t,DP_Orphans_Maternal,DP_Orphans_NON_AIDS,1,max(DP_A0,a-1)))/2);

      DP.SetNewOrphAgeLastBDay(t,DP_PaternalAIDS_ALB,a-1, (DP.GetOrphansNew(t,DP_Orphans_Paternal,DP_Orphans_AIDS,0,a)
                                                         + DP.GetOrphansNew(t,DP_Orphans_Paternal,DP_Orphans_AIDS,1,a)
                                                         +DP.GetOrphansNew(t,DP_Orphans_Paternal,DP_Orphans_AIDS,0,max(DP_A0,a-1))
                                                         + DP.GetOrphansNew(t,DP_Orphans_Paternal,DP_Orphans_AIDS,1,max(DP_A0,a-1)))/2 );

      DP.SetNewOrphAgeLastBDay(t,DP_PaternalNonAIDS_ALB,a-1,( DP.GetOrphansNew(t,DP_Orphans_Paternal,DP_Orphans_NON_AIDS,0,a)
                                                         + DP.GetOrphansNew(t,DP_Orphans_Paternal,DP_Orphans_NON_AIDS,1,a)
                                                         +DP.GetOrphansNew(t,DP_Orphans_Paternal,DP_Orphans_NON_AIDS,0,max(DP_A0,a-1))
                                                         + DP.GetOrphansNew(t,DP_Orphans_Paternal,DP_Orphans_NON_AIDS,1,max(DP_A0,a-1)))/2 );
    end;
  end;
  }

  for t := GB.GetCalcYearIdx(p) to DP.GetFinalIndex do
  begin
    for a := DP_A0 to DP_MaxChildAge do
    begin
      DP.SetNewOrphans(DP_Maternal_AIDS,a-1, t, DP.GetOrphansNew(t,DP_Orphans_Maternal,DP_Orphans_AIDS,0,a)
                                                        +DP.GetOrphansNew(t,DP_Orphans_Maternal,DP_Orphans_AIDS,1,a));
      DP.SetNewOrphans(DP_Maternal_NON_AIDS,a-1, t, DP.GetOrphansNew(t,DP_Orphans_Maternal,DP_Orphans_NON_AIDS,0,a)
                                                         + DP.GetOrphansNew(t,DP_Orphans_Maternal,DP_Orphans_NON_AIDS,1,a) );
      DP.SetNewOrphans(DP_Paternal_AIDS,a-1, t, DP.GetOrphansNew(t,DP_Orphans_Paternal,DP_Orphans_AIDS,0,a)
                                                         + DP.GetOrphansNew(t,DP_Orphans_Paternal,DP_Orphans_AIDS,1,a) );
      DP.SetNewOrphans(DP_Paternal_NON_AIDS,a-1, t, DP.GetOrphansNew(t,DP_Orphans_Paternal,DP_Orphans_NON_AIDS,0,a)
                                                         + DP.GetOrphansNew(t,DP_Orphans_Paternal,DP_Orphans_NON_AIDS,1,a) );
      DP.SetNewOrphans(DP_Maternal_Orphans,a-1,t,  DP.GetOrphansNew(t,DP_Orphans_Maternal,DP_Orphans_Total,0,a)
                                                        -DP.GetOrphansNew(t,DP_Orphans_Maternal,DP_Orphans_Total,1,a));
      DP.SetNewOrphans(DP_Paternal_Orphans,a-1,t,  DP.GetOrphansNew(t,DP_Orphans_Paternal,DP_Orphans_Total,0,a)
                                                        -DP.GetOrphansNew(t,DP_Orphans_Paternal,DP_Orphans_Total,1,a));
    end;
  end;

  for t := GB.GetCalcYearIdx(p) to DP.GetFinalIndex do
    for a := DP_A0 to DP_MaxChildAge do
      for c :=  DP_Maternal_AIDS to  DP_Paternal_Orphans  do
         DP.SetNewOrphans(DP_Total_Orphans,a-1,t, DP.GetNewOrphans(DP_Total_Orphans,a-1,t) + DP.GetNewOrphans(c,a-1,t));

  {for t := GB.GetCalcYearIdx(p) to DP.GetFinalIndex do
    for s := DP_Total_AIDS_NonAIDS to DP_Non_AIDS do
      for c := DP_MaternalTotal to DP_PaternalTotal do
        DP.SetTotalNewOrphans(c,s,t, DP.GetOrphansNew(t, c, s, 0, DP_AllSingAges ) + DP.GetOrphansNew(t, c, s, 1, DP_AllSingAges ));
  }
 (*.............END NEW ORPHANS.............*)

 DualOrphanCalcs();
 TotalOrphanCalcs();

  for s:= DP_Orphans_Maternal to DP_Orphans_Paternal do
  begin
    for t := GB.GetCalcYearIdx(p) to DP.GetFinalIndex do
    begin
      tempVal := 0;
      for a := DP_A0 to DP_MaxChildAge do
        tempVal:= tempVal + DP.GetBigPop(t, DP_BothSexes, a);
    end;
  end;
end;

function TDPProj.PaternalOrphanCalcs : boolean;
var
  orphsNoHIV, orphsHIVNoART, orphsHIVOnART, temp,
  a, dVal, cumSurvR, tempSurvR, tempVal,
  pNoHIVcNoHIV, pNoHIVcHIV, pHIVNoARTcNoHIV,pHIVNoARTcHIV,
  pHIVOnARTcNoHIV, pHIVOnARTcHIV,alphaVal,betaVal, fVal,
  PInfMInf, Ks, HIVMfert, probMTCT, childSurv, tempfertNum,
  tempfertDen, tempfert, prevNum, prevDen                      : Double;
  AART, AatInf, orpA , parentA, dA, dCD4, a5, aPast, i, j      : Byte;
  tBorn, t, tempA,d, tSurv, tProxy                             : Integer;
  f                                                            : Textfile;
  MaleASFR                                                     : Array[ DP_A15_19..DP_A80_Up] of Double;
  MaleTFR                                                      : Array of Double;
  FemalePrev                                                   : array of single;
  MalePrev                                                     : array of single;
begin
  SetLength(MaleTFR,GB.GetNumDynYears(p));
  SetLength(FemalePrev,GB.GetNumDynYears(p));
  SetLength(MalePrev,GB.GetNumDynYears(p));

  Result := false;

         //Constant values for Male fertility
  MaleASFR	[DP_A15_19]	:=	0.00111221013656935	*	5	*	100;
  MaleASFR	[DP_A20_24]	:=	0.0181909386334066	*	5	*	100;
  MaleASFR	[DP_A25_29]	:=	0.0357020127056169	*	5	*	100;
  MaleASFR	[DP_A30_34]	:=	0.0369100610966247	*	5	*	100;
  MaleASFR	[DP_A35_39]	:=	0.0318622561819123	*	5	*	100;
  MaleASFR	[DP_A40_44]	:=	0.0254888371283555	*	5	*	100;
  MaleASFR	[DP_A45_49]	:=	0.0192506224379379	*	5	*	100;
  MaleASFR	[DP_A50_54]	:=	0.0138770855271757	*	5	*	100;
  MaleASFR	[DP_A55_59]	:=	0.00938047697522717	*	5	*	100;
  MaleASFR	[DP_A60_64]	:=	0.0054987622371946	*	5	*	100;
  MaleASFR	[DP_A65_69]	:=	0.00219993086976762	*	5	*	100;
  MaleASFR	[DP_A70_74]	:=	0.000475268889966518	*	5	*	100;
  MaleASFR	[DP_A75_79]	:=	5.14947E-05	*	5	*	100;
  MaleASFR	[DP_A80_Up]	:=	0;

      //Estimate Male TFR
  for t := GB.GetCalcYearIdx(p) to DP.GetFinalIndex do
  begin
    temp := 0 ;
    for a5 :=  DP_A15_19 to DP_A75_79 do
    begin
      temp :=  temp + DP.GetPop(t, DP_Male, a5) * MaleASFR [a5] / 5/ 100;
    end;
   if temp <= 0
     then MaleTFR[t] := 0
     else MaleTFR [t] :=  DP.GetBirths(t, DP_AllAges,DP_BothSexes) / (temp);
  end;

   (* Concordance Matrix *)
   (* Female HIV Prevalence *)
  for t := GB.GetCalcYearIdx(p)  to DP.GetFinalIndex do
  begin
    prevNum := 0;
    prevDen := 0;
    for a5 := DP_A15_19 to DP_A45_49  do
    begin
      prevNum := prevNum + DP.GetHIV(t, DP_Female, a5) ;
      prevDen := prevDen + DP.GetPop(t, DP_Female, a5);
    end;
    if prevDen <= 0 then FemalePrev[t] := 0
                    else FemalePrev[t] :=  prevNum / prevDen;

    prevNum := 0;
    prevDen := 0;
    for a5 := DP_A15_19 to DP_A75_79  do
    begin
      prevNum := prevNum + DP.GetHIV(t, DP_Male, a5) ;
      prevDen := prevDen + DP.GetPop(t, DP_Male, a5);
    end;
    if prevDen <= 0 then MalePrev[t] := 0
                    else MalePrev[t] :=  prevNum / prevDen;
  end;

  alphaVal := -0.838694494324521;  //regression coefficient for estimating pr
  betaVal := 9.13404971128;
  Ks := 0.64; //relative reduction in fertility in women from HIV

  (*Calculate duration of infection for time during no ART by CD4 count*)
  CreateInfDur (1, DP_Male)      ;    // transition matrix does not change over time so setting t = 1
  (*Calculate duration of infection on ART by CD4 count*)
  setDurARTIni(DP_Male) ; //

  (* Calculate Paternal AIDS, Non-AIDS, and Total Orphans   *)
  for t := GB.GetCalcYearIdx(p)-DP_MaxChildAge to DP.GetFinalIndex do
  begin
    // tProxy is used to access quantities like ASFR that are unavailable before
    // the projection begins
    tProxy := max(GB.GetCalcYearIdx(p), t);
    for orpA := DP_A0 to DP_MaxChildAge do
    begin
      tBorn := t - orpA + 1;

      cumSurvR := 1;
      //Estimate survival of child to age orpA
      for aPast := DP_A0 to orpA do
      begin
        // Use survival rates for the projection start year for years of the
        // child's life that precede the projection start year
        tSurv := max(GB.GetCalcYearIdx(p), tBorn - 1 + aPast);
        cumSurvR := cumSurvR * (Get_Survival_Rate (DP, tSurv, DP_Male,   aPast-1) +
                                Get_Survival_Rate (DP, tSurv, DP_Female, aPast-1))/2;
      end;
      // clamp child birth year to the projection interval
      if tBorn < GB.GetCalcYearIdx(p) then tBorn := GB.GetCalcYearIdx(p);

      pNoHIVcNoHIV := 0;
      pNoHIVcHIV := 0;
      pHIVNoARTcNoHIV := 0;
      pHIVNoARTcHIV := 0;
      pHIVOnARTcNoHIV := 0;
      pHIVOnARTcHIV := 0;
      for parentA := (DP_A15 + orpA ) to (DP_A79) do
      begin
        dA := trunc((parentA - orpA) / 5) + 1;   //age of father at birth of child
        if MaleASFR[dA] / 5 > 0 then
        begin
          //Determine probabilities of children getting infected at birth or through BF
          probMTCT:= 0;
          childSurv := 0;
          //set probMTCT and childSurv values
          if DP.GetNeedPMTCT(tBorn) > 0 then
          begin
            case orpA of
              DP_A0..DP_A14           : probMTCT :=  DP.GetPerinatalTransmission(min(tBorn + 2,DP.GetFinalIndex))/100;
              DP_A15..DP_MaxChildAge  : probMTCT :=  DP.GetPerinatalTransmission(min(tBorn + 2,DP.GetFinalIndex))/100;
            end;
            if probMTCT > 0 then
            begin
              if orpA > DP_A14 then
              //if statements needed for countries with low HIV
              begin
                if DP.GetNeedPMTCT(tBorn + orpA - DP_A14) > 0 then
                  childSurv :=  (DP.GetHIVChildSingleAge(tProxy,DP_A14) / (DP.GetNeedPMTCT(tBorn + orpA - DP_A14) * probMTCT))
                else
                  childSurv := cumSurvR;
                for i := DP_A15 to orpA  do
                  if  DP.GetHIVChildSingleAge(max(t-1,1),DP_A13) > 0   then
                    childSurv := childSurv * DP.GetHIVChildSingleAge(tProxy,DP_A14)/ DP.GetHIVChildSingleAge(max(t-1,1),DP_A13);
              end
              else
              begin
                childSurv := DP.GetHIVChildSingleAge(tProxy,orpA) / (DP.GetNeedPMTCT(tBorn) * probMTCT);
              end;
            end;
          end; //  end of set probMTCT and childSurv values

          childSurv := min (childSurv, cumSurvR);
          //probability mother infected given father infected
          fVal := Exp(alphaVal + betaVal * FemalePrev[tBorn])/(1 + Exp(alphaVal + betaVal * FemalePrev[tBorn]));
          //probability mother infected given father uninfected
          PInfMinf := 0;
          if (MalePrev[tBorn] > 0) and (MalePrev[tBorn] < 1) then // RLG 2017-08-21: added <1 check to avoid divide-by-zero during fitting
            PInfMinf :=  MalePrev[tBorn]*((FemalePrev[tBorn]/MalePrev[tBorn]) - fVal)/(1- MalePrev[tBorn]);

          //Fertility factor when mother is infected based on ART or noART
          //HIVMfert := MaleTFR[tBorn]* (MaleASFR[dA] / 5 / 100);
          tempfert := 1;
          tempfertNum := 0;
          tempfertDen := 0;
          if DP.GetHIV(tBorn,DP_Female, DP_AllAges) > 0 then
          begin
            for i := DP_A15_19 to DP_A45_49 do
            begin
               tempfertNum := tempfertNum +
                              ((DP.GetHIV(tBorn,DP_Female, i) - DP.GetOnART(tBorn,DP_Female,i))* (DP.GetASFR(tBorn, i) /100) *
                                DP.GetHIVTFR(DP_Data, i, tBorn) + (DP.GetASFR(tBorn, i) /100 )* DP.GetOnART(tBorn,DP_Female,i));
               tempfertDen := tempfertDen + DP.GetHIV(tBorn,DP_Female, i) *  (DP.GetASFR(tBorn, i) /100 );
            end;
            if tempfertDen > 0 then
             tempfert := tempfertNum / tempfertDen
            else
              tempfert := 0;
          end;
          HIVMfert := MaleTFR[tBorn]* (MaleASFR[dA] / 5 / 100)  * tempfert ;

        //1.Orphans from Non-AIDS  fathers
         //1a. mother uninfected  given father uninfected
          if parentA = DP_A80 then
            dVal := DP.GetDeaths(tProxy, DP_A80_Up, DP_Male)* MaleTFR[tBorn]* (MaleASFR[dA] / 5 / 100) *(1- PInfMinf)
          else
            dVal := DP.GetDeathsBySingleAge(tProxy, parentA, DP_Male)* MaleTFR[tBorn]* (MaleASFR[dA] / 5 / 100) *(1- PInfMinf);
          pNoHIVcNoHIV := pNoHIVcNoHIV + dVal *  cumSurvR;

        //1b. mother infected  given father uninfected: reduced fertility if mother not on PMTCT
          if parentA = DP_A80 then
            dVal := DP.GetDeaths(tProxy, DP_A80_Up, DP_Male)* PInfMinf * HIVMfert
          else
            dVal := DP.GetDeathsBySingleAge(tProxy, parentA, DP_Male)* PInfMinf * HIVMfert;
          pNoHIVcNoHIV := pNoHIVcNoHIV + dVal *  (1- probMTCT) * cumSurvR;
          pNoHIVcHIV := pNoHIVcHIV + dVal* probMTCT * childSurv;

        //2. Orphans from HIV fathers  not on ART
          for dCD4 := DP_CD4_GT500 to DP_CD4_LT50  do
          begin
            AatInf := getAgeAtInfection (DP_Male, parentA-1, dCD4);

            //2a.Births before infection
            if (parentA -orpA) < (AatInf - 1) then
            begin
              //2.a.i. mother uninfected  given father uninfected at time of birth
              dVal := DP.GetAIDSDeathsCD4NoART (tProxy, parentA, DP_Male, dCD4, DP_D_HIVPos)* MaleTFR[tBorn]* (MaleASFR[dA] / 5 / 100) *(1- PInfMinf);
              pHIVNoARTcNoHIV := pHIVNoARTcNoHIV + dVal *  cumSurvR;

              //2.a.ii. mother infected  given father uninfected at time of birth
              dVal := DP.GetAIDSDeathsCD4NoART (tProxy, parentA, DP_Male, dCD4, DP_D_HIVPos) * PInfMinf  * HIVMfert;
              pHIVNoARTcNoHIV := pHIVNoARTcNoHIV + dVal *  (1- probMTCT) * cumSurvR;
              pHIVNoARTcHIV := pHIVNoARTcHIV + dVal* probMTCT * childSurv;
            end
           //2.b. Births during infection
            else
            begin
              //2.b.i. mother uninfected  given father infected at time of birth
              dVal := DP.GetAIDSDeathsCD4NoART (tProxy, parentA, DP_Male, dCD4, DP_D_HIVPos)* MaleTFR[tBorn]* (MaleASFR[dA] / 5 / 100) *(1- fVal);
              pHIVNoARTcNoHIV := pHIVNoARTcNoHIV + dVal *  cumSurvR;

              //2.b.ii. mother infected  given father infected at time of birth
              dVal := DP.GetAIDSDeathsCD4NoART (tProxy, parentA, DP_Male, dCD4, DP_D_HIVPos)* fVal * HIVMfert;
              pHIVNoARTcNoHIV := pHIVNoARTcNoHIV + dVal *  (1- probMTCT) * cumSurvR;
              pHIVNoARTcHIV := pHIVNoARTcHIV + dVal* probMTCT * childSurv;
            end;
          end;
        //End of Orphans from HIV  fathers  NOT on ART

        //3. Orphans from HIV  fathers  on ART
          for dCD4 := DP_CD4_GT500 to DP_CD4_LT50  do
          begin
            for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
            begin
              AART :=  parentA;
              case d of
                DP_D_ARTlt6m   :   AART := parentA;
                DP_D_ART6to12m :   AART := parentA - 1;
                DP_D_ARTgt12m  :   if DP.GetAIDSDeathsCD4OnART(tProxy, parentA, DP_Male, dCD4, DP_D_ARTgt12m) > 0 then
                                AART := max(15, parentA -  max(0,round((durARTMatrix[DP_male, dCD4, tProxy] *
                                           DP.GetAIDSDeathsCD4OnART(tProxy, parentA, DP_Male, dCD4, DP_D_All) - 0.5 *
                                           DP.GetAIDSDeathsCD4OnART(tProxy, parentA, DP_Male, dCD4, DP_D_ARTlt6m) - 1 *
                                           DP.GetAIDSDeathsCD4OnART(tProxy, parentA, DP_Male, dCD4, DP_D_ART6to12m)) /
                                           DP.GetAIDSDeathsCD4OnART(tProxy, parentA, DP_Male, dCD4, DP_D_ARTgt12m))));
              end;
              AatInf := max(15,getAgeAtInfection (DP_male, AART-1 , dCD4));

              (*3a. Births  before infection  *)
              if (parentA -orpA ) < (AatInf - 1) then
              begin
                (*3.a.i. mother uninfected  given father uninfected at time of birth  *)
                dVal := DP.GetAIDSDeathsCD4OnART(tProxy, parentA, DP_Male, dCD4, d)* MaleTFR[tBorn]* (MaleASFR[dA] / 5 / 100) *(1- PInfMinf);
                pHIVOnARTcNoHIV := pHIVOnARTcNoHIV + dVal *  cumSurvR;

                (*3.a.ii. mother infected  given father uninfected at time of birth*)
                dVal := DP.GetAIDSDeathsCD4OnART(tProxy, parentA, DP_Male, dCD4, d)* PInfMinf *HIVMfert;
                pHIVOnARTcNoHIV := pHIVOnARTcNoHIV + dVal *  (1- probMTCT) * cumSurvR;
                pHIVOnARTcHIV := pHIVOnARTcHIV + dVal* probMTCT * childSurv;
              end
            (*3.b. Births after infection *)
              else
              begin
                //3.b.i. mother uninfected  given father infected at time of birth
                dVal := DP.GetAIDSDeathsCD4OnART(tProxy, parentA, DP_Male, dCD4, d)* MaleTFR[tBorn]* (MaleASFR[dA] / 5 / 100) *(1- fVal);
                pHIVOnARTcNoHIV := pHIVOnARTcNoHIV + dVal *  cumSurvR;

                //3.b.ii. mother infected  given father infected at time of birth
                dVal := DP.GetAIDSDeathsCD4OnART(tProxy, parentA, DP_Male, dCD4, d)  * fVal *HIVMfert;
                pHIVOnARTcNoHIV := pHIVOnARTcNoHIV + dVal *  (1- probMTCT) * cumSurvR;
                pHIVOnARTcHIV := pHIVOnARTcHIV + dVal* probMTCT * childSurv;
              end;//end of IF before infection
            end; // end d
          end; //end of CD4 FOR
           //End of Orphans from HIV  father  on ART
        end; //end of IF to determine if parentA within fertiliy age
      end;  // end of parentA

      if orpA -1 = 0 then
      begin
        (* Dividing by 2 assuming parent dies mid year and so half of full year births *)
        pNoHIVcNoHIV := pNoHIVcNoHIV / 2;
        pNoHIVcHIV   := pNoHIVcHIV / 2;
        pHIVNoARTcNoHIV := pHIVNoARTcNoHIV / 2;
        pHIVNoARTcHIV :=  pHIVNoARTcHIV / 2;
        pHIVOnARTcNoHIV := pHIVOnARTcNoHIV / 2;
        pHIVOnARTcHIV := pHIVOnARTcHIV / 2;
      end;
      //new orphans of each age
      DP.SetOrphansNew(t, DP_Orphans_Paternal, DP_Orphans_NON_AIDS, 0, orpA, pNoHIVcNoHIV);
      DP.SetOrphansNew(t, DP_Orphans_Paternal, DP_Orphans_NON_AIDS, 1, orpA, pNoHIVcHIV);
      DP.SetOrphansNew(t, DP_Orphans_Paternal, DP_Orphans_AIDS, 0, orpA , pHIVNoARTcNoHIV + pHIVOnARTcNoHIV);
      DP.SetOrphansNew(t, DP_Orphans_Paternal, DP_Orphans_AIDS, 1, orpA , pHIVNoARTcHIV + pHIVOnARTcHIV);
      DP.SetOrphansNew(t, DP_Orphans_Paternal, DP_Orphans_Total, 0, orpA , pNoHIVcNoHIV+ pHIVNoARTcNoHIV + pHIVOnARTcNoHIV);
      DP.SetOrphansNew(t, DP_Orphans_Paternal, DP_Orphans_Total, 1, orpA , pNoHIVcHIV + pHIVNoARTcHIV + pHIVOnARTcHIV);

      //New orphans summed over all ages
      for i := DP_Orphans_Total to DP_Orphans_NON_AIDS  do
      begin
        for j := 0 to 1 do
          begin
            DP.SetOrphansNew(t, DP_Orphans_Paternal, i, j, DP_AllSingAges,
                      DP.GetOrphansNew(t, DP_Orphans_Paternal, i, j, DP_AllSingAges)
                      + DP.GetOrphansNew(t, DP_Orphans_Paternal, i, j, orpA));
        end;
      end;

      //if orphans age is zero then total orphans of age 0 = new orphans of age zero
      if orpA -1 = 0 then
      begin
        DP.SetOrphansTotal(t, DP_Orphans_Paternal, DP_Orphans_NON_AIDS, 0, orpA, pNoHIVcNoHIV);
        DP.SetOrphansTotal(t, DP_Orphans_Paternal, DP_Orphans_NON_AIDS, 1, orpA, pNoHIVcHIV);
        DP.SetOrphansTotal(t, DP_Orphans_Paternal, DP_Orphans_AIDS, 0, orpA , pHIVNoARTcNoHIV + pHIVOnARTcNoHIV);
        DP.SetOrphansTotal(t, DP_Orphans_Paternal, DP_Orphans_AIDS, 1, orpA , pHIVNoARTcHIV + pHIVOnARTcHIV);
        DP.SetOrphansTotal(t, DP_Orphans_Paternal, DP_Orphans_Total, 0, orpA , pNoHIVcNoHIV+ pHIVNoARTcNoHIV + pHIVOnARTcNoHIV);
        DP.SetOrphansTotal(t, DP_Orphans_Paternal, DP_Orphans_Total, 1, orpA , pNoHIVcHIV + pHIVNoARTcHIV + pHIVOnARTcHIV);
      end
      else
      begin
        (*Total orphans = NewOrphans + Orphans in age orpA-1 in Total Orphans * Survival to age orpA*)

       //Total uninfected orphans from uninfected fathers
        DP.SetOrphansTotal(t, DP_Orphans_Paternal, DP_Orphans_NON_AIDS, 0, orpA,
                                    pNoHIVcNoHIV + DP.GetOrphansTotal(t-1, DP_Orphans_Paternal, DP_Orphans_NON_AIDS, 0, orpA - 1) *
                                    (Get_Survival_Rate(DP, tProxy, DP_Male, orpA -1)+Get_Survival_Rate (DP, tProxy, DP_FeMale, orpA -1))/2);

        //Total uninfected orphans from infected parent
         DP.SetOrphansTotal(t, DP_Orphans_Paternal, DP_Orphans_AIDS, 0, orpA,
                                pHIVNoARTcNoHIV + pHIVOnARTcNoHIV +
                                DP.GetOrphansTotal(t-1, DP_Orphans_Paternal, DP_Orphans_AIDS, 0, orpA - 1) *
                                (Get_Survival_Rate (DP, tProxy, DP_Male, orpA-1 )+Get_Survival_Rate (DP, tProxy, DP_FeMale, orpA -1))/2);


        tempVal := (Get_Survival_Rate (DP, tProxy, DP_Male, orpA -1)+Get_Survival_Rate (DP, tProxy, DP_FeMale, orpA -1))/2;
        tempA := max(orpA, 14);

        if tempA <= DP_A2 then
        begin
          if (DP.GetNeedPMTCT(max(t - orpA+1,1)) * DP.GetPerinatalTransmission(tProxy + 2)/100)> 0 then
            tempVal := DP.GetHIVChildSingleAge(tProxy,orpA) / (DP.GetNeedPMTCT(max(t - orpA+1,1)) * DP.GetPerinatalTransmission(tProxy + 2)/100) ;
        end
        else
        begin
           if DP.GetHIVChildSingleAge(max(t-1,1),tempA-1) > 0 then
            tempVal := DP.GetHIVChildSingleAge(tProxy,tempA)/DP.GetHIVChildSingleAge(max(t-1,1),tempA-1);
        end;

        tempVal := min(tempVal, (Get_Survival_Rate(DP, tProxy, DP_Male, orpA -1)+Get_Survival_Rate(DP, tProxy, DP_FeMale, orpA -1))/2);

        //Total infected orphans from uninfected fathers
        DP.SetOrphansTotal(t, DP_Orphans_Paternal, DP_Orphans_NON_AIDS, 1, orpA,  +  pNoHIVcHIV
        + (DP.GetOrphansTotal(t-1, DP_Orphans_Paternal, DP_Orphans_NON_AIDS, 1, orpA - 1) *  tempVal));


        //Total infected orphans from infected fathers
        DP.SetOrphansTotal(t, DP_Orphans_Paternal, DP_Orphans_AIDS, 1, orpA,  +  pHIVNoARTcHIV + pHIVOnARTcHIV
        + (DP.GetOrphansTotal(t-1, DP_Orphans_Paternal, DP_Orphans_AIDS, 1, orpA - 1) *  tempVal));

        for i := 0 to 1 do
        DP.SetOrphansTotal(t, DP_Orphans_Paternal, DP_Orphans_Total, i, orpA,
                          DP.GetOrphansTotal(t, DP_Orphans_Paternal, DP_Orphans_NON_AIDS, i, orpA)
                          + DP.GetOrphansTotal(t, DP_Orphans_Paternal, DP_Orphans_AIDS, i, orpA));
     end;

     for i := DP_Orphans_Total to DP_Orphans_NON_AIDS  do
       for j := 0 to 1 do
         DP.SetOrphansTotal(t, DP_Orphans_Paternal, i,j, DP_AllSingAges,
                          DP.GetOrphansTotal(t, DP_Orphans_Paternal, i,j, DP_AllSingAges) +
                          DP.GetOrphansTotal(t, DP_Orphans_Paternal, i,j, orpA));

    end;    //end of orpA
  end; //end of t loop
  Result := true;

  Finalize(MaleTFR);
  Finalize(FemalePrev);
  Finalize(MalePrev);
end;

(*Start of Maternal orphans*)
function TDPProj.MaternalOrphanCalcs : boolean;
var
  orphsNoHIV, orphsHIVNoART, orphsHIVOnART, temp,
  a, dVal, cumSurvR, MHR1, MHR2, tempSurvR, tempVal,
  pNoHIVcNoHIV, pHIVNoARTcNoHIV,pHIVNoARTcHIV, births,
  pHIVOnARTcNoHIV, pHIVOnARTcHIV, probMTCT, childSurv      : Double;
  AART, AatInf, orpA, aPast, parentA, dA, dCD4, a5, i, j   : Byte;
  tBorn, t, tempA, d, tSurv, tProxy                        : Integer;
  f                                                        : Textfile;
begin
  Result := false;
  (*Calculate duration of infection for time during no ART by CD4 count*)
  CreateInfDur (1, DP_Female);    // transition matrix does not change over time so setting t = 1
  (*Calculate duration of infection on ART by CD4 count*)
  setDurARTIni(DP_Female); //
  MHR1 := 3.2; MHR2 := 6.2; //mortality hazard ratio in before and after death of mother

  (* Calculate Maternal AIDS, Non-AIDS, and Total Orphans   *)
  for t := GB.GetCalcYearIdx(p)-DP_MaxChildAge to DP.GetFinalIndex do
  begin
    // tProxy is used to access quantities like ASFR that are unavailable before
    // the projection begins
    tProxy := max(GB.GetCalcYearIdx(p), t);
    for orpA := DP_A0 to DP_MaxChildAge do
    begin
      tBorn := t - orpA + 1; // year child born

      cumSurvR := 1;
      //Estimate survival of child to age orpA by adjusting with mortality hazard ratio in year before and after mothers death
      if orpA <= DP_A1 then
      begin
        // Child survival rates in years lived before projection began are
        // assumed to match rates in the projection start year
        tSurv := max(GB.GetCalcYearIdx(p), tBorn - 1 + orpA);
        cumSurvR := (1 - MHR2 * (1 - (Get_Survival_Rate(DP, tSurv, DP_Male,   orpA-1) +
                                      Get_Survival_Rate(DP, tSurv, DP_Female, orpA-1))/2));
        if orpA = DP_A1 then
        begin
          tSurv := max(GB.GetCalcYearIdx(p), tBorn - 1 + DP_A0);
          cumSurvR := cumSurvR * (1 - MHR1 * (1 - (Get_Survival_Rate(DP, tSurv, DP_Male,   DP_A0-1) +
                                                   Get_Survival_Rate(DP, tSurv, DP_Female, DP_A0-1))/2));
        end;
      end
      else
      begin
        // Iterate through previous years of life
        for aPast := DP_A0 to orpA - 2 do
        begin
          // tSurv is the year that child age equals aPast
          tSurv := max(GB.GetCalcYearIdx(p), tBorn - 1 + aPast);
          cumSurvR := cumSurvR * (Get_Survival_Rate(DP, tSurv, DP_Male,   aPast-1) +
                                  Get_Survival_Rate(DP, tSurv, DP_Female, aPast-1))/2;
        end;
        //Year before mothers death childs age equals aPast
        aPast := orpA - 2 + 1;
        tSurv := max(GB.GetCalcYearIdx(p), tBorn - 1 + aPast);
        cumSurvR := cumSurvR *(1 - MHR1 * (1 - (Get_Survival_Rate(DP, tSurv, DP_Male,   aPast-1) +
                                                Get_Survival_Rate(DP, tSurv, DP_Female, aPast-1))/2));
        //Year at mothers death
        tSurv := max(GB.GetCalcYearIdx(p), tBorn - 1 + orpA);
        cumSurvR := cumSurvR *(1 - MHR2 * (1 - (Get_Survival_Rate(DP, tSurv, DP_Male,   orpA-1) +
                                                Get_Survival_Rate(DP, tSurv, DP_Female, orpA-1))/2));
      end;
      // clamp child birth year to the projection interval
      if tBorn < GB.GetCalcYearIdx(p) then tBorn := GB.GetCalcYearIdx(p);

      births := 0;
      pNoHIVcNoHIV := 0;
      pHIVNoARTcNoHIV := 0;
      pHIVNoARTcHIV := 0;
      pHIVOnARTcNoHIV := 0;
      pHIVOnARTcHIV := 0;
      for parentA := (DP_A15 + orpA - 1) to (DP_A79) do
      begin
        dA := trunc((parentA - orpA) / 5 ) + 1;   //age of mother at birth of child
        if (DP.GetASFR(tBorn, dA)) / 5 > 0 then
        begin
            //Determine probabilities of children getting infected at birth or through BF
          probMTCT:= 0;
          childSurv := 0;
         //set probMTCT and childSurv values
          if DP.GetNeedPMTCT(tBorn) > 0 then
          begin
            case orpA of
              DP_A0                   :   probMTCT := DP.GetMTCTRate6Wks(tProxy)/100;
              DP_A1                   :   if DP.GetNeedPMTCT(tProxy) > 0 then
                                            probMTCT := DP.GetNewInfantInfections(tProxy)/(DP.GetNeedPMTCT(tProxy));
              DP_A2..DP_A14           :   probMTCT :=  DP.GetPerinatalTransmission(min(tBorn + 2,DP.GetFinalIndex))/100;
              DP_A15..DP_MaxChildAge  :   probMTCT :=  DP.GetPerinatalTransmission(min(tBorn + 2,DP.GetFinalIndex))/100;
            end;
            if probMTCT > 0 then
            begin
              if orpA > DP_A14 then
              begin
                if DP.GetNeedPMTCT(tBorn + orpA - DP_A14) > 0 then
                  childSurv :=  (DP.GetHIVChildSingleAge(tProxy,DP_A14) / (DP.GetNeedPMTCT(tBorn + orpA - DP_A14) * probMTCT))
                else
                  childSurv := cumSurvR;
                for i := DP_A15 to orpA  do
                  if  DP.GetHIVChildSingleAge(max(t-1,1),DP_A13) > 0   then
                    childSurv := childSurv * DP.GetHIVChildSingleAge(tProxy,DP_A14)/ DP.GetHIVChildSingleAge(max(t-1,1),DP_A13);
              end
              else
              begin
                childSurv := DP.GetHIVChildSingleAge(tProxy,orpA) / (DP.GetNeedPMTCT(tBorn) * probMTCT);
              end;
            end;
          end; //  end of set probMTCT and childSurv values

           //1. Orphans from Non-AIDS  mothers
             //AIDS not included in  DP.GetDeathsBySingleAge according  to procedure in DPData
           dVal := (DP.GetDeathsBySingleAge(tProxy, parentA, DP_Female)* DP.GetTFR(tBorn)) * (DP.GetASFR(tBorn, dA)/ 5) / 100;
                     pNoHIVcNoHIV := pNoHIVcNoHIV + dVal *  cumSurvR;

             // births := births + dVal;
          //2. Orphans from HIV mothers  not on ART

           for dCD4 := DP_CD4_GT500 to DP_CD4_LT50  do
           begin
              AatInf := getAgeAtInfection (DP_Female, parentA-1, dCD4);

                    //2a. Births before infection
              if (parentA -orpA) < (AatInf - 1) then
              begin
                dVal := DP.GetAIDSDeathsCD4NoART(tProxy, parentA, DP_Female, dCD4, DP_D_HIVPos) * DP.GetTFR(tBorn) *
                        DP.GetASFR(tBorn, dA) / 5 / 100;

                // births := births + dVal;
                pHIVNoARTcNoHIV:= pHIVNoARTcNoHIV + dVal *  cumSurvR;
              end
            // 2b. Births during infection
              else
              begin
                dVal := DP.GetAIDSDeathsCD4NoART(tProxy, parentA, DP_Female, dCD4, DP_D_HIVPos)  * DP.GetTFR(tBorn) *
                        DP.GetASFR(tBorn, dA) / 5 * DP.GetHIVTFR(DP_Data, dA, tBorn) / 100;
                births := births + dVal * (probMTCT );
               // births := births + dVal * (1-probMTCT );
                pHIVNoARTcHIV := pHIVNoARTcHIV + dVal * probMTCT * childSurv;

                pHIVNoARTcNoHIV := pHIVNoARTcNoHIV + dVal * (1-probMTCT ) * cumSurvR;
              end;
           end;
           //End of Orphans from HIV  mothers  NOT on ART

           //3. Orphans from HIV  mothers  on ART
           for dCD4 := DP_CD4_GT500 to DP_CD4_LT50 do
           begin
             for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
             begin
                AART :=  parentA;
                case d of
                  DP_D_ARTlt6m   : AART :=   parentA;
                  DP_D_ART6to12m : AART :=   parentA - 1;
                  DP_D_ARTgt12m  : if DP.GetAIDSDeathsCD4OnART(tProxy, parentA, DP_Female, dCD4, DP_D_ARTgt12m) > 0 then
                                  AART := max(15, parentA - max(0, round((durARTMatrix[DP_Female, dCD4, tProxy] *
                                              DP.GetAIDSDeathsCD4OnART(tProxy, parentA, DP_Female, dCD4, DP_D_All) - 0.5 *
                                              DP.GetAIDSDeathsCD4OnART(tProxy, parentA, DP_Female, dCD4, DP_D_ARTlt6m) - 1 *
                                              DP.GetAIDSDeathsCD4OnART(tProxy, parentA, DP_Female, dCD4, DP_D_ART6to12m)) /
                                              DP.GetAIDSDeathsCD4OnART(tProxy, parentA, DP_Female, dCD4, DP_D_ARTgt12m))));
                end;

                AatInf := getAgeAtInfection (DP_Female, AART -1, dCD4);
                (*3a. Births  before infection*)
                if (parentA -orpA) < (AatInf - 1) then
                begin
                  dVal :=   DP.GetAIDSDeathsCD4OnART(tProxy, parentA, DP_FeMale, dCD4,d) * DP.GetTFR(tBorn) *
                            DP.GetASFR(tBorn, dA) / 5 / 100;
                   //births := births + dVal;

                  pHIVOnARTcNoHIV:= pHIVOnARTcNoHIV + dVal *  cumSurvR;
                end
                //3b. Births during infection but before ART
                else
                begin
                  if (parentA -orpA) >  AART then
                  begin
                    dVal := DP.GetAIDSDeathsCD4OnART(tProxy, parentA, DP_Female, dCD4,d)  * DP.GetTFR(tBorn) *
                             DP.GetASFR(tBorn, dA)  / 5 / 100;
                    births := births + dVal * (probMTCT );
                       // births := births + dVal * (1-probMTCT );
                    pHIVOnARTcHIV := pHIVOnARTcHIV + dVal * probMTCT * childSurv;
                    pHIVOnARTcNoHIV := pHIVOnARTcNoHIV + dVal * (1-  probMTCT) * cumSurvR;
                  end
                  //3c. Births during ART
                  else
                  begin
                    dVal := DP.GetAIDSDeathsCD4OnART(tProxy, parentA, DP_Female, dCD4,d)  * DP.GetTFR(tBorn) *
                            DP.GetASFR(tBorn, dA) / 5 * DP.GetHIVTFR(DP_Data, dA, tBorn)/ 100;
                    births := births + dVal * (probMTCT);
                    //  births := births + dVal * (1-probMTCT );

                    pHIVOnARTcHIV := pHIVOnARTcHIV + dVal * probMTCT * childSurv;
                    pHIVOnARTcNoHIV:= pHIVOnARTcNoHIV + dVal * (1 - probMTCT) * cumSurvR;
                  end;//end of IF
                end;
             end; //end d
           end;//end dCD4
            //End of Orphans from HIV  mothers  on ART
        end; //end of IF to determine if parentA within fertiliy age
      end;  // end of parentA

      if orpA -1 = 0 then
      begin
        (*Dividing by 2 assuming parent dies mid year and so half of full year births*)
        //births := births /2;
        pNoHIVcNoHIV := pNoHIVcNoHIV / 2;
        pHIVNoARTcNoHIV := pHIVNoARTcNoHIV / 2;
        pHIVNoARTcHIV :=  pHIVNoARTcHIV / 2;
        pHIVOnARTcNoHIV := pHIVOnARTcNoHIV / 2;
        pHIVOnARTcHIV := pHIVOnARTcHIV / 2;
      end;

      //new orphans of each age
      DP.SetOrphansNew(t, DP_Orphans_Maternal, DP_Orphans_NON_AIDS, 0, orpA, pNoHIVcNoHIV);
      DP.SetOrphansNew(t, DP_Orphans_Maternal, DP_Orphans_AIDS, 0, orpA , pHIVNoARTcNoHIV + pHIVOnARTcNoHIV);
      DP.SetOrphansNew(t, DP_Orphans_Maternal, DP_Orphans_AIDS, 1, orpA , pHIVNoARTcHIV + pHIVOnARTcHIV);
      DP.SetOrphansNew(t, DP_Orphans_Maternal, DP_Orphans_Total, 0, orpA , pNoHIVcNoHIV + pHIVNoARTcNoHIV + pHIVOnARTcNoHIV);
      DP.SetOrphansNew(t, DP_Orphans_Maternal, DP_Orphans_Total, 1, orpA ,  pHIVNoARTcHIV + pHIVOnARTcHIV);

      //New orphans summed over all ages
      for i := DP_Orphans_Total to DP_Orphans_NON_AIDS  do
        for j := 0 to 1 do
            DP.SetOrphansNew(t, DP_Orphans_Maternal, i, j, 0,
                           DP.GetOrphansNew(t, DP_Orphans_Maternal, i, j, 0) +
                           DP.GetOrphansNew(t, DP_Orphans_Maternal, i, j, orpA));


      //if orphans age is zero then total orphans of age 0 = new orphans of age zero
      if orpA -1 = 0 then
      begin
        DP.SetOrphansTotal(t, DP_Orphans_Maternal, DP_Orphans_NON_AIDS, 0, orpA, pNoHIVcNoHIV);
        DP.SetOrphansTotal(t, DP_Orphans_Maternal, DP_Orphans_AIDS, 0, orpA , pHIVNoARTcNoHIV + pHIVOnARTcNoHIV);
        DP.SetOrphansTotal(t, DP_Orphans_Maternal, DP_Orphans_AIDS, 1, orpA , pHIVNoARTcHIV + pHIVOnARTcHIV);
        DP.SetOrphansTotal(t, DP_Orphans_Maternal, DP_Orphans_Total, 0, orpA , pNoHIVcNoHIV + pHIVNoARTcNoHIV + pHIVOnARTcNoHIV);
        DP.SetOrphansTotal(t, DP_Orphans_Maternal, DP_Orphans_Total, 1, orpA ,  pHIVNoARTcHIV + pHIVOnARTcHIV);
      end
      else
      begin
        (*Total orphans = NewOrphans + Orphans in age orpA-1 in Total Orphans * Survival to age orpA*)

       //Total uninfected orphans from uninfected mothers
        DP.SetOrphansTotal(t, DP_Orphans_Maternal, DP_Orphans_NON_AIDS, 0, orpA,
                          pNoHIVcNoHIV + DP.GetOrphansTotal(t-1, DP_Orphans_Maternal, DP_Orphans_NON_AIDS, 0, orpA - 1) *
                         (Get_Survival_Rate(DP, tProxy, DP_male, orpA-1)+Get_Survival_Rate(DP, tProxy, DP_Female, orpA-1))/2);

        //Total uninfected orphans from infected mothers
         DP.SetOrphansTotal(t, DP_Orphans_Maternal, DP_Orphans_AIDS, 0, orpA,
                           pHIVNoARTcNoHIV + pHIVOnARTcNoHIV +
                          (DP.GetOrphansTotal(t-1, DP_Orphans_Maternal, DP_Orphans_AIDS, 0, orpA - 1) *
                          (Get_Survival_Rate(DP, tProxy, DP_male, orpA-1)+Get_Survival_Rate(DP, tProxy, DP_Female, orpA-1))/2));

        tempVal := (Get_Survival_Rate(DP, tProxy, DP_male, orpA-1)+Get_Survival_Rate(DP, tProxy, DP_Female, orpA-1))/2;

        tempA := orpA;
        if orpA >14 then
          tempA := 14;
        //if t > 2 then
        //begin
        if DP.GetHIVChildSingleAge(max(t-1,1),tempA-1) > 0 then
          tempVal := DP.GetHIVChildSingleAge(tProxy,tempA)/DP.GetHIVChildSingleAge(max(t-1,1),tempA-1);
        //end;
        if tempVal > (Get_Survival_Rate(DP, tProxy, DP_male, orpA-1)+Get_Survival_Rate(DP, tProxy, DP_Female, orpA-1))/2 then
          tempVal := (Get_Survival_Rate(DP, tProxy, DP_male, orpA-1)+Get_Survival_Rate(DP, tProxy, DP_Female, orpA-1))/2;

        //Total infected orphans from infected mothers
        DP.SetOrphansTotal(t, DP_Orphans_Maternal, DP_Orphans_AIDS, 1, orpA, pHIVNoARTcHIV + pHIVOnARTcHIV
          + (DP.GetOrphansTotal(t-1, DP_Orphans_Maternal, DP_Orphans_AIDS, 1, orpA - 1) *  tempVal));

        for i := 0 to 1 do
          DP.SetOrphansTotal(t, DP_Orphans_Maternal, DP_Orphans_Total, i, orpA,
                            DP.GetOrphansTotal(t, DP_Orphans_Maternal, DP_Orphans_NON_AIDS, i, orpA) +
                            DP.GetOrphansTotal(t, DP_Orphans_Maternal, DP_Orphans_AIDS, i, orpA));
      end;

      for i := DP_Orphans_Total to DP_Orphans_NON_AIDS  do
      begin
        for j := 0 to 1 do
        begin
        DP.SetOrphansTotal(t, DP_Orphans_Maternal, i,j, DP_AllSingAges,
         DP.GetOrphansTotal(t, DP_Orphans_Maternal, i,j, DP_AllSingAges)
         +DP.GetOrphansTotal(t, DP_Orphans_Maternal, i,j, orpA)) ;
        end;
      end;
    end; //end of orpA
  end; //end of t loop
  Result := true;
end;

function TDPProj.DualOrphanCalcs : Boolean;
var
  t, a,i                    : Integer;
  TotalOrphansAgeLastDay,
  TotalNewOrphansAgeLastDay : array [DP_Maternal_AIDS..DP_Dual_Total] of array of single;
  FemalePrevPre18           : array[DP_First_Index-19..-1] of double;
  FemalePrev18Plus          : array of double;
  FPrevDen, FPrevNum,
  Temp2, Temp3, Temp4,
  ma,pa,Value               : double;
  C                         : array [0..17,DP_Total_AIDS_NonAIDS..DP_Non_AIDS] of array of real;
begin
  //Result := false;
  for i := DP_Maternal_AIDS to DP_Dual_Total do
  begin
    SetLength(TotalOrphansAgeLastDay[i], GB.GetNumDynYears(p));
    SetLength(TotalNewOrphansAgeLastDay[i], GB.GetNumDynYears(p));
  end;

  SetLength(FemalePrev18Plus, GB.GetNumDynYears(p));

  for a := 0 to 17 do
    for i := DP_Total_AIDS_NonAIDS to DP_Non_AIDS do
      SetLength(C[a,i], GB.GetNumDynYears(p));

  for t := GB.GetCalcYearIdx(p) to DP.GetFinalIndex do
    for i := DP_Maternal_AIDS to DP_Dual_Total do
    begin
      TotalOrphansAgeLastDay[i,t] := 0.0;
      TotalNewOrphansAgeLastDay[i,t] := 0.0;
    end;

  for t := DP_First_Index - 19 to -1 do
    FemalePrevPre18[t] := 0;

  for t := GB.GetCalcYearIdx(p) to DP.GetFinalIndex do
  begin
    FPrevDen := 0;
    FPrevNum := 0;
    for a := DP_A15_19 to DP_A60_64 do
    begin
      FPrevDen := FPrevDen + DP.GetPop(t, DP_Female, a);
      FPrevNum := FPrevNum + DP.GetHIV(t, DP_Female, a) ;
    end;
    if FPrevDen = 0 then
      FemalePrev18Plus[t] := 0
    else
      FemalePrev18Plus[t] := FPrevNum /  FPrevDen;
  end;

  for t := GB.GetCalcYearIdx(p) - 1 downto GB.GetCalcYearIdx(p) - 18 do
  begin
    if t >= 0 then
      FemalePrev18Plus[t] := FemalePrev18Plus[t+1]/100
    else if t = -1 then
      FemalePrevPre18[t] := FemalePrev18Plus[t+1]/100
    else
      FemalePrevPre18[t] := FemalePrevPre18[t+1]/100;
  end;

    (**************** DUAL AIDS ORPHANS ************************)
  (* Calculate coeffiecients *)
  for t := GB.GetCalcYearIdx(p) to GB.GetCalcYearIdx(p) + 4 do
    for a := 0 to 17 do
    begin
      C[a, DP_Total_AIDS_NonAIDS, t] := 0;
      C[a, DP_aids, t] := 0;
      C[a, DP_non_aids, t] := 0;
    end;

  for t := GB.GetCalcYearIdx(p)  to DP.GetFinalIndex do
  begin
    for a := 0 to 17 do
    begin
      C[a, DP_Total_AIDS_NonAIDS, t] := 0;
      C[a, DP_non_aids, t] := 0;
      if t > GB.GetCalcYearIdx(p) then//+ 4 then
      begin
        for i := 0 to 9 do
        begin
          if (t = GB.GetCalcYearIdx(p) )then //+ 5) then
          begin
            if (a <= 14) then
            begin
             { Temp3 = Female Prevalence 12 years earlier }
              if t > 12 then Temp3 := FemalePrev18Plus[t-12]
                       else Temp3 := 0;

             { Temp2 = expected dual if risk of mother and
                        father dying are independent }
              Temp2 := ((DP.GetOrphans(DP_Paternal_AIDS, a, t)+ DP.GetOrphans(DP_Paternal_non_AIDS, a, t)) /
                       DP.GetBigPop(t, DP_BothSexes, a)) * ((DP.GetOrphans(DP_Maternal_Aids, a, t)
                       + DP.GetOrphans(DP_Maternal_non_AIDS, a, t)) /DP.GetBigPop(t, DP_BothSexes, a));

              if temp2 < 0.0000001 then
                temp2 := 0.0000001;
              C[a, DP_Total_AIDS_NonAIDS, t] := C[a, DP_Total_AIDS_NonAIDS, t]
                                   + (exp(1.3478
                                   + 0.7539 * ln(Temp2 * 100)
                                   - 0.02108 * (a + 0.5 + i/10)
                                   - 0.1506 * ln(max(0.5,(FemalePrev18Plus[t-1] * 100)))
                                   + 0.007225 * (a + 0.5 + i/10) * ln(max(0.5,(FemalePrev18Plus[t-1] * 100)))
                                   + 0.09887 * ln(max(0.5,(Temp3 * 100)))) / 10);

              C[a, DP_non_aids, t] := C[a, DP_Total_AIDS_NonAIDS, t];
            end
            else
            begin
              (* Set value to previous age *)
              (* Done only for ages 15-17 *)
              C[a, DP_Total_AIDS_NonAIDS, t] := C[a-1, DP_Total_AIDS_NonAIDS, t];
              C[a, DP_non_aids, t] := C[a-1, DP_non_aids, t];
            end;
          end
          else
          begin
            if (a <= 14) then
            begin
             { Temp3 = Female Prevalence 12 years earlier }
              if t > 12 then
                Temp3 := FemalePrev18Plus[t-12]
              else
                Temp3 := 0;

             { Temp2 = expected dual orphans if risk of mother and
                        father dying are independent }
              if DP.GetBigPop(t, DP_BothSexes, a) = 0  then
                Temp2 := 0
              else
                Temp2 := ((DP.GetOrphans(DP_Paternal_AIDS, a, t) + DP.GetOrphans(DP_Paternal_non_AIDS, a, t)) /
                         DP.GetBigPop(t, DP_BothSexes, a))* ((DP.GetOrphans(DP_Maternal_Aids, a, t) +
                         DP.GetOrphans(DP_Maternal_non_AIDS, a, t)) / DP.GetBigPop(t, DP_BothSexes, a));
              if temp2 <= 0 then
                temp2 := 0.00000001;

              C[a, DP_Total_AIDS_NonAIDS, t] := C[a, DP_Total_AIDS_NonAIDS, t]
                                   + (exp(1.3478
                                   + 0.7539 * ln(Temp2 * 100)
                                   - 0.02108 * (a + 0.5 + i/10)
                                   - 0.1506 * ln(max(0.5,(FemalePrev18Plus[t-1] * 100)))
                                   + 0.007225 * (a + 0.5 + i/10) * ln(max(0.5,(FemalePrev18Plus[t-1] * 100)))
                                   + 0.09887 * ln(max(0.5,(Temp3 * 100)))) / 10);

              if (DP.GetOrphans(DP_Paternal_AIDS, a, t) + DP.GetOrphans(DP_Maternal_Aids, a, t)) < 10
              then
                C[a, DP_non_aids, t] := C[a, DP_Total_AIDS_NonAIDS, t]
              else
                C[a, DP_non_aids, t] := C[a, DP_non_aids, t-1];
            end
            else
            begin
              (* Set value to previous age *)
              (* Done only for ages 15-17 *)
              C[a, DP_Total_AIDS_NonAIDS, t] := C[a-1, DP_Total_AIDS_NonAIDS, t];
              C[a, DP_non_aids, t] := C[a-1, DP_non_aids, t];
            end;
          end; (* END OF IF *)
        end; (* END OF I LOOP *)
      end; (* END OF IF STATEMENT *)
    end; (* END OF A LOOP *)
  end; (* END OF T LOOP *)

  for t := GB.GetCalcYearIdx(p)  to DP.GetFinalIndex do
  begin
    for a := 0 to 17 do
    begin
      if c[a, DP_Total_AIDS_NonAIDS, t] <> 0 then
      begin
        (* Calculate Total Dual *)
       Temp2 := c[a,DP_Total_AIDS_NonAIDS,t] / 100 * DP.GetBigPop(t, DP_BothSexes, a);
       DP.SetOrphans(DP_Dual_Total, a, t,temp2);

       (* Calc Total Dual *)
       TotalOrphansAgeLastDay[DP_Dual_Total,t] := TotalOrphansAgeLastDay[DP_Dual_Total,t] +
                                                  DP.GetOrphans(DP_Dual_Total, a, t);
      end;

      if c[a, DP_non_aids, t] <> 0 then
      begin
      (* Calculate Dual Non AIDS *)
       (* if DP.GetBigPop(t, DP_BothSexes, a) = 0 then
          Temp3 := 0
        else
          Temp3 := c[a, DP_non_aids, t] *
                 DP.GetAgeLastBDay(t, DP_Paternal_Non_AIDS,a) *
                DP.GetAgeLastBDay(t, DP_Maternal_Non_AIDS,a)
                 /DP.GetBigPop(t, DP_BothSexes, a);*)

        if  (DP.GetOrphans(DP_Maternal_Orphans, a, t) * DP.GetOrphans(DP_Paternal_Orphans, a, t))  = 0 then
          Temp3 := 0
        else
        (*Population X Risk_MaternalNonAIDS_Orphans X  Risk_PaternalNonAIDS_Orphans  X Increased_Risk_Dual_Orphans*)
          Temp3 := DP.GetBigPop(t, DP_BothSexes, a) * DP.GetOrphans(DP_Maternal_Non_AIDS, a, t) *
                   DP.GetOrphans(DP_Paternal_Non_AIDS, a, t) * (c[a,DP_Total_AIDS_NonAIDS,t] / 100 ) /
                   (DP.GetOrphans(DP_Maternal_Orphans, a, t) * DP.GetOrphans(DP_Paternal_Orphans, a, t)) ;

          DP.SetOrphans(DP_Dual_Non_AIDS, a, t,temp3);

         (* Calc Total Dual Non AIDS *)
          TotalOrphansAgeLastDay[DP_Dual_Non_AIDS,t] := TotalOrphansAgeLastDay[DP_Dual_Non_AIDS,t] +
                                                        DP.GetOrphans(DP_Dual_Non_AIDS, a, t);
      end;

      if (DP.GetOrphans(DP_Dual_Total, a, t) - DP.GetOrphans(DP_Dual_Non_AIDS, a, t)) >= 0 then
      begin
        (* Calculate Dual AIDS *)
        temp4 := DP.GetOrphans(DP_Dual_Total, a, t) - DP.GetOrphans(DP_Dual_Non_AIDS, a, t);
        if temp4 < 0 then
          temp4 := 0;
        if temp4 > DP.GetOrphans(DP_Paternal_AIDS, a, t) then
          temp4 := DP.GetOrphans(DP_Paternal_AIDS, a, t);
        if temp4 > DP.GetOrphans(DP_Maternal_AIDS, a, t) then
          temp4 := DP.GetOrphans(DP_Maternal_AIDS, a, t);

        DP.SetOrphans(DP_Dual_AIDS, a, t,temp4);

        (* Calc Total Dual AIDS *)
        TotalOrphansAgeLastDay[DP_Dual_AIDS,t] := TotalOrphansAgeLastDay[DP_Dual_AIDS,t] +
                                                  DP.GetOrphans(DP_Dual_AIDS, a, t);
     end;

      if (DP.GetOrphans(DP_Dual_Total, a, t) - DP.GetOrphans(DP_Dual_Non_AIDS, a, t)) >= 0 then
      begin
        (* New Dual AIDS Orphans *)
        if t = GB.GetCalcYearIdx(p) then //DP_FIRST_INDEX  then
          DP.SetNewOrphans(DP_Dual_AIDS, a, t, 0)
        else
        begin
          if a = 0 then
            DP.SetNewOrphans(DP_Dual_AIDS, a, t, DP.GetOrphans(DP_Dual_AIDS, a, t))
          else
          begin
            (* Calc Average of infected + Uninfected / Total Orphans and apply *)
            (* to Dual *)
            if DP.GetOrphans(DP_Maternal_AIDS, a, t) <> 0 then
              ma := (DP.GetOrphansTotal(t,DP_Orphans_Maternal,DP_Orphans_AIDS,0,a) +
                     DP.GetOrphansTotal(t,DP_Orphans_Maternal,DP_Orphans_AIDS,1,a)) /
                     DP.GetOrphans(DP_Maternal_AIDS, a, t)
            else
              ma := 0;
            pa := ma;

            (* NewOrphans(t,a,type) = Orphans(t,a,type) - orphans(t-1,a-1,type) *
            average * (Get_Survival_Rate(p,t,DP_male,a)+ Get_Survival_Rate(p,t,DP_Female,a))/2 *)

            Value := DP.GetOrphans(DP_Dual_AIDS, a, t) - DP.GetOrphans(DP_Dual_AIDS, a-1, t-1) *
                     (ma+pa)/2 *(Get_Survival_Rate(DP,t,DP_male,a) + Get_Survival_Rate(DP,t,DP_Female,a))/2;
            DP.SetNewOrphans(DP_Dual_AIDS, a, t, Value);
          end;

          TotalNewOrphansAgeLastDay[DP_Dual_AIDS,t] := TotalNewOrphansAgeLastDay[DP_Dual_AIDS,t] +
                                                       DP.GetNewOrphans(DP_Dual_AIDS, a, t);
        end;

        (* New Dual NON-AIDS Orphans *)

        if t = GB.GetCalcYearIdx(p) then //DP_FIRST_INDEX  then
          DP.SetNewOrphans(DP_Dual_Non_AIDS, a, t, 0)
        else
        begin
          if a = 0 then
            DP.SetNewOrphans(DP_Dual_Non_AIDS, a, t, DP.GetOrphans(DP_Dual_Non_AIDS, a, t))
          else
          begin
            (* NewOrphans(t,a,type) = Orphans(t,a,type) - orphans(t-1,a-1,type) *
            ( Get_Survival_Rate(p,t,DP_male,a)+ Get_Survival_Rate(p,t,DP_Female,a))/2 *)
            Value := DP.GetOrphans(DP_Dual_Non_AIDS, a, t) -
                     DP.GetOrphans(DP_Dual_Non_AIDS, a-1, t-1) *
                    (Get_Survival_Rate(DP,t,DP_male,a) + Get_Survival_Rate(DP,t,DP_Female,a))/2;
            DP.SetNewOrphans(DP_Dual_Non_AIDS, a, t, Value);
          end;
        end;

        TotalNewOrphansAgeLastDay[DP_Dual_Non_AIDS,t] := TotalNewOrphansAgeLastDay[DP_Dual_Non_AIDS,t] +
                                                         DP.GetNewOrphans(DP_Dual_Non_AIDS, a, t);

        (* New Dual Total Orphans *)
        temp4 := DP.GetNewOrphans(DP_Dual_AIDS, a, t) + DP.GetNewOrphans(DP_Dual_Non_AIDS, a, t);
        DP.SetNewOrphans(DP_Dual_Total, a, t, Temp4);

        TotalNewOrphansAgeLastDay[DP_Dual_Total,t] := TotalNewOrphansAgeLastDay[DP_Dual_Total,t] +
                                                      DP.GetNewOrphans(DP_Dual_Total, a, t);
      end;
    end;(* END OF A LOOP *)
  end; (* END OF T LOOP *)

  (* Dual Orphans *)

  for t := GB.GetCalcYearIdx(p)  to DP.GetFinalIndex do
  begin
    (* Dual AIDS *)
    DP.SetTotalOrphans(DP_DualTotal, DP_AIDS,t, TotalOrphansAgeLastDay[DP_Dual_AIDS,t]);

    (* Dual Non AIDS *)
    DP.SetTotalOrphans(DP_DualTotal, DP_Non_AIDS,t, TotalOrphansAgeLastDay[DP_Dual_Non_AIDS,t]);

    (* Dual Total *)
    DP.SetTotalOrphans(DP_DualTotal, DP_Total_AIDS_NonAIDS,t, TotalOrphansAgeLastDay[DP_Dual_Total,t]);

    (* New Dual AIDS *)
    //DP.SetTotalNewOrphans(DP_DualTotal, DP_AIDS,t, TotalNewOrphansAgeLastDay[DP_Dual_AIDS,t]);

    (* New Dual Non AIDS *)
    //DP.SetTotalNewOrphans(DP_DualTotal, DP_Non_AIDS,t, TotalNewOrphansAgeLastDay[DP_Dual_Non_AIDS,t]);

    (* New Dual Total *)
    //DP.SetTotalNewOrphans(DP_DualTotal, DP_Total_AIDS_NonAIDS,t, TotalNewOrphansAgeLastDay[DP_Dual_Total,t]);
  end;

 Result := true;
 for i := DP_Maternal_AIDS to DP_Dual_Total do
  begin
    Finalize(TotalOrphansAgeLastDay[i]);
    Finalize(TotalNewOrphansAgeLastDay[i]);
  end;

  Finalize(FemalePrev18Plus);

 for a := 0 to 17 do
  for i := DP_Total_AIDS_NonAIDS to DP_Non_AIDS do
    Finalize(C[a,i]);
end;

function TDPProj.TotalOrphanCalcs : Boolean;
var
  a, t                         : Byte;
  temp, temp2, temp3, value    : Double;
  p5                           : FiveAgeArray;
  sa                           : SingleAgeArray;
begin
  Result := false;
  (* Calc Total New Orphans *)

  for t := GB.GetCalcYearIdx(p)  to DP.GetFinalIndex do
  begin
    if t = GB.GetCalcYearIdx(p) then //DP_FIRST_INDEX  then
      for a := 0 to DP_MAX_AGE do
        DP.SetNewOrphans(DP_Total_orphans, a, t,0)
    else
    begin
      (* Calc single age for AIDS deaths *)
      for a := 0 to DP_MAX_AGE do
        p5[a] := DP.GetAIDSDeaths(t, a, Dp_BothSexes);
      Calc_Single_Pop_Ages(p5, sa);

      temp := 0;
      for a := 0 to 17 do
      begin
        if a = 0 then
        begin
          DP.SetNewOrphans(DP_Total_orphans, a, t,
                            DP.GetAgeLastBDay(t, DP_MaternalAIDS_ALB, a) +
                            DP.GetAgeLastBDay(t, DP_MaternalNonAIDS_ALB, a) +
                            DP.GetAgeLastBDay(t, DP_PaternalAIDS_ALB, a) +
                            DP.GetAgeLastBDay(t, DP_PaternalNONAIDS_ALB, a) -
                            DP.GetOrphans(DP_Dual_Total, a,t));
          temp := Temp + DP.GetNewOrphans(DP_Total_orphans, a, t);
        end
        else
        begin
          (* NewOrphans(t,a,type) = Orphans(t,a,type) - orphans(t-1,a-1,type) *
          ( Get_Survival_Rate(p,t,DP_male,a)+ Get_Survival_Rate(p,t,DP_Female,a))/2 *)
          temp2 := DP.GetAgeLastBDay(t, DP_MaternalAIDS_ALB, a) +
                   DP.GetAgeLastBDay(t, DP_MaternalNonAIDS_ALB, a) +
                   DP.GetAgeLastBDay(t, DP_PaternalAIDS_ALB, a) +
                   DP.GetAgeLastBDay(t, DP_PaternalNONAIDS_ALB, a) -
                   DP.GetOrphans(DP_Dual_Total, a,t);

          temp3 := DP.GetAgeLastBDay(t-1, DP_MaternalAIDS_ALB, a-1) +
                   DP.GetAgeLastBDay(t-1, DP_MaternalNonAIDS_ALB, a-1) +
                   DP.GetAgeLastBDay(t-1, DP_PaternalAIDS_ALB, a-1) +
                   DP.GetAgeLastBDay(t-1, DP_PaternalNONAIDS_ALB, a-1) -
                   DP.GetOrphans(DP_Dual_Total, a-1,t-1);

          Value := temp2 - temp3 * (Get_Survival_Rate(DP,t,DP_male,a)+
                                    Get_Survival_Rate(DP,t,DP_Female,a))/2 +
                                    sa[a] * 0.5;
          DP.SetNewOrphans(DP_Total_orphans, a, t, Value);
          temp := temp + value;
        end;
        DP.SetAllNewOrphans(t,DP_TotalOrphans,temp);
      end; (* End of A Loop *)
    end;
  end;

  for t := GB.GetCalcYearIdx(p) to DP.GetFinalIndex do
  begin
    (* All AIDS Orphans *)
    temp := Get_DPDefined_Age(DP, t, DP_BothSexes, 0, 17);
    if temp <> 0 then
      Value :=  (DP.GetTotalOrphans(DP_MaternalTotal, DP_AIDS, t)/ temp +
                 DP.GetTotalOrphans(DP_PaternalTotal, DP_AIDS, t)/ temp -
                 DP.GetTotalOrphans(DP_DualTotal, DP_AIDS, t)/ temp +
                 DP.GetTotalOrphans(DP_PaternalTotal, DP_non_AIDS, t) / temp *
                 DP.GetTotalOrphans(DP_MaternalTotal, DP_AIDS, t)/ temp +
                 DP.GetTotalOrphans(DP_MaternalTotal, DP_non_AIDS, t) / temp *
                 DP.GetTotalOrphans(DP_PaternalTotal, DP_AIDS, t)/ temp ) * temp
    else
      Value := 0;
    DP.SetAllOrphans(t,DP_AllAIDSOrphans,value);

    (* All Non-AIDS Orphans *)
    Value := (DP.GetTotalOrphans(DP_MaternalTotal, DP_Non_AIDS, t)+
              DP.GetTotalOrphans(DP_PaternalTotal, DP_Non_AIDS, t)-
              DP.GetTotalOrphans(DP_DualTotal, DP_Non_AIDS, t));
    DP.SetAllOrphans(t,DP_AllNonAIDSOrphans,value);

    (* Total Orphans *)
    Value := (DP.GetTotalOrphans(DP_MaternalTotal, DP_Total_AIDS_NonAIDS, t)+
              DP.GetTotalOrphans(DP_PaternalTotal, DP_Total_AIDS_NonAIDS, t)-
              DP.GetTotalOrphans(DP_DualTotal, DP_Total_AIDS_NonAIDS, t));
    DP.SetAllOrphans(t,DP_TotalOrphans,value);
  end; (* END OF T LOOP *)
end;

(*Estimates transition time between CD4 counts for each age *)
procedure TDPProj.CreateInfDur (t, s:Byte);
var
  age, j, a                   : Byte;
  i                           : Integer;
  value, sum                  : Double;
  addValues                   : array [0..1000] of Double;//a matrix to simulate 1000 exponentially distributed random number to estimate transition time
  RandVal                     : extended;
begin

  for age := 15 to DP_MaxSingleAges do
  begin
    a := DP_CD4_15_24;
    case age of
      15..24               : a := DP_CD4_15_24;
      25..34               : a := DP_CD4_25_34;
      35..44               : a := DP_CD4_35_44;
      45..DP_MaxSingleAges : a := DP_CD4_45_54;
    end;

    for i := 0 to 1000 do
      addValues[i] := 0;

    for j := DP_CD4_GT500 to DP_CD4_LT50 do
    begin
      sum:= 0;
      for i := 0 to 1000 do
      begin
        RandVal := random;
        if RandVal = 0 then
          addValues[i] := 0
        else if j <  DP_CD4_LT50   then
          addValues[i] := addValues[i] + (-ln(RandVal) * DP.GetAdultAvgDurByCD4(DP_Data, s, a, j) )
        else
          addValues[i] := addValues[i] + (-ln(RandVal) * DP.GetAdultMortByCD4NoART(DP_Data, s, a, DP_CD4_LT50) );
        sum:= sum +   addValues[i]  ;
      end;

      value := age + sum / 1000;
      case round(value) of
        15..24                : a := DP_CD4_15_24;
        25..34                : a := DP_CD4_25_34;
        35..44                : a := DP_CD4_35_44;
        45..DP_MaxSingleAges  : a := DP_CD4_45_54;
      end;

      durMatrix[s, age - 15 , j] := round(value);
    end;
  end;
end;


function TDPProj.getAgeAtInfection (s, A, CD4 : Byte): Byte  ;
var
  i       : Byte;
  found   : Boolean;
begin
  found := false;
  i := 0;
  while (i <= DP_MaxSingleAges - 14) and (NOT found) do
  begin
    if durMatrix[s, i, CD4] >= A then
      found := true;
    i := i + 1;
  end;
  Result:=  i - 2 + 15;
end;

procedure TDPProj.setDurARTIni(s : Byte) ;
var
  t, c                : integer;
  sumProdMatrix       : array[DP_CD4_GT500 .. DP_CD4_LT50] of array of Double;
  cumSumART           : array[DP_CD4_GT500 .. DP_CD4_LT50] of array of Double;
begin
  for c := DP_CD4_GT500 to DP_CD4_LT50 do
  begin
    SetLength(sumProdMatrix[c], GB.GetNumDynYears(p));
    SetLength(cumSumART[c], GB.GetNumDynYears(p));
  end;

  for c := DP_CD4_GT500 to DP_CD4_LT50 do
    for t := DP_FIRST_INDEX  to DP.GetFinalIndex do
      durARTMatrix[s,c,t] :=  0;

  for c := DP_CD4_GT500 to DP_CD4_LT50 do
  begin
    for t := DP_FIRST_INDEX  to DP.GetFinalIndex do
    begin;
      sumProdMatrix[c,t] := 0;
      durARTMatrix[s,c,t] :=  0;
      cumSumART[c,t] := 0;
    end;
    if DP.GetOnART(DP_FIRST_INDEX, s, DP_AllAges) > 0 then
      durARTMatrix[s, c, DP_FIRST_INDEX] := round (DP.GetOnART(DP_FIRST_INDEX, s, DP_AllAges)/ DP.GetOnART(DP_FIRST_INDEX, s, DP_AllAges));

    cumSumART[c,DP_FIRST_INDEX] := DP.GetOnART(DP_FIRST_INDEX, s, DP_AllAges);
  end;

  for t := DP_FIRST_INDEX + 1 to DP.GetFinalIndex do
  begin
    for c :=   DP_CD4_GT500 to DP_CD4_LT50 do
    begin
      sumProdMatrix[c,t] :=  sumProdMatrix[c,t-1] + DP.GetAIDSDeathsCD4OnART(t-1, DP_AllSingAges, s, c, DP_D_All) * durARTMatrix[s,c,t-1];

      cumSumART[c,t] :=  cumSumART[c,t-1] + DP.GetCD4Distribution(t, s, c, DP_OnART);

      durARTMatrix[s, c, t] := 0;
      if DP.GetCD4Distribution(t, s, c, DP_OnART) > 0 then
        durARTMatrix[s,c,t] := round(max(0,(cumSumART[c,t] - sumProdMatrix[c,t])) /  DP.GetCD4Distribution(t, s, c, DP_OnART));
    end;
  end;

  for c := DP_CD4_GT500 to DP_CD4_LT50 do
  begin
    Finalize(sumProdMatrix[c]);
    Finalize(cumSumART[c]);
  end;
end;

procedure TDPProj.CalcBFtransmission(p, t, m1, m2: integer);
var
  d,r                 : byte;
//  Percent_NoProph,
  Percent_OptA,
  Percent_OptB,
//  Percent_ART,
//  DropOut_OptA,
//  DropOut_OptB,
  DropOut_ART,
  Excess            : double;
begin
  BFTR := 0;
  Percent_OptA := TreatPercent[DP_OptA_BF];
  Percent_OptB := TreatPercent[DP_OptB_BF];

{ Use drop out rate for <12 months if m1 is < 7, since m1 is 2-month periods, otherwise use 12+ drop out rate}
  if (m1 < 7) then
    DropOut_ART := DP.GetARVRegimen(DP_AnnDropPostnatalProph, DP_ART0_12MthsBF, DP_Percent, t)
  else
    DropOut_ART := DP.GetARVRegimen(DP_AnnDropPostnatalProph, DP_ARTGT12MthsBF, DP_Percent, t);
//  DropOut_ART := DP.GetARVRegimen(DP_AnnDropPostnatalProph, DP_ART0_12MthsBF, DP_Percent, t);
  { Adjust transmission rates for Opt A and Opt B if
    some have CD4 counts under 350 }
  OptA_Transmission_Rate := DP.GetTransEffAssump(DP_Data, DP_OptionA, DP_BreastfeedingGE350);
  OptB_Transmission_Rate := DP.GetTransEffAssump(DP_Data, DP_OptionB, DP_BreastfeedingGE350);
  if PropGE350 > 0 then
  begin
    if (Percent_OptA + Percent_OptB - TreatPercent[DP_TripleARTBefPreg]
          - TreatPercent[DP_TripleARTDurPreg]) > PropGE350 then
    begin
      Excess := (Percent_OptA + Percent_OptB - TreatPercent[DP_TripleARTBefPreg]
               - TreatPercent[DP_TripleARTDurPreg]) - PropGE350;
      OptA_Transmission_Rate :=
        (PropGE350 * DP.GetTransEffAssump(DP_Data, DP_OptionA,
          DP_BreastfeedingGE350) + Excess * 1.45 / 0.46 * DP.GetTransEffAssump
          (DP_Data, DP_OptionA, DP_BreastfeedingGE350)) / (PropGE350 + Excess);
      OptB_Transmission_Rate :=
        (PropGE350 * DP.GetTransEffAssump(DP_Data, DP_OptionB,
          DP_BreastfeedingGE350) + Excess * 1.45 / 0.46 * DP.GetTransEffAssump
          (DP_Data, DP_OptionB, DP_BreastfeedingGE350)) / (PropGE350 + Excess);
    end;
  end;

  for d := m1 to m2 do
  begin
    { Initialize no prophylaxis to HIV+ PW who have not transmitted }
    TreatPercent[DP_NoProphylaxis] := 1 - PTR3/100 - BFTR;
    { For those on ART }
    for r := DP_SingleDoseNev to DP_ARTStartDurPreg_Late do
    begin
      if d > 1 then TreatPercent[r-3] := TreatPercent[r-3] * (1 - DropOut_ART / 100 * 2);
      BFTR := BFTR + TreatPercent[r-3] * DP.GetTransEffAssump(DP_Data, r,DP_BreastfeedingLT350) / 100  * 2
                     * (1 - DP.GetInfantFeedingOptions(d, DP_InPMTCT,t) / 100);
      TreatPercent[DP_NoProphylaxis] := TreatPercent[DP_NoProphylaxis] - TreatPercent[r-3];
    end;
    { For those not on ART }
    if TreatPercent[DP_NoProphylaxis] < 0 then TreatPercent[DP_NoProphylaxis] := 0;
    BFTR := BFTR + TreatPercent[DP_NoProphylaxis]
                   * (1 - DP.GetInfantFeedingOptions(d, DP_NotInPMTCT,t) / 100)
                   * ( PropLT350 * DP.GetTransEffAssump(DP_Data, DP_NoProphylaxis,DP_BreastfeedingLT350) / 100 * 2
                      + (1-PropLT350)* DP.GetTransEffAssump(DP_Data, DP_BreastfeedingGE350,DP_BreastfeedingGE350) / 100 * 2);
    { Exclude transmission during first 6 weeks since it is already included in the perinatal rate }
    if d = 1 then BFTR := BFTR /4;
  end;
end; {CalcBFtransmission}

Destructor TDPProj.Destroy;
var
  a, c, s : byte;
begin
  Finalize(TotalOldAIDS);

  for a := DP_A0 to DP_A14 do
    for c := 1 to DP_CD4_MaxCategories do
    begin
      Finalize(Count_dist[a,c]);
      Finalize(Percent_dist[a,c]);
    end;

  for s := DP_Male to DP_Female do
    for c := DP_CD4_GT500 to DP_CD4_LT50 do
      Finalize(durARTMatrix[s,c]);
end;

(* Capture pop at year prior calc year and set to DPCalcStateData variable.  This code
   is only used when the calc state data for a country is being set. It is assumed that
   the projection to determine the calc state pop values starts at 1970 and goes to
   the final year of the primary projection.*)
procedure TDPProj.CaptureCalcStateData(pop1 : DP_TPopType);
var
  a,s,c,d,r         : byte;
  CalcStateData     : DP_TCalcStateObject;
  CaptureYearIdx    : integer;
  AllZeroPop        : boolean;
  PopAdj            : array [0..DP_MaxSingleAges,DP_Male..DP_Female] of double;
begin
  CalcStateData := GBSession.GetProjection(GB.GetPrimaryProj).DP.GetCalcStateData;
 (* CalcStateIdx is equal to the first year of the primary projection *)
  CaptureYearIdx := GBTempProjIdx(GB, GB_FirstYrIdxOfPrimProj);

  AllZeroPop := GBSession.GetProjection(GB.GetPrimaryProj).DP.GetCalcYearPopALlZero;

  if (GB.GetInsertUsersPopBoolean) and (not AllZeroPop) then
  begin
    for a := 0 to DP_MaxSingleAges do
      for s := DP_Male to DP_Female do
      begin
        (* We need to know the adjustment made to the calculated population so we can
           apply the same adjustment to the population and all HIV levels at time t = 1
           so that births and deaths can be calculated properly (births, for example, uses
           the current year and prior year divided by 2). The values in the pop1 array at
           this point are the calculated values in the temporary projection based on the
           user population from the First Year Population editor (which are most likely
           different than what the user put in the editor since they are calculated) and
           BigPop has actual the user entered value from the first year population editor for t = 2 *)
        if pop1[CaptureYearIdx+1,s,a,DP_H_All,DP_D_All] = 0 then
          PopAdj[a,s] := 1
        else
          PopAdj[a,s] := GBSession.GetProjection(GB.GetPrimaryProj).DP.GetBigPop(2,s,a)
                         / pop1[CaptureYearIdx+1,s,a,DP_H_All,DP_All];
      end;
  end
  else
  begin
    for a := 0 to DP_MaxSingleAges do
    for s := DP_Male to DP_Female do
      PopAdj[a,s] := 1;
  end;

  for a := 0 to DP_MaxSingleAgesChildren do
    for s := DP_Male to DP_Female do
      for c:= DP_H_All to DP_H_MAX do
        for d := DP_D_All to DP_D_Max do
          CalcStateData.CPop[a,s,c,d] := pop1[CaptureYearIdx,s,a,c,d] * PopAdj[a,s];

  for a := 15 to DP_MaxSingleAges do
    for s := DP_Male to DP_Female do
      for c:= DP_H_All to DP_H_MAX do
        for d := DP_D_All to DP_D_Max do
          CalcStateData.Pop[a,s,c,d] := pop1[CaptureYearIdx,s,a,c,d] * PopAdj[a,s];

  for a := DP_A0 to DP_A14 do
    for c := 1 to DP_CD4_MaxCategories do
    begin
      CalcStateData.CountDist[a,c] := Count_Dist[a,c,CaptureYearIdx];
      CalcStateData.PercentDist[a,c] := Percent_Dist[a,c,CaptureYearIdx];
    end;

  for a := DP_AllAges to DP_MAX_AGE do
    for s := DP_BothSexes to DP_Female do
    begin
      if s = DP_BothSexes then
      begin
        CalcStateData.Births[a,s] := (DP.GetBirths(CaptureYearIdx, a, DP_Male) * PopAdj[a,DP_Male]) +
                                     (DP.GetBirths(CaptureYearIdx, a, DP_Female) * PopAdj[a,DP_Female]);

        CalcStateData.Deaths[a,s] := (DP.GetDeaths(CaptureYearIdx, a, DP_Male) * PopAdj[a,DP_Male]) +
                                     (DP.GetDeaths(CaptureYearIdx, a, DP_Female) * PopAdj[a,DP_Female]);

        CalcStateData.AIDSDeaths[a,s] := (DP.GetAIDSDeaths(CaptureYearIdx, a, DP_Male) * PopAdj[a,DP_Male]) +
                                         (DP.GetAIDSDeaths(CaptureYearIdx, a, DP_Female) * PopAdj[a,DP_Female]);

        CalcStateData.AIDSDeathsART[a,s] := (DP.GetAIDSDeathsART(CaptureYearIdx, a, DP_Male) * PopAdj[a,DP_Male]) +
                                            (DP.GetAIDSDeathsART(CaptureYearIdx, a, DP_Female) * PopAdj[a,DP_Female]);

        CalcStateData.AIDSDeathsNoART[a,s] := (DP.GetAIDSDeathsNoART(CaptureYearIdx, a, DP_Male) * PopAdj[a,DP_Male]) +
                                              (DP.GetAIDSDeathsNoART(CaptureYearIdx, a, DP_Female) * PopAdj[a,DP_Female]);

        CalcStateData.NewInfections[a,s] := (DP.GetNewInfections(CaptureYearIdx, a, DP_Male) * PopAdj[a,DP_Male]) +
                                            (DP.GetNewInfections(CaptureYearIdx, a, DP_Female) * PopAdj[a,DP_Female]) ;
      end
      else
      begin
        CalcStateData.Births[a,s] := DP.GetBirths(CaptureYearIdx, a, s) * PopAdj[a,s];
        CalcStateData.Deaths[a,s] := DP.GetDeaths(CaptureYearIdx, a, s) * PopAdj[a,s];
        CalcStateData.AIDSDeaths[a,s] := DP.GetAIDSDeaths(CaptureYearIdx, a, s) * PopAdj[a,s];
        CalcStateData.AIDSDeathsART[a,s] := DP.GetAIDSDeathsART(CaptureYearIdx, a, s) * PopAdj[a,s];
        CalcStateData.AIDSDeathsNoART[a,s] := DP.GetAIDSDeathsNoART(CaptureYearIdx, a, s) * PopAdj[a,s];
        CalcStateData.NewInfections[a,s] := DP.GetNewInfections(CaptureYearIdx, a, s) * PopAdj[a,s];
      end;
    end;

  CalcStateData.AdultDeaths := DP.GetAdultDeaths(CaptureYearIdx);

  for s := DP_BothSexes to DP_Female do
    for a := DP_AllAGes to DP_MAX_AGE + 1 do
    begin
      if s = DP_BothSexes then
      begin
        CalcStateData.NewlyNeedART[s,a] := (DP.GetNewlyNeedART(CaptureYearIdx, DP_Male, a) * PopAdj[a,DP_Male]) +
                                           (DP.GetNewlyNeedART(CaptureYearIdx, DP_Female, a) * PopAdj[a,DP_Female]);

        CalcStateData.UnmetNeedART[s,a] := (DP.GetUnmetNeedART(CaptureYearIdx, DP_Male, a) * PopAdj[a,DP_Male]) +
                                           (DP.GetUnmetNeedART(CaptureYearIdx, DP_Female, a) * PopAdj[a,DP_Female]);

        CalcStateData.NeedART[s,a] := (DP.GetNeedART(CaptureYearIdx, DP_Male, a) * PopAdj[a,DP_Male]) +
                                      (DP.GetNeedART(CaptureYearIdx, DP_Female, a) * PopAdj[a,DP_Female]);

        CalcStateData.OnART[s,a] := (DP.GetOnART(CaptureYearIdx, DP_Male, a) * PopAdj[a,DP_Male]) +
                                    (DP.GetOnART(CaptureYearIdx, DP_Female, a) * PopAdj[a,DP_Female]);
      end
      else
      begin
        CalcStateData.NewlyNeedART[s,a] := DP.GetNewlyNeedART(CaptureYearIdx, s, a) * PopAdj[a,s];
        CalcStateData.UnmetNeedART[s,a] := DP.GetUnmetNeedART(CaptureYearIdx, s, a) * PopAdj[a,s];
        CalcStateData.NeedART[s,a] := DP.GetNeedART(CaptureYearIdx, s, a) * PopAdj[a,s];
        CalcStateData.OnART[s,a] := DP.GetOnART(CaptureYearIdx, s, a) * PopAdj[a,s];
      end;
    end;

  for s := DP_BothSexes to DP_Female do
  begin
    CalcStateData.AdjLE[s] := DP.GetAdjLE(CaptureYearIdx, s);
    CalcStateData.HIVPopUnder1[s] := DP.GetHIVPopUnder1(CaptureYearIdx, s);
    CalcStateData.ChildCTXNeed1To4[s] := DP.GetChildCTXNeed1To4(CaptureYearIdx, s);
    CalcStateData.ChildARTNeed1To4[s] := DP.GetChildARTNeed1To4(CaptureYearIdx, s);
    CalcStateData.HIVPop1517[s] := DP.GetHIVPop1517(CaptureYearIdx, s);
    CalcStateData.ARTNeedUnder2[s] := DP.GetARTNeedUnder2(CaptureYearIdx, s);
  end;

  for s := DP_Male to DP_Female do
    CalcStateData.ARTTwoYearsAgoBySex[s] := ARTTwoYearsAgoBySexCalcState[s];

  for r := DP_NumNeedCot to DP_NumOnArt do
    for s := DP_BothSexes to DP_Female do
      CalcStateData.ChildARTCalc[r,s] := DP.GetChildARTCalc(r,s,CaptureYearIdx);

  for r := DP_PerChildHIVPosCot to DP_PerChildHIVRecART do
    CalcStateData.ChildTreatInputs[r] := DP.GetChildTreatInputs(r,CaptureYearIdx);

  for a := 0 to DP_MaxSingleAges do
    for s := DP_Male to DP_Female do
    begin
      CalcStateData.DeathsBySingleAge[s,a] := DP.GetDeathsBySingleAge(CaptureYearIdx,a,s) * PopAdj[a,s];
      CalcStateData.AIDSDeathsBySingleAge[s,a] := DP.GetAIDSDeathsBySingleAge(CaptureYearIdx,a,s) * PopAdj[a,s];
    end;

  CalcStateData.OnPMTCT := DP.GetOnPMTCT(CaptureYearIdx);
  CalcStateData.NeedPMTCT := DP.GetNeedPMTCT(CaptureYearIdx);
  CalcStateData.MTCTRate6Wks := DP.GetMTCTRate6Wks(CaptureYearIdx);
  CalcStateData.NewInfantInfections := DP.GetNewInfantInfections(CaptureYearIdx);
  CalcStateData.HIVPregWomen := DP.GetHIVPregWomen(CaptureYearIdx);
  CalcStateData.HIVPregWomenCD4LT350 := DP.GetHIVPregWomenCD4LT350(CaptureYearIdx);

  (* Need survival rates for 1970 *)
  for s := DP_Male to DP_Female do
    for a := 1 to 5 do
      CalcStateData.SurvRate[s,a] := DP.GetSurvRate(DP_FIRST_INDEX,s,a);

  (* Need births the calcyear from the primary projection and the 0-4 population
     of the temporary projection to make sure TFR is accurate when using
     the wizard to create a subnational projection *)
  for s := DP_Male to DP_Female do
  begin
    CalcStateData.NatBirths[s] := DP.GetBirths(CaptureYearIdx + 1, DP_AllAges, s) {* PopAdj[DP_AllAges,s]};

    for a := DP_A0 to DP_A4 do
      CalcStateData.NatPop0_4[s] := CalcStateData.NatPop0_4[s] +
                                    DP.GetBigPop(CaptureYearIdx + 1,s,a-1);
  end;

  (* CalcState was estabilshed above by the temporary projection to be used in the
    primary projection.  Now we need to transfer the CalcState to the primary projection. *)
  GBSession.GetProjection(GB.GetPrimaryProj).DP.SetCalcStateData(CalcStateData);

  (* if the population for t = 2 is all zeros, then this is the first time that CalcState is
     being done for the projection (this is a new projection) so we need to populate BigPop
     for t = 2 so the user will see values in the First Year Popultion editor. We need to display
     the values for year 2, the CalcYear in the editor. *)
 (* Save pop by single age *)
 if AllZeroPop then
 begin
  for a := 0 to DP_MaxSingleAges do
    for s := DP_BothSexes to DP_Female do
    begin
      if s = DP_BothSexes then
      begin
        GBSession.GetProjection(GB.GetPrimaryProj).DP.SetBigPop(1, DP_BothSexes, a,
                    pop1[CaptureYearIdx,DP_Male,a,DP_H_All,DP_D_All]+pop1[CaptureYearIdx,DP_Female,a,DP_H_All,DP_D_All]);
        GBSession.GetProjection(GB.GetPrimaryProj).DP.SetBigPop(2, DP_BothSexes, a,
                    pop1[CaptureYearIdx+1,DP_Male,a,DP_H_All,DP_D_All]+pop1[CaptureYearIdx+1,DP_Female,a,DP_H_All,DP_D_All]);
      end
      else
      begin
        GBSession.GetProjection(GB.GetPrimaryProj).DP.SetBigPop(1, s, a, pop1[CaptureYearIdx,s,a,DP_H_All,DP_D_All]);
        GBSession.GetProjection(GB.GetPrimaryProj).DP.SetBigPop(2, s, a, pop1[CaptureYearIdx+1,s,a,DP_H_All,DP_D_All]);
      end;
    end; {s loop}
 end
 else
 begin
   for a := 0 to DP_MaxSingleAges do
    for s := DP_BothSexes to DP_Female do
    begin
      if s = DP_BothSexes then
      begin
        GBSession.GetProjection(GB.GetPrimaryProj).DP.SetBigPop(1, DP_BothSexes, a,
                    (pop1[CaptureYearIdx,DP_Male,a,DP_H_All,DP_D_All] * PopAdj[a,DP_Male]) +
                    (pop1[CaptureYearIdx,DP_Female,a,DP_H_All,DP_D_All] * PopAdj[a,DP_Female]));
      end
      else
      begin
        GBSession.GetProjection(GB.GetPrimaryProj).DP.SetBigPop(1, s, a, pop1[CaptureYearIdx,s,a,DP_H_All,DP_D_All] * PopAdj[a,s]);
      end;
    end; {s loop}
 end;
end;

(* Used to retrieve the pop values for the calcstate from the values previously captured *)
procedure TDPProj.RetrieveCalcStateData(pop1 : DP_TPopType);
var
  a,r,s,c,d,a1         : byte;
  CalcStateData        : DP_TCalcStateObject;
  temp1                : double;
const
  FirstYearIdx = 1;
begin
  CalcStateData := DP.GetCalcStateData;

  for a := 0 to DP_MaxSingleAgesChildren do
    for s := DP_BothSexes to DP_Female do
      for c:= DP_H_All to DP_H_MAX do
        for d := DP_D_All to DP_D_Max do
          if s = DP_BothSexes then
            pop1[FirstYearIdx,s,a,c,d] := CalcStateData.CPop[a,DP_Male,c,d] + CalcStateData.CPop[a,DP_Female,c,d]
          else
            pop1[FirstYearIdx,s,a,c,d] := CalcStateData.CPop[a,s,c,d];

  for a := 15 to DP_MaxSingleAges do
    for s := DP_BothSexes to DP_Female do
      for c:= DP_H_All to DP_H_MAX do
        for d := DP_D_All to DP_D_Max do
        begin
          if s = DP_BothSexes then
            pop1[FirstYearIdx,s,a,c,d] := CalcStateData.Pop[a,DP_Male,c,d] + CalcStateData.Pop[a,DP_Female,c,d]
          else
            pop1[FirstYearIdx,s,a,c,d] := CalcStateData.Pop[a,s,c,d];
        end;

  for a := DP_A0 to DP_A14 do
    for c := 1 to DP_CD4_MaxCategories do
    begin
      Count_Dist[a,c,FirstYearIdx] := CalcStateData.CountDist[a,c];
      Percent_Dist[a,c,FirstYearIdx] := CalcStateData.PercentDist[a,c];
    end;

  for a := DP_AllAges to DP_MAX_AGE do
    for s := DP_BothSexes to DP_Female do
    begin
      DP.SetBirths(FirstYearIdx,a,s,CalcStateData.Births[a,s]);
      DP.SetDeaths(FirstYearIdx,a,s,CalcStateData.Deaths[a,s]);
      DP.SetAIDSDeaths(FirstYearIdx,a,s,CalcStateData.AIDSDeaths[a,s]);
      DP.SetAIDSDeathsART(FirstYearIdx,a,s,CalcStateData.AIDSDeathsART[a,s]);
      DP.SetAIDSDeathsNoART(FirstYearIdx,a,s,CalcStateData.AIDSDeathsNoART[a,s]);
      DP.SetNewInfections(FirstYearIdx,a,s,CalcStateData.NewInfections[a,s]);
    end;

  DP.SetAdultDeaths(FirstYearIdx,CalcStateData.AdultDeaths);

  for s := DP_BothSexes to DP_Female do
    for a := DP_AllAGes to DP_MAX_AGE + 1 do
    begin
      DP.SetNewlyNeedART(FirstYearIdx,s,a,CalcStateData.NewlyNeedART[s,a]);
      DP.SetUnmetNeedART(FirstYearIdx,s,a,CalcStateData.UnmetNeedART[s,a]);
      DP.SetNeedART(FirstYearIdx,s,a,CalcStateData.NeedART[s,a]);
      DP.SetOnART(FirstYearIdx,s,a,CalcStateData.OnART[s,a]);
    end;

  for s := DP_BothSexes to DP_Female do
  begin
    DP.SetAdjLE(FirstYearIdx,s,CalcStateData.AdjLE[s]);
    DP.SetHIVPopUnder1(FirstYearIdx,s,CalcStateData.HIVPopUnder1[s]);
    DP.SetChildCTXNeed1To4(FirstYearIdx,s,CalcStateData.ChildCTXNeed1To4[s]);
    DP.SetChildARTNeed1To4(FirstYearIdx,s,CalcStateData.ChildARTNeed1To4[s]);
    DP.SetHIVPop1517(FirstYearIdx,s,CalcStateData.HIVPop1517[s]);
    DP.SetARTNeedUnder2(FirstYearIdx,s,CalcStateData.ARTNeedUnder2[s]);
  end;

  for s := DP_Male to DP_Female do
    ARTTwoYearsAgoBySexCalcState[s] := CalcStateData.ARTTwoYearsAgoBySex[s];

  for r := DP_NumNeedCot to DP_NumOnArt do
    for s := DP_BothSexes to DP_Female do
      DP.SetChildARTCalc(r,s,FirstYearIdx,CalcStateData.ChildARTCalc[r,s]);

  for r := DP_PerChildHIVPosCot to DP_PerChildHIVRecART do
    DP.SetChildTreatInputs(r,FirstYearIdx,CalcStateData.ChildTreatInputs[r]);

  for a := 0 to DP_MaxSingleAges do
    for s := DP_Male to DP_Female do
    begin
      DP.SetDeathsBySingleAge(FirstYearIdx,a,s,CalcStateData.DeathsBySingleAge[s,a]);
      DP.SetAIDSDeathsBySingleAge(FirstYearIdx,a,s,CalcStateData.AIDSDeathsBySingleAge[s,a]);
    end;

  DP.SetOnPMTCT(FirstYearIdx,CalcStateData.OnPMTCT);
  DP.SetNeedPMTCT(FirstYearIdx,CalcStateData.NeedPMTCT);
  DP.SetMTCTRate6Wks(FirstYearIdx,CalcStateData.MTCTRate6Wks);
  DP.SetNewInfantInfections(FirstYearIdx,CalcStateData.NewInfantInfections);
  DP.SetHIVPregWomen(FirstYearIdx,CalcStateData.HIVPregWomen);
  DP.SetHIVPregWomenCD4LT350(FirstYearIdx,CalcStateData.HIVPregWomenCD4LT350);

  for s := DP_Male to DP_Female do
    for a := 1 to 5 do
      DP.SetSurvRate1970(s,a,CalcStateData.SurvRate[s,a]);

  for a := 0 to DP_MaxSingleAges do
    for s := DP_BothSexes to DP_Female do
      if s = DP_BothSexes then
        DP.SetBigPop(FirstYearIdx, DP_BothSexes, a,
                      pop1[FirstYearIdx,DP_Male,a,DP_H_All,DP_D_All]+pop1[FirstYearIdx,DP_Female,a,DP_H_All,DP_D_All])
      else
        DP.SetBigPop(FirstYearIdx, s, a, pop1[FirstYearIdx,s,a,DP_H_All,DP_D_All]);

  (* Calculate HIV population for AIM State *)
  (* Determine HIV Pop for all ages except 80 *)
  for s := DP_Male to DP_Female do
  begin
    for a := DP_A0_4 to DP_A75_79 do
    begin
      temp1 := 0;
      for a1 := (a-1)*5 to (a-1)*5+4 do
      begin
        for c:= DP_CD4_GT500 to DP_CD4_LT50 do
          for d := DP_D_HIVPos to DP_D_Max do
            temp1 := temp1 + pop1[FirstYearIdx,s,a1,c,d];
      end; {a1 loop}
      DP.SetHIV(FirstYearIdx,s,a,temp1);
    end; {a loop}
    (* Add 80 year olds *)
    temp1 := 0;
    for c:= DP_CD4_GT500 to DP_CD4_LT50 do
      for d := DP_D_HIVPos to DP_D_Max do
        temp1 := temp1 + pop1[FirstYearIdx,s,80,c,d];
    DP.SetHIV(FirstYearIdx,s,DP_A80_Up,temp1);
  end; {s loop}
  (* Determine totals for males and females for all ages *)
  for s := DP_Male to DP_Female do
  begin
    temp1 := 0;
    for a := DP_A0_4 to DP_A80_up do
    begin
      DP.SetHIV(FirstYearIdx,DP_BothSexes,a,DP.GetHIV(FirstYearIdx,DP_Male,a)
                                               + DP.GetHIV(FirstYearIdx,DP_Female,a));
      temp1 := temp1 + DP.GetHIV(FirstYearIdx,s,a);
    end; {a loop}
    DP.SetHIV(FirstYearIdx,s,DP_AllAges,temp1);
  end; {s loop}
  (* Determine total for Both Sexes for all ages *)
  DP.SetHIV(FirstYearIdx,DP_BothSexes,DP_AllAges,DP.GetHIV(FirstYearIdx,DP_Male,DP_AllAges)
                                               + DP.GetHIV(FirstYearIdx,DP_Female,DP_AllAges));
end;

{$REGION 'DPPROJA'}
  procedure TDPProj.Do_AIDS_calculations(p, t: integer);
  var
    sr, sum, sum4,
    sum5, sum6, sum7, sum8,
    sum9, sum10, NHB, adultHIVnum,
    NewHIVnum, pop5, pop25, prev25,
    temp_cd4, adj,
    PMTCTFraction, v1, v2, v3, v4,
    v5, v6, AIDSdeaths, CtxPer, ArtPer,
    BothPer, NonePer, TxEff, InfTxEff,
    ChildTxEff, PrevAdj, ARTpercent,
    FractionProgressing, progressors,
    NeedForFirstLineART, NeedForSecondLineART,
    NeedForFirstLineARTmale,
    NeedForFirstLineARTfemale,
    UnmetNeedForFirstLineART,
    UnmetNeedForFirstLineARTmale,
    UnmetNeedForFirstLineARTfemale,
    NeedForFirstLine1To4,
    OnFirstLineArt, OnSecondLineArt,
    NewFLArt,
    OnFirstLineArtmale,
    OnFirstLineArtfemale,
    OnCTX, NewCTX,
    NewInfFromBFGT12,
    UnmetNeedForCTX,
    FYARTsurv, PercentExposed, IdentifiedNeed,
    IdentifiedNeedFLART,
    NeededART, ARTLastYear,
    MonthsOfExposure, ExclusiveBF, MixedBF,
    MaleHIVChildren, FemaleHIVChildren,
    Excess, NotOnTx,
    TotalPopAge, nART,
    MedianDurationBF,
    ArtEff,
    ChildrenOnCTX, CTXnum,
    ProportionAverted                                     : single;
    SumOnART_Age0,
    SumOnART0_4, SumOnART5_9, SumOnART10_14,
    ARTsum, MedCD4, mort, ltfu, peddrop, pedsum,
    value, ART_Excess,
    Number_Eligible, sum_start,
    sum_end, ARTsurvivors,
    Count_eligible, Percent_eligible,
    SumART, SumART_15_49, SumART_BS,
    SumART_15_49_BS,
    NewART, RemainingNewART,
    StartART, totalE,
    NotOnART, OnART                                       : double;
    FeedOpt,
    a, s, d, h, t1, a1, a2, a7, a5,
    ysi, count, d1, ag, a9,
    agroup, r, c1, c2, c, s1                              : byte;
    yr                                                    : integer;
    InFile,
    Outfile                                               : TextFile;
    tstring                                               : string;
    ARTTwoYearsAgoBySex                                   : array[DP_Male..DP_Female] of single;
    FeedPercent,
    FeedPercentIn,
    FeedPercentOut                                        : array[DP_MixFeed_HIVMoth .. DP_RepFeed_HIVMoth] of single;
    temps, temps1, temps2, temps3                         : string;
    AdultsNotOnART                                        : array[DP_Male..DP_Female,DP_CD4_GT500..DP_CD4_LT50] of double;
    EligibleSpecialPop                                    : array[DP_CD4_GT500..DP_CD4_LT50, DP_Male..DP_Female] of double;
    Prop3                                                 : array [DP_CD4_GT500 .. DP_CD4_LT50] of double;
    lambda                                                : array [1 .. 4, DP_CD4_GT500 .. DP_CD4_LT50, DP_Male .. DP_Female] of double;
    plambda                                               : array [DP_CD4_0_2 .. DP_CD4_5_14, DP_CD4_Per_GT30 .. DP_CD4_Per_LT5,
                                                                   DP_Male .. DP_Female, DP_P_Perinatal .. DP_P_BF12] of double;
    mu                                                    : array [1 .. 4, DP_CD4_GT500 .. DP_CD4_LT50, DP_Male .. DP_Female] of double;
    pmu                                                   : array [DP_CD4_0_2 .. DP_CD4_5_14, DP_CD4_Per_GT30 .. DP_CD4_Per_LT5,
                                                                   DP_Male .. DP_Female, DP_P_Perinatal .. DP_P_BF12] of double;
    alpha                                                 : array [1 .. 4, DP_CD4_GT500 .. DP_CD4_LT50, DP_Male .. DP_Female, DP_D_ARTlt6m .. DP_D_ARTgt12m] of double;
    palpha                                                : array [0..14, DP_CD4_Per_GT30 .. DP_CD4_Per_LT5,
                                                                   DP_Male .. DP_Female, DP_D_ARTlt6m .. DP_D_ARTgt12m] of double;
    beta1, beta2                                          : array [1 .. 4, DP_CD4_GT500 .. DP_CD4_LT50, DP_Male .. DP_Female] of double;
    CD_perc_at_infection                                  : array [DP_CD4_Per_GT30 .. DP_CD4_Per_LT5] of double;
    DP_CD4_lower_limit                                    : array [DP_CD4_GT500-1 .. DP_CD4_LT50] of double;
    EligibleAdults                                        : array [DP_BothSexes .. DP_Female] of array of double;
    entrants, exits                                       : array [DP_CD4_GT500 .. DP_CD4_LT50, DP_D_HIVPos .. DP_D_Max] of double;
    NewPatients, Eadults                                  : array [0 .. DP_CD4_LT50, DP_BothSexes .. DP_Female] of double;
    NewChildART                                           : array [0 .. 14, DP_BothSexes .. DP_Female, DP_CD4_Per_GT30 .. DP_CD4_Per_LT5,
                                                                   DP_D_All .. DP_P_BF12] of double;
    Age0to14CD4Distr, NewARTCD4Distr                      : array [DP_CD4_GT500 .. DP_CD4_LT50] of double;
    valueMinMort                                          : array [15 ..DP_MaxSingleAges,  DP_Male .. DP_Female] of double;
    AdultAIDSDeaths                                       : array [15..DP_MaxSingleAges, DP_Male..DP_Female, DP_NoTreat..DP_OnART] of double;
    ARTPropMortOrLTFU                                     : array [DP_Male..DP_Female, 15..DP_MaxSingleAges, DP_CD4_GT500 .. DP_CD4_LT50, DP_D_ARTlt6m .. DP_D_ARTgt12m] of double;
    age10                                                 : array [15..DP_MaxSingleAges] of byte;
    ChildAIDSDeathsOnART                                  : array [0..14, DP_Male..DP_Female] of double;
    ChildTiming                                           : array [DP_CD4_Per_GT30..DP_CD4_Per_LT5, DP_P_Perinatal..DP_P_BF12] of double;
  const
    DP_TimeStep = 10;

    procedure CalcInfectionsByAge(proj : byte; t, t1 : integer);
    const
      dt = 1.0 / DP_TimeStep;
    var
      a, a5, b, s, c, d, q, _minDebut          : byte;
      hiv15p, inciLastYr, prevLastYr,
      propSuppressed, effectCondom,
      effectART, effectCircum, _MaleCircumicision: double;
      condoms,condoms_t0,condoms_nu,
      contact, contact_t0, contactFactors      : double;
      RiskFactor, RiskFactor_t0,RiskFactor_t   : double;
      ImpactFactor,ImpactFactor_t0,ImpactFactor_t,ImpactFactor_t0_t:double;
      violence,violence_t0,IPV_RR_condoms:double;

      //tally : array[DP_Male..DP_Female] of double;

      _SexDebut                                : array[DP_Male..DP_Female] of integer;
      effectPrEP, effectPrEP_t0                : array[DP_Male..DP_Female,DP_15_24..DP_45_54] of double;
      dstPrEPFactor,dstPrEPFactor_t0,
      dstFactors,
      _ProbAcquire, _STICoFac                  : array[DP_Male..DP_Female] of double;
      _PrevHasSTI                              : array[DP_A15_19.. DP_A45_49,DP_Male..DP_Female] of double;
      _POP, _HIVNeg,
      srcFactor, _ForceInfect, _NewInfAge      : array[DP_Male..DP_Female,15..DP_MaxSingleAges] of double;
      force_infect:double;
      _HIVPos, _HIVART, stageMult, propPrimary : array[DP_Male..DP_Female,15..DP_MaxSingleAges,DP_CD4_GT500..DP_CD4_LT50] of double;
      age05                                    : array[15..DP_MaxSingleAges] of byte;

      SexDebutDelayFactor,InterGenSexReducFactor:double;
      t0:byte;
      R_P,R_I:array[DP_Male..DP_Female,DP_15_24..DP_45_54] of double;
      KeyPopFactors:double;

      HA : THAData;
      total_newinf : double; //to be removed after testing
    begin
      HA := GetHAData(proj);
      t0 := HA.GetBehaveChangeBaseYearIdx;

      if (t1 = 1) then
      begin
        // We only need to calculate transmission factors for condom use
        // on the first time step of each year, since these remain constant
        // throughout the year
        HAObj.CalcCondomFactor(t);
      end;

      //Biomedical intervention scale-up is controlled by
      //HA.GetApplyViralLoadSupp
      //HA.GetApplyPrEPSusceptibleReduc
      //HA.GetApplyMaleCircumcision
      //HA.GetApplyViralLoadSupp

      //Prioritization by location for PrEP is controlled by
      //GetGeoPrioritization
      //and for Behav and Struct Interventions in:
      //HA_BehaveChangeFactorCalc(proj : byte);

      //correction for geographic-prioritization
      for s := DP_Male to DP_Female do
      for a := DP_15_24 to DP_45_54 do
      begin
        R_P[s,a]:=0;//pop propr, not relative prop in this use
        R_I[s,a]:=1;
      end;

      if(HA.GetGeoPrioritization=TRUE) then
      begin
         //for s := DP_Male to DP_Female do
         for s in [DP_Female] do
         for a := DP_15_24 to DP_45_54 do
         begin
            R_P[s,a]:=HA.GetAgePopSize(a, t)/100;
            R_I[s,a]:=HA.GetAgePopIncRR(a, t);
         end;
      end;
      //////////////////


      (* Initialize direct address tables for looking up age groups based on
      single year of age *)
      for a := 15 to 19 do age05[a] := DP_A15_19;
      for a := 20 to 24 do age05[a] := DP_A20_24;
      for a := 25 to 29 do age05[a] := DP_A25_29;
      for a := 30 to 34 do age05[a] := DP_A30_34;
      for a := 35 to 39 do age05[a] := DP_A35_39;
      for a := 40 to 44 do age05[a] := DP_A40_44;
      for a := 45 to DP_MaxSingleAges do age05[a] := DP_A45_49;

      _SexDebut[DP_Male ] := floor(HA.GetSexDebut(DP_Male ));
      _SexDebut[DP_Female] := floor(HA.GetSexDebut(DP_Female));
      _minDebut := Min(_SexDebut[DP_Male], _SexDebut[DP_Female]);

      for s := DP_Male to DP_Female do
        for a := 15 to DP_MaxSingleAges do
        begin
          inciLastYr := DP.GetNewInfectionsBySingleAge(t-1,s,a+1);
          for c := DP_CD4_GT500 to DP_CD4_LT50 do
          begin
            prevLastYr := pop1[t-1,s,a,c,DP_D_All];
            if (prevLastYr > 0) then
            begin
              propPrimary[s,a,c] := (DP.GetAdultDistNewInfectionsCD4(DP_Data, s, age10[a], c) / 100)
                                  * (HA.GetMonthsInPrimaryStage / 12.0)
                                  * inciLastYr / prevLastYr;
              if (propPrimary[s,a,c] > 1.0) then propPrimary[s,a,c] := 1.0;
            end
            else propPrimary[s,a,c] := 0.0;
          end;

          for c := DP_CD4_GT500 to DP_CD4_200_249 do
            stageMult[s,a,c] := HA.GetInfectiousness(DP_HA_PrimaryInfectionStage) * propPrimary[s,a,c]
                              + HA.GetInfectiousness(DP_HA_AsymptomaticStage) * (1.0 - propPrimary[s,a,c]);
          for c := DP_CD4_100_199 to DP_CD4_LT50 do
            stageMult[s,a,c] := HA.GetInfectiousness(DP_HA_PrimaryInfectionStage) * propPrimary[s,a,c]
                              + HA.GetInfectiousness(DP_HA_SymptomaticStage) * (1.0 - propPrimary[s,a,c]);
        end;

      effectCondom := HA.GetCondomEffect / 100;
      effectCircum := HA.GetCircumcisionEffect / 100;
      effectART := HA.GetViralSuppressionEffect / 100;

      // arbitrarily introduce HIV to the 16-26 year old age group at epidemic start year
      if ((t = HA.GetEpiYearIdx) and (t1 = 1)) then
      begin
        for a := 16 to 26 do
          for s := DP_Male to DP_Female do
          begin
            pop1[t,s,a,DP_CD4_GT500,DP_D_HIVPos] := pop1[t,s,a,DP_CD4_GT500,DP_D_HIVPos] + 1;
            pop1[t,s,a,DP_H_No, DP_D_HIVNeg] := pop1[t,s,a,DP_H_No, DP_D_HIVNeg] - 1;
          end;
      end;

      (* Aggregate the adult population *)
      for s := DP_Male to DP_Female do
        for a := 15 to DP_MaxSingleAges do
        begin
          _HIVNeg[s,a] := pop1[t,s,a,DP_H_No,DP_D_HIVNeg];
          for c := DP_CD4_GT500 to DP_CD4_LT50 do
          begin
            (* People who have been on ART for <6 months are counted in HIVPos
               instead of HIVART because we assume they are not virologically
               suppressed yet *)
            _HIVPos[s,a,c] := pop1[t,s,a,c,DP_D_HIVPos] + pop1[t,s,a,c,DP_D_ARTlt6m];
            _HIVART[s,a,c] := pop1[t,s,a,c,DP_D_ART6to12m] + pop1[t,s,a,c,DP_D_ARTgt12m];
          end;

          _POP[s,a] := _HIVNeg[s,a];
          for c := DP_CD4_GT500 to DP_CD4_LT50 do
            _POP[s,a] := _POP[s,a] + _HIVPos[s,a,c] + _HIVART[s,a,c];
        end;

      for s := DP_Male to DP_Female do
        for a := 15 to DP_MaxSingleAges do
          _NewInfAge[s,a] := 0.0;

      if (t1 = 1) then
      begin
        for s := DP_Male to DP_Female do
        begin
          for a := DP_A0_4 to DP_A80_up do
            DP.SetNewInfections(t, a, s, 0);
          for a := DP_A0 to DP_A80 do
            DP.SetNewInfectionsBySingleAge(t, s, a, 0);
        end;
      end;

      (* Exit early if there are no HIV+ individuals *)
      hiv15p := 0.0;
      for s := DP_Male to DP_Female do
        for a := 15 to DP_MaxSingleAges do
          hiv15p := hiv15p + (_POP[s,a] - _HIVNeg[s,a]);
      if (hiv15p <= 0) then exit;

      (* Calculate force of infection and new infections by age *)
      total_newinf := 0;
      HA_PartnerRateCalc(proj, t);

      for s := DP_Male to DP_Female do
      begin
        _STICoFac[s] := HA.GetSTICoFac(s);
        _ProbAcquire[s]:= HA.GetHIVAcquisitionProb(s) / 100;
      end;

      for s := DP_Male to DP_Female do
      for a := DP_A15_19 to DP_A45_49 do
        _PrevHasSTI[a,s] := HA.GetHasSTI(a, s, t) / 100;

      (* Initialize PrEP incidence reductions *)
      for s := DP_Male to DP_Female do
        for a := DP_15_24 to DP_45_54 do
        begin
          effectPrEP_t0[s,a] := 0.0;
          effectPrEP[s,a] := 0.0;
        end;


      if(t>=t0) then
      for a := DP_15_24 to DP_45_54 do
      begin
        s := DP_Male;
        for q in [DP_HA_PrEP_Oral, DP_HA_PrEP_Inject] do
        begin
             effectPrEP[s,a] := effectPrEP[s,a]
                             + HA.GetPrEPCoverage(s,a,q,t) / 100
                             * HA.GetPrEPEfficacy(s,a,q)/ 100
                             * HA.GetPrEPAdherence(s,a,q,t) / 100;

            if((HA.GetApplyPrEPSusceptibleReduc=FALSE)and(t>=t0)) then
              effectPrEP[s,a] := effectPrEP[s,a]
                               + HA.GetPrEPCoverage(s,a,q,t0) / 100
                               * HA.GetPrEPEfficacy(s,a,q)/ 100
                               * HA.GetPrEPAdherence(s,a,q,t0) / 100;
        end;

        s := DP_Female;
        for q := DP_HA_PrEP_FirstMethod to DP_HA_PrEP_FinalMethod do
        begin
             effectPrEP[s,a] := effectPrEP[s,a]
                              + HA.GetPrEPCoverage(s,a,q,t) / 100
                              * HA.GetPrEPEfficacy(s,a,q)/ 100
                              * HA.GetPrEPAdherence(s,a,q,t) / 100;

            if((HA.GetApplyPrEPSusceptibleReduc=FALSE)and(t>=t0)) then
              effectPrEP[s,a] := effectPrEP[s,a]
                              + HA.GetPrEPCoverage(s,a,q,t0) / 100
                              * HA.GetPrEPEfficacy(s,a,q)/ 100
                               * HA.GetPrEPAdherence(s,a,q,t0) / 100;
        end;

        for q in [DP_HA_PrEP_Oral, DP_HA_PrEP_Inject] do
        begin
             effectPrEP_t0[s,a] := effectPrEP_t0[s,a]
                             + HA.GetPrEPCoverage(s,a,q,t0) / 100
                             * HA.GetPrEPEfficacy(s,a,q)/ 100
                             * HA.GetPrEPAdherence(s,a,q,t0) / 100;
        end;

        s := DP_Female;
        for q := DP_HA_PrEP_FirstMethod to DP_HA_PrEP_FinalMethod do
        begin
             effectPrEP_t0[s,a] := effectPrEP_t0[s,a]
                              + HA.GetPrEPCoverage(s,a,q,t0) / 100
                              * HA.GetPrEPEfficacy(s,a,q)/ 100
                              * HA.GetPrEPAdherence(s,a,q,t0) / 100;
        end;

      end;

      for s := DP_Male to DP_Female do
        for a := 15 to DP_MaxSingleAges do
        begin
          _ForceInfect[s,a] := 0.0;
          //recorded for inter-gen-sex display in HA
          for b := 15 to DP_MaxSingleAges do
            HA.SetForceInfect(s,a,b,t,0);
        end;


      (* Cache the transmission term based on HIV+ source partner factors *)
      for a := 15 to DP_MaxSingleAges do
      begin
        for s := DP_Male to DP_Female do
        begin

          propSuppressed := HA.GetARTViralSuppression(s, age10[a], t) / 100;
          if((HA.GetApplyViralLoadSupp=FALSE)and(t>=t0)) then
              propSuppressed := HA.GetARTViralSuppression(s, age10[a], t0) / 100;

          srcFactor[s,a] := 0.0;
          for c := DP_CD4_GT500 to DP_CD4_LT50 do
          begin
            srcFactor[s,a] := srcFactor[s,a]
                            + _HIVPos[s,a,c] * stageMult[s,a,c]
                            + _HIVART[s,a,c] * stageMult[s,a,c] * (1.0 - propSuppressed)
                            + _HIVART[s,a,c] * (1.0 - effectART) * propSuppressed;
          end;
        end;
        srcFactor[DP_Male,a] := _ProbAcquire[DP_Female] * srcFactor[DP_Male,a];
        srcFactor[DP_Female,a] := _ProbAcquire[DP_Male] * srcFactor[DP_Female,a];
      end;

      for a := _minDebut to DP_MaxSingleAges do
      begin

        _MaleCircumicision := HA.GetMaleCircumByAge(age05[a], t) / 100;
        if((HA.GetApplyMaleCircumcision=FALSE)and(t>=t0)) then
          _MaleCircumicision := HA.GetMaleCircumByAge(age05[a], t0) / 100;

        { Factors that affect HIV acquisition risk regardless of partner age }
        dstPrEPFactor[DP_Female] := (1.0 - effectPrEP[DP_Female,age10[a]]);
        dstFactors[DP_Female]      := (1.0 + _PrevHasSTI[age05[a],DP_Female] * _STICoFac[DP_Female]);

        dstPrEPFactor[DP_Male] := (1.0 - effectPrEP[DP_Male,age10[a]]);
        dstFactors[DP_Male] := (1.0 - _MaleCircumicision * effectCircum)
                             * (1.0 + _PrevHasSTI[age05[a],DP_Male] * _STICoFac[DP_Male]);

        { Factors that affect HIV acquisition risk regardless of partner age }
        //need at time t0 for GeoPrior
        dstPrEPFactor_t0[DP_Female] := (1.0 - effectPrEP_t0[DP_Female,age10[a]]);
        dstPrEPFactor_t0[DP_Male] := (1.0 - effectPrEP_t0[DP_Male,age10[a]]);

        { Sum force of infection acting on HIV- females from HIV+ males }
        if (a >= _SexDebut[DP_Female]) then
        begin

          //correction for delay in sexual debut
          SexDebutDelayFactor:= 1;
          if(a<_SexDebut[DP_Female]+HA.GetBehaveChangeDebutDelayFactor(DP_Female, t)) then
          begin
            if(_SexDebut[DP_Female]+HA.GetBehaveChangeDebutDelayFactor(DP_Female, t) - a >= 1) then
              SexDebutDelayFactor:= 0
            else
              SexDebutDelayFactor:= frac(HA.GetBehaveChangeDebutDelayFactor(DP_Female, t));
          end;

          for b := _SexDebut[DP_Male] to DP_MaxSingleAges do // Loop over HIV+ partners' ages
          begin
            // Inner loop. These calculations repeat 43,560 times per year

            // correction for reduction in intergenerational sex for AGYW
            if ((15 <= a) and (a <= 24))
            then InterGenSexReducFactor := HA.GetBehaveChangePartnerChoiceFactor(a, b, t)
            else InterGenSexReducFactor := 1;

            if _POP[DP_Male,b] > 0 then
            begin

              //at time t0
              violence_t0 := HA.GetPartnerViolence(age10[a], t0) / 100;
              condoms_t0 := HA.GetCondomByAge(age10[b], t0) / 100;
              contact_t0 := HA.GetPartnerChoice(a, b);
              contact_t0 := contact_t0 * HA.GetPartnerAgeTrend(a);

              RiskFactor_t0:= (1.0 - condoms_t0 * effectCondom) * dstPrEPFactor_t0[DP_Female];

              condoms := HA.GetCondomByAge(age10[b], t) / 100;
              if((HA.GetApplyCondomUse=FALSE)and(t>=t0)) then
                condoms := HA.GetCondomByAge(age10[b], t0) / 100;
              condoms := 1.0 - (1.0 - condoms) * HA.GetBehaveChangeCondomsFactor(age10[a],age10[b],t);

              //Effect of violence prevention, by acting on reducing violence, which reduces condom non-use
              //conduse = 1 - ((1-p)*c_nu + p*c_nu*RR),
              //where c_nu is condom non-use, p % experience violence and RR for condon non-use for them
              if(t>t0) then
              begin
                violence_t0 := HA.GetPartnerViolence(age10[a], t0) / 100;
                violence := max(0,violence_t0 * HA.GetBehaveChangeViolenceFactor(age10[a], t));
                HA.SetPartnerViolence(age10[a], t, violence*100 );

                IPV_RR_condoms:=HA.GetIPVCondomNonUseRR(age10[a], t);
                condoms_nu:=1-condoms;
                if(violence<violence_t0) then
                  condoms := condoms + max(0,0*(violence-violence_t0)*(condoms_nu) + (violence_t0-violence)*(condoms_nu)*IPV_RR_condoms);

                condoms := min(condoms,1);
              end;

              if(b-a>=5) then//inter gen sex reduction with men five years older
                contact := HA.GetPartnerChoice(a, b) * InterGenSexReducFactor
              else
                contact := HA.GetPartnerChoice(a, b);
              contact := contact * HA.GetPartnerAgeTrend(a) * SexDebutDelayFactor;
              contact := contact * HA.GetBehaveChangePartnerFactor(age10[a],age10[b],t);

              RiskFactor_t:= (1.0 - condoms * effectCondom)* dstPrEPFactor[DP_Female];
              contactFactors:=  contact * dstFactors[DP_Female];
              //Key population factors
              KeyPopFactors := HA.GetKPBehaveChangeCondomsFactor(HA_KP_FSW, age10[a], t)
                             * HA.GetKPBehaveChangeCondomsFactor(HA_KP_FIDU, age10[a], t)
                             * HA.GetKPIDUNeedleShareFactor(HA_KP_FIDU, age10[a], t)
                             * HA.GetKPPrEPFactor(HA_KP_FSW, age10[a], t)
                             * HA.GetKPPrEPFactor(HA_KP_FIDU, age10[a], t);

             RiskFactor:=RiskFactor_t;

             if((HA.GetGeoPrioritization=TRUE) and (t>t0) and (R_P[DP_Female,age10[a]]>0)) then
             begin
               //Estimate the fraction of impact that can be attributed to prioritizing to high risk groups
               //incidence ~ (1-p)*RiskFactor + p*RiskFactor*RR),
               //e.g. p % AGYW in high incidence districts

               //This correction approach in meant to represent all interventions, i.e. the effect of
               //prioritization is not done seperately on each interventions

               //baseline risk
               Impactfactor_t0 := RiskFactor_t0*(1-R_P[DP_Female,age10[a]])+RiskFactor_t0*R_P[DP_Female,age10[a]]*R_I[DP_Female,age10[a]];
               //reduced risk following implementation of interventions to all populations
               ImpactFactor_t  := RiskFactor_t*(1-R_P[DP_Female,age10[a]])+RiskFactor_t*R_P[DP_Female,age10[a]]*R_I[DP_Female,age10[a]];
               //reduced risk following implementation of interventions to high-risk sub-population
               Impactfactor_t0_t := RiskFactor_t0*(1-R_P[DP_Female,age10[a]])+RiskFactor_t*R_P[DP_Female,age10[a]]*R_I[DP_Female,age10[a]];

               ImpactFactor:=1;
               if(Impactfactor_t0-Impactfactor_t <> 0) then
               begin
                  RiskFactor:=RiskFactor_t0;
                  if(Impactfactor_t0-Impactfactor_t <> 0) then
                  begin
                    ImpactFactor:= (Impactfactor_t0-Impactfactor_t0_t)/(Impactfactor_t0-Impactfactor_t);
                    if( (RiskFactor_t0>RiskFactor_t) and (ImpactFactor>=0) and (ImpactFactor<=1)) then
                      RiskFactor:= RiskFactor_t0-ImpactFactor*(RiskFactor_t0-RiskFactor_t);
                  end;
               end;

             end;

              force_infect:=(ContactFactors * RiskFactor * KeyPopFactors) * srcFactor[DP_Male,b] / _POP[DP_Male,b];
              _ForceInfect[DP_Female,a] := _ForceInfect[DP_Female,a] + force_infect;

              if(t1=DP_TimeStep) then
                HA.SetForceInfect(DP_Female,a,b,t,(1-exp(-dt*force_infect)));

            end;
          end;
        end;

        { Sum force of infection acting on HIV- males from HIV+ females }
        if (a >= _SexDebut[DP_Male]) then
        begin

          //correction for delay in sexual debut
          SexDebutDelayFactor:= 1;
          if(a<_SexDebut[DP_Male]+HA.GetBehaveChangeDebutDelayFactor(DP_Male, t)) then
          begin
          if(_SexDebut[DP_Male]+HA.GetBehaveChangeDebutDelayFactor(DP_Male, t) - a >= 1) then
            SexDebutDelayFactor:= 0
          else
            SexDebutDelayFactor:= frac(HA.GetBehaveChangeDebutDelayFactor(DP_Male, t));
          end;

          if _POP[DP_Male,a] > 0 then
          begin
            for b := _SexDebut[DP_Female] to DP_MaxSingleAges do // Loop over HIV+ partners' ages
            begin
              condoms := HA.GetCondomByAge(age10[a], t) / 100;
              if((HA.GetApplyCondomUse=FALSE)and(t>=t0)) then
                condoms := HA.GetCondomByAge(age10[a], t0) / 100;
              condoms := 1.0 - (1.0 - condoms) * HA.GetBehaveChangeCondomsFactor(age10[b],age10[a],t);

              contact := HA.GetPartnerChoice(b, a);
              contact := contact * HA.GetPartnerAgeTrend(b) * SexDebutDelayFactor;
              contact := contact * HA.GetBehaveChangePartnerFactor(age10[b],age10[a],t);
              //Key population factors
              KeyPopFactors :=  HA.GetKPBehaveChangeCondomsFactor(HA_KP_MSM, age10[a], t)
                              * HA.GetKPBehaveChangeCondomsFactor(HA_KP_MIDU, age10[a], t)
                              * HA.GetKPIDUNeedleShareFactor(HA_KP_MIDU, age10[a], t)
                              * HA.GetKPPrEPFactor(HA_KP_MSM, age10[a], t)
                              * HA.GetKPPrEPFactor(HA_KP_MIDU, age10[a], t);

              force_infect:= contact * KeyPopFactors * dstPrEPFactor[DP_Male] * dstFactors[DP_Male] * (1.0 - condoms * effectCondom) * srcFactor[DP_Female,b] / _POP[DP_Male,a];
              _ForceInfect[DP_Male,a] := _ForceInfect[DP_Male,a] + force_infect;

              if(t1=DP_TimeStep) then
                HA.SetForceInfect(DP_Male,a,b,t, (1- exp(-dt*force_infect)));
            end;
          end;
        end;
      end;

      (* Transfer new infections to population *)
      for a := _minDebut to DP_MaxSingleAges do
      begin
        a5 := a Div 5 + 1;
        for s := DP_Male to DP_Female do
        begin
          _NewInfAge[s,a] := _HIVNeg[s,a] * (1.0 - exp(-dt * _ForceInfect[s,a]));
          total_newinf := total_newinf + _NewInfAge[s,a];

          DP.SetNewInfections(t, a5, s, DP.GetNewInfections(t,a5,s) + _NewInfAge[s,a]);
          DP.SetNewInfectionsBySingleAge(t, s, a+1, DP.GetNewInfectionsBySingleAge(t, s, a+1) + _NewInfAge[s,a]);
          for c := DP_CD4_GT500 to DP_CD4_LT50 do
          pop1[t,s,a,c,DP_D_HIVPos] := pop1[t,s,a,c,DP_D_HIVPos]
          + _NewInfAge[s,a]
          * DP.GetAdultDistNewInfectionsCD4(DP_Data, s, age10[a], c) / 100;
          pop1[t,s,a,DP_H_No,DP_D_HIVNeg] := pop1[t,s,a,DP_H_No,DP_D_HIVNeg] - _NewInfAge[s,a];
        end;
      end;

//    { BGN DEBUGGING CODE, used to test transmission calculator. We don't want this code in beta }
//     for s := GB_Male to GB_Female do
//     begin
//        tally[s] := 0.0;
//        for a := 15 to DP_MaxSingleAges do
//          tally[s] := tally[s] + _HIVNeg[s,a] * _ForceInfect[s,a]; // _NewInfAge[s,a];
//     end;
//     tally[GB_Male] := tally[GB_Male] * 1.0;
//     { END DEBUGGING CODE }

    end;{CalcInfectionsByAge}

    procedure CalcInitGoalsPopVars(t : integer);
    var
      t1,s,a,c,d  : byte;
    begin
      for t1 := 1 to DP_TimeStep do
      begin
        for s := DP_Male to DP_Female do
          for a := DP_AllAges to DP_MAX_AGE do
          begin
            for c := DP_CD4_GT500 to DP_CD4_LT50 do
            begin
              DP.SetAIDSCd4MortNoART(t, t1, a, s, c, 0);
              DP.SetAIDSCd4MortART(t, t1, a, s, c, 0);
              DP.SetAIDSCd4Trans(t, t1, a, s, c, 0);
              for d := DP_D_HIVPos to DP_D_ARTgt12m do
                DP.SetAIDSAging15(t,s,c,d,0);
            end;
            DP.SetAIDSAging15(t,s,DP_H_No,DP_D_HIVNeg,0);
          end;
      end;

      for s := DP_Male to DP_Female do
        for a := DP_A15_19 to DP_A50_54 do
          for c := DP_Pop_H_NoART to DP_Pop_H_OnART do
            DP.SetPopSexAgeHIV(t, s, a, c, 0);

      for s := DP_Male to DP_Female do
        DP.SetPop15to49(t, s, 0);

      DP.SetPop15to49ART(t, 0);

      for s := DP_Male to DP_Female do
        DP.SetGoalsNewHIVbySex(s,0);

       (* Sum Array to get all durations and status *)
      if t > GB.GetCalcYearIdx(p) + 1 then                             //Goals issue here
        SumPopulationByHIVStatus(t-1);
    end;

    procedure CalcGoalsReceiveART(t : integer);
    var
      s,a,c,d     : byte;
      temp_total  : double;
    begin
      temp_total := 0;
      for a := 15 to 49 do
        for s := DP_Male to DP_Female do
          for c := DP_CD4_GT500 to DP_CD4_LT50 do
            for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
                temp_total := temp_total + pop1[t,s,a,c,d];

      DP.SetPop15to49ART(t, temp_total);
    end;

    procedure CalcGoalsPopVars(t : integer);
    var
      s,a,c       : byte;
      temp_total  : double;
    begin
      for s := DP_Male to DP_Female do
      begin
        //population 15-49
        temp_total := 0;
        for a := 15 to 49 do
          temp_total := temp_total + pop1[t-1,s,a,DP_H_All,DP_D_All];
        DP.SetPop15to49(t, s, temp_total);

        for a := 15 to 54 do
        begin
          a5 := trunc(a / 5) + 1;
          temp_total := 0;
          // HIV+, not on ART
          for c := DP_CD4_GT500 to DP_CD4_LT50 do
            temp_total := temp_total + pop1[t-1,s,a,c,DP_D_HIVPos];
          DP.SetPopSexAgeHIV(t, s, a5, DP_Pop_H_NoART, DP.GetPopSexAgeHIV(t, s, a5, DP_Pop_H_NoART)+temp_total);

          temp_total := 0;
          // HIV+, on ART
          for c := DP_CD4_GT500 to DP_CD4_LT50 do
            temp_total := temp_total + pop1[t-1,s,a,c,DP_D_All] - pop1[t-1,s,a,c,DP_D_HIVPos];
          DP.SetPopSexAgeHIV(t, s, a5, DP_Pop_H_OnART, DP.GetPopSexAgeHIV(t, s, a5, DP_Pop_H_OnART)+temp_total);
        end;{a loop}
      end;{s loop}
    end;

    procedure CalcGoalsNewEntrants(t: integer);
    var
      s,c,d       : byte;
      temp_age    : double;
    begin
      for s := DP_Male to DP_Female do
      begin
        //age 15, HIV-
        temp_age := pop1[t,s,15,DP_H_No,DP_D_All];
        DP.SetAIDSAging15(t,s,DP_H_No,DP_D_HIVNeg,temp_age);

        for c := DP_CD4_GT500 to DP_CD4_LT50 do
        begin
          temp_age := pop1[t,s,15,c,DP_D_HIVPos];
          DP.SetAIDSAging15(t,s,c,DP_D_HIVPos,temp_age);

          for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
          begin
            temp_age := pop1[t,s,15,c,d];
            DP.setAIDSAging15(t,s,c,d,temp_age);
          end;
        end;
      end;{s loop}
    end;

    procedure CalcGoalsAgingExits(t : integer);
    var
      s,c,d     : byte;
      temp_age  : double;
    begin
      for s := DP_Male to DP_Female do
      begin
        temp_age := pop1[t,s,50,DP_H_No,DP_D_HIVNeg];
        DP.SetAIDSAging50(t, s, DP_H_No, DP_D_HIVNeg, temp_age);

        for c := DP_CD4_GT500 to DP_CD4_LT50 do
        begin
          //no ART
          temp_age := pop1[t,s,50,c,DP_D_HIVPos];
          DP.SetAIDSAging50(t, s, c, DP_D_HIVPos, temp_age);

          //ART
          for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
          begin
            temp_age := pop1[t,s,50,c,d];
            DP.SetAIDSAging50(t, s, c, d, temp_age);
          end;
        end;{c loop}
      end;{s loop}
    end;

    procedure CalcNewHIVGoals(t1 : integer);
    var
      s, a5, a, c  : byte;
    begin
      { Calculate new adult infections }

      (* Sum Array to get all durations and status *)
      SumPopulationByHIVStatus(t);

      { Set final age }
      LastAge := 49;

//DP.SetDistOfHIV(t,DP_Male,DP_A15_19,0.299);
//DP.SetDistOfHIV(t,DP_Male,DP_A20_24,0.781);
//DP.SetDistOfHIV(t,DP_Male,DP_A30_34,0.965);
//DP.SetDistOfHIV(t,DP_Male,DP_A35_39,0.856);
//DP.SetDistOfHIV(t,DP_Male,DP_A40_44,0.691);
//DP.SetDistOfHIV(t,DP_Male,DP_A45_49,0.588);
//DP.SetDistOfHIV(t,DP_Male,DP_A50_54,0.475);
//DP.SetDistOfHIV(t,DP_Male,DP_A55_59,0.378);
//DP.SetDistOfHIV(t,DP_Male,DP_A60_64,0.285);
//DP.SetDistOfHIV(t,DP_Male,DP_A65_69,0.206);
//DP.SetDistOfHIV(t,DP_Male,DP_A70_74,0.158);
//DP.SetDistOfHIV(t,DP_Male,DP_A75_79,0.139);
//DP.SetDistOfHIV(t,DP_Male,DP_A80_up,0.0);
//
//DP.SetDistOfHIV(t,DP_Female,DP_A15_19,1.192);
//DP.SetDistOfHIV(t,DP_Female,DP_A20_24,1.166);
//DP.SetDistOfHIV(t,DP_Female,DP_A30_34,0.840);
//DP.SetDistOfHIV(t,DP_Female,DP_A35_39,0.697);
//DP.SetDistOfHIV(t,DP_Female,DP_A40_44,0.590);
//DP.SetDistOfHIV(t,DP_Female,DP_A45_49,0.505);
//DP.SetDistOfHIV(t,DP_Female,DP_A50_54,0.436);
//DP.SetDistOfHIV(t,DP_Female,DP_A55_59,0.380);
//DP.SetDistOfHIV(t,DP_Female,DP_A60_64,0.334);
//DP.SetDistOfHIV(t,DP_Female,DP_A65_69,0.296);
//DP.SetDistOfHIV(t,DP_Female,DP_A70_74,0.263);
//DP.SetDistOfHIV(t,DP_Female,DP_A75_79,0.236);
//DP.SetDistOfHIV(t,DP_Female,DP_A80_up,0.0);

        {Split incidence age pattern into single ages }
        for s := DP_Male to DP_Female do
        begin
          for a := 0 to DP_MaxSingleAges do
            sa[a] := 0;
          for a := DP_A0_4 to DP_MAX_AGE do
              p5[a] := DP.GetDistOfHIV(t,s,a);

          Calc_Single_Ages(GB,DP,p5, sa, s);
          for a := 0 to DP_MaxSingleAges do
          begin
            PrevRatio[a,s] := sa[a];
            if PrevRatio[a,s] < 0 then
            begin
              PrevRatio[a,s] := 0;
              if (a>14) and (a<60) then
                DP_AddLogEntry(GB, DP_DPAIM,'DPNewHIV','~2-10','PrevRatio[a,s]', DP_LTZ, PrevRatio[a,s],
                               'a = ' + IntToStr(a) + ',s = ' + IntToStr(s));
            end;
          end;
        end;{s loop}

        NewHIV[DP_Male]:= DP.GetGoalsNewHIVbySex(DP_Male);
        NewHIV[DP_Female]:= DP.GetGoalsNewHIVbySex(DP_Female);

        NewHIV[DP_BothSexes] := NewHIV[DP_Female] + NewHIV[DP_Male];
        NewHIV[DP_BothSexes] := NewHIV[DP_Male] + NewHIV[DP_Female];
        { Distribute new HIV infections by age and sex }

        for s := DP_Male to DP_Female do
        begin
          if(t1=1) then
          begin
            for a := DP_A0_4 to DP_A80_up do
              DP.SetNewInfections(t,a,s,0);
            for a := DP_A0 to DP_A80 do
              DP.SetNewInfectionsBySingleAge(t,s,a,0);
          end;

          for a := 0 to 14 do
            NewHIVbyAge[a] := 0;

          (*: added start*)

          { Calculate new infections if incidence were 1 in 25-29 age group }
          sum := 0;
          for a := 15 to LastAge do
            sum := sum + pop1[t,s,a,DP_H_No,DP_D_HIVNeg] * PrevRatio[a,s];
          { Scale new infections to match total required }
          if sum >0 then AdjFactor := NewHIV[s] / sum else AdjFactor := 1;

          for a := 15 to DP_MaxSingleAges do
            NewHIVbyAge[a] := AdjFactor * PrevRatio[a,s] * pop1[t,s,a,DP_H_No,DP_D_HIVNeg];
           (*: added end*)

          for a := 15 to DP_MaxSingleAges do
          begin
            a10 := age10[a];
           //: NewHIVbyAge[a] := NewHIV[s] * PrevRatio[a,s];
            a5 := trunc(a/5) + 1;
            DP.SetNewInfections(t, a5, s, DP.GetNewInfections(t,a5,s) + NewHIVbyAge[a]);
            DP.SetNewInfectionsBySingleAge(t, s, a+1, DP.GetNewInfectionsBySingleAge(t, s, a+1) + NewHIVbyAge[a]);
            (* :
            pop1[a,t,s,DP_CD4_GT500,1] := pop1[a,t,s,DP_CD4_GT500,1]
                                           + (1 - DP.GetProp350_500(DP_Data,s,a10))* NewHIVbyAge[a];
            pop1[a,t,s,DP_CD4_350_500,1] := pop1[a,t,s,DP_CD4_350_500,1]
                                           + (DP.GetProp350_500(DP_Data,s,a10))* NewHIVbyAge[a];
            *)
            (* begin*)
            for c := DP_CD4_GT500 to DP_CD4_LT50 do
            pop1[t,s,a,c,DP_D_HIVPos] := pop1[t,s,a,c,DP_D_HIVPos]
                               + NewHIVbyAge[a]
                               * DP.GetAdultDistNewInfectionsCD4(DP_Data,s,a10,c) / 100;
            (* end*)
            pop1[t,s,a,DP_H_No,DP_D_HIVNeg] := pop1[t,s,a,DP_H_No,DP_D_HIVNeg] - NewHIVbyAge[a];
          end;{a loop}
        end; {s loop}
    end;

    procedure AdjustMortForViralSuppression(p, t : integer);
    var
      HA            : THAData;
      a, c, s       : byte;
      scale, v0, vt : double;
    const
      vmax = 95.0; // cap on viral suppression levels
      smax = 0.5;  // maximum reduction in ART mortality
    begin
      HA := GetHAData(p);

      { Goals ASM models decreases in mortality on ART as viral suppression improves }
      for a := DP_CD4_15_24 to DP_CD4_45_54 do
        for s := GB_Male to GB_Female do
        begin
          scale := 1.0;
          v0 := HA.GetARTViralSuppression(s, a, 1);
          vt := HA.GetARTViralSuppression(s, a, t);
          if v0 < vmax then
          begin
            scale := min(1.0, 1.0 - smax * (vt - v0) / (vmax - v0));
            for c := DP_CD4_GT500 to DP_CD4_LT50 do
            begin
              alpha[a,c,s,DP_D_ARTlt6m  ] := scale * DP.GetAdultMortByCD4WithART0to6( DP_Data, s, a, c);
              alpha[a,c,s,DP_D_ART6to12m] := scale * DP.GetAdultMortByCD4WithART7to12(DP_Data, s, a, c);
              alpha[a,c,s,DP_D_ARTgt12m ] := scale * DP.GetAdultMortByCD4WithARTgt12( DP_Data, s, a, c);
            end;
          end;
        end;
    end;

    function EligiblePregnantWomen(a, t1 : integer): double;
    var
      c1             : byte;
      PW, sumabove, prev,
      sumbelow, Need : double;
    begin
      if pop1[t1,DP_Female,a,DP_H_All, DP_D_All] > 0 then
        prev := (pop1[t1,DP_Female,a,DP_H_All, DP_D_All]
                 - pop1[t1,DP_Female,a,DP_H_No, DP_D_HIVNeg])
                 / pop1[t1,DP_Female,a,DP_H_All, DP_D_All]
      else
        prev := 0;

      a5 := trunc(a / 5) + 1;
      { Add HIV+ Pregnant women above eligibility threshhold }
      PW := 0;
      for c1 := DP_CD4_GT500 to DP_CD4_LT50 do
        PW := PW + pop1[t1, DP_Female, a, c1, DP_D_HIVPos];
      if (DP.GetHIVTFR(DP_Data,a5,t1)>0) and (DP.GetASFR(t1, DP_AllAges)>0)
      then
        PW := PW * DP.GetTFR(t1) *
                DP.GetHIVTFR(DP_Data,a5,t1) /
               (DP.GetHIVTFR(DP_Data,a5,t1) * prev + (1.0 - prev)) *
               (DP.GetASFR(t1, a5) /
                DP.GetASFR(t1, DP_AllAges)) / 5;
      { Calculate proportion of PW above eligibility threshold }
      sumAbove := 0;
      sumBelow := 0;
      for c1 := DP_CD4_GT500 to DP_CD4_LT50 do
        if DP_CD4_lower_limit[c1] < DP.GetCD4ThreshHoldAdults(t1)
          then SumBelow := SumBelow + pop1[t1, DP_Female, a, c1, DP_D_HIVPos]
          else SumAbove := SumAbove + pop1[t1, DP_Female, a, c1, DP_D_HIVPos];
      if (SumAbove + SumBelow) > 0 then
        Need := PW * SumAbove / (SumAbove + SumBelow)
      else
        Need := 0;
      Result := Need;
    end;


    procedure CalcNeedForART(pop1 : DP_TPopType);
    var
      s, a5, a, d, c, r      : byte;
      Need, Need15, sumAbove : double;
    begin
      for s := DP_Male to DP_Female do
        for c := DP_CD4_GT500 to DP_CD4_LT50 do
          EligibleSpecialPop[c,s] := 0;

      { Put single age data on need for Tx into 5-year age arrays }
      for s := DP_Male to DP_Female do
      begin
        Need15 := 0;
        for a5 := DP_A0_4 to DP_A80_up do
          DP.SetNeedART(t, s, a5, 0);
        for a := 15 to 80 do
        begin
          a5 := trunc(a / 5) + 1;
          Need := 0;
          for c := DP_CD4_GT500 to DP_CD4_LT50 do
          begin
            { Add those on ART }
            for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
              Need := Need + pop1[t, s, a, c, d];
            { Add those not on ART but eligible by CD4 count }
            if DP_CD4_lower_limit[c] < CD4_Threshold then {DP.GetCD4ThreshHoldAdults(t) then}
              Need := Need + pop1[t, s, a, c, DP_D_HIVPos]
            else {Add any special populations }
            begin
              { TH/HIV, SDC, SW, MSM, IDU, Other }
              for r := DP_EligTreatTB_HIV to DP_EligTreatOtherPop do
                if DP.GetPopsEligTreat(r).Eligible
                   and (DP.GetPopsEligTreat(r).Year <= (DP.GetFirstYear + t - 1)) then
                begin
                  if not(((r = DP_EligTreatMSM) and (s = DP_Female)) or ((r = DP_EligTreatSexWorkers) and (s = DP_male))) then
                  begin
                    Need := Need + DP.GetPopsEligTreat(r).PercentHIV / 100
                                 * pop1[t, s, a, c, DP_D_HIVPos];
                    EligibleSpecialPop[c,s] := EligibleSpecialPop[c,s]
                               + DP.GetPopsEligTreat(r).PercentHIV / 100
                                  * pop1[t, s, a, c, DP_D_HIVPos];
                  end;
                end;
            end;
          end; {c loop}
          { HIV+ Pregnant women }

          if (s = DP_Female) and (a>=15) and (a<=49) and
             (DP.GetPopsEligTreat(DP_EligTreatPregnantWomen).Eligible)
             and (DP.GetPopsEligTreat(DP_EligTreatPregnantWomen).Year
                  <= (DP.GetFirstYear + t - 1)) then
          begin
            Need := Need + EligiblePregnantWomen(a,t);
            { Distribute eligible pregnant women to CD4 categories above
              eligibility threshold }
            SumAbove := 0;
            for c := DP_CD4_GT500 to DP_CD4_LT50 do
              if DP_CD4_lower_limit[c] >= CD4_Threshold then {DP.GetCD4ThreshHoldAdults(t) then}
                SumAbove := SumAbove + pop1[t, s, a, c, DP_D_HIVPos];
            if SumAbove > 0 then
              for c := DP_CD4_GT500 to DP_CD4_LT50 do
                if DP_CD4_lower_limit[c] >= CD4_Threshold then {DP.GetCD4ThreshHoldAdults(t) then}
                  EligibleSpecialPop[c,s] := EligibleSpecialPop[c,s]
                                          + EligiblePregnantWomen(a,t)
                                          * pop1[t,s,a,c,DP_D_HIVPos] / SumAbove;
          end;
          DP.SetNeedART(t, s, a5, DP.GetNeedART(t, s, a5) + Need);
          if a >= 15 then
            Need15 := Need15 + Need;
        end; {a loop}
        DP.SetARTCalc_15Plus(DP_TotNeedART, s, t, Need15);
      end; {s loop}


      { Sum up for all ages }
      for s := DP_Male to DP_Female do
      begin
        sum5 := 0;
        for a5 := DP_A0_4 to DP_A80_up do
        begin
          sum5 := sum5 + DP.GetNeedART(t, s, a5);
          DP.SetNeedART(t, s, DP_AllAges, sum5);
        end;
      end; {s loop}

      { Sum up for both sexes }
      for a5 := DP_AllAges to DP_A80_up do
        DP.SetNeedART(t, DP_BothSexes, a5, DP.GetNeedART(t, DP_Male, a5)
            + DP.GetNeedART(t, DP_Female, a5));
      DP.SetARTCalc_15Plus(DP_TotNeedART, DP_BothSexes, t,
           DP.GetARTCalc_15Plus(DP_TotNeedART, DP_Male, t)
           + DP.GetARTCalc_15Plus(DP_TotNeedART, DP_Female, t));
    end; { CalcNeedForART }

    procedure CalcNeedForARTEndYear;
    var
      s, a5, a, d, c, r : byte;
      Need              : double;
    begin

      { Put single age data on need for Tx into 5-year age arrays }
      for s := DP_Male to DP_Female do
      begin
        for a5 := DP_A15_19 to DP_A80_up do
          DP.SetNeedARTDec31(t, s, a5, 0);
        for a := 15 to 80 do
        begin
          a5 := trunc(a / 5) + 1;
          Need := 0;
          for c := DP_CD4_GT500 to DP_CD4_LT50 do
          begin
            { Add those on ART }
            for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
              Need := Need + pop1[t, s, a, c, d];
            { Add those not on ART but eligible by CD4 count }
            if DP_CD4_lower_limit[c] < CD4_Threshold then {DP.GetCD4ThreshHoldAdults(t) then}
              Need := Need + pop1[t, s, a, c, DP_D_HIVPos]
            else {Add any special populations }
            begin
              { TH/HIV, SDC, SW, MSM, IDU, Other }
              for r := DP_EligTreatTB_HIV to DP_EligTreatOtherPop do
                if DP.GetPopsEligTreat(r).Eligible
                  and (DP.GetPopsEligTreat(r).Year <= (DP.GetFirstYear + t - 1)) then

                  if not(((r = DP_EligTreatMSM) and (s = DP_Female)) or ((r = DP_EligTreatSexWorkers) and (s = DP_male))) then
                    begin
                      Need := Need +
                               DP.GetPopsEligTreat(r).PercentHIV / 100
                               * pop1[t, s, a, c, DP_D_HIVPos];
                    end;
            end;
          end; {c loop}
          { HIV+ Pregnant women }
          if (s = DP_Female) and (a>=15) and (a<=49) and
             (DP.GetPopsEligTreat(DP_EligTreatPregnantWomen).Eligible)
             and (DP.GetPopsEligTreat(DP_EligTreatPregnantWomen).Year
                  <= (DP.GetFirstYear + t - 1)) then
               Need := Need + EligiblePregnantWomen(a,t);
          DP.SetNeedARTDec31(t, s, a5, DP.GetNeedARTDec31(t, s, a5) + Need);
        end; {a loop}
      end; {s loop}

      { Sum up for all ages }
      for s := DP_Male to DP_Female do
      begin
        sum5 := 0;
        for a5 := DP_A15_19 to DP_A80_up do
        begin
          sum5 := sum5 + DP.GetNeedARTDec31(t, s, a5);
          DP.SetNeedARTDec31(t, s, DP_AllAges, sum5);
        end;
      end; {s loop}
      DP.SetNeedARTDec31(t, DP_BothSexes, DP_AllAges,
                        DP.GetNeedARTDec31(t, DP_Male, DP_AllAges)
                        + DP.GetNeedARTDec31(t, DP_Female, DP_AllAges));

      { Sum up for both sexes }
      for a5 := DP_A15_19 to DP_A80_up do
        DP.SetNeedARTDec31(t, DP_BothSexes, a5, DP.GetNeedARTDec31(t, DP_Male, a5)
            + DP.GetNeedARTDec31(t, DP_Female, a5));
    end; { CalcNeedForARTEndYear }

    procedure CalcNewARTPatients(s : byte);
    var
      c, c2, a, a1      : byte;
      sum, sum1,
      propByMort,
      propByElig        : double;
      MortRate,
      Prop1, Prop2      : array [DP_CD4_GT500 .. DP_CD4_LT50] of double;
      EligByAge         : array [DP_CD4_15_24 .. DP_CD4_45_54, DP_CD4_GT500 .. DP_CD4_LT50, DP_Male .. DP_Female] of double;
    begin
      { Calculate eligible adults by CD4 category }
      for c := DP_CD4_GT500 to DP_CD4_LT50 do
        Eadults[c,s] := 0;
      for c := DP_CD4_GT500 to DP_CD4_LT50 do
      begin
        if DP_CD4_lower_limit[c] < CD4_Threshold then {DP.GetCD4ThreshHoldAdults(t) then}
          for a := 15 to 80 do
            Eadults[c,s] := Eadults[c,s] + pop1[t, s, a, c, DP_D_HIVPos]
        else
          Eadults[c,s] := Eadults[c,s] + EligibleSpecialPop[c,s];
      end;

      if DP.GetNewARTPatAllocMethod = DP_AdvOpt_ART_LowCD4First then
      begin
        { Allocate ART to lowest CD4 cell counts first }
        sum := NewART;
        for c := DP_CD4_LT50 downto DP_CD4_GT500 do
        begin
          if Eadults[c,s] < sum then
          begin
            NewPatients[c,s] := Eadults[c,s];
            if NewART > 0 then DP.SetNewARTCD4(t,c, Eadults[c,s] / NewART);
            sum := sum - Eadults[c,s];
          end
          else
          begin
            NewPatients[c,s] := sum;
            if NewART > 0 then DP.SetNewARTCD4(t, c, sum / NewART);
            sum := 0;
          end;
        end;
      end
      else
      begin

        { If values are given for median CD4 count at ART initiation then adjust to match }
//        if DP.GetMedCD4CountInit(t) > 0 then
//        begin
//          { Sum total eligible adults }
//          sum := 0;
//          for c := DP_CD4_LT50 downto DP_CD4_GT500 do
//            sum := sum + Eadults[c,s];
//          { Calculate proportions }
//          if sum > 0 then
//            for c := DP_CD4_LT50 downto DP_CD4_GT500 do
//              Prop1[c] := Eadults[c,s] / sum;
//
//          { First check for median above 500 }
//          if DP.GetMedCD4CountInit(t) > DP_CD4_lower_limit[DP_CD4_GT500] then
//          begin
//            c1 := DP_CD4_GT500;
//            { Divide assuming a mythical top for the CD4 > 500 compartment of 1000.
//              Adjust the proportion in the lower part of the 500+ compartment
//              to only contain the fraction below the CD4 median}
//            Prop1[DP_CD4_GT500] := Prop1[DP_CD4_GT500]
//                                   * (DP.GetMedCD4CountInit(t) - DP_CD4_lower_limit[DP_CD4_GT500])
//                                     / (1000 - DP_CD4_lower_limit[DP_CD4_GT500]);
//           { Sum proportions }
//           sum1 := 0;
//           for c := DP_CD4_LT50 downto DP_CD4_GT500 do
//             sum1 := sum1 + Prop1[c];
//           { Normalize to 0.5 }
//           for c := DP_CD4_LT50 downto DP_CD4_GT500 do
//             Prop1[c] := Prop1[c] / sum1 * 0.5;
//           { Now add 50% to the top category to deal with the fact that is has the additional 50% }
//           Prop1[DP_CD4_GT500] := Prop1[DP_CD4_GT500] + 0.5;
//          end
//          else
//          begin
//            { Find CD4 category containing the median }
//            c1 := DP_CD4_LT50;
//            while (DP_CD4_lower_limit[c1-1] < DP.GetMedCD4CountInit(t)) do
//              c1 := c1 - 1;
//            { Divide into two halves, adjust each half to contain 50% }
//            { Sum proportions for categories below the median }
//            sum1 := 0;
//            for c := DP_CD4_LT50 downto c1+1 do
//              sum1 := sum1 + Prop1[c];
//            sum1 := sum1 + Prop1[c1]
//                         * (DP.GetMedCD4CountInit(t) - DP_CD4_lower_limit[c1])
//                         / (DP_CD4_lower_limit[c1-1] - DP_CD4_lower_limit[c1]);
//            { Normalize to 50% }
//            if sum1 > 0 then
//              for c := DP_CD4_LT50 downto c1+1 do
//                Prop1[c] := Prop1[c] / sum1 * 0.5;
//
//            sum1 := 0;
//            for c := DP_CD4_LT50 downto c1+1 do
//              sum1 := sum1 + Prop1[c];
//            { Sum2 is the leftover from the lower half }
//            sum2 := 0.5 - sum1;
//
//            { Now allocate to categories above the median }
//            sum3 := 0;
//            for c := c1-1 downto DP_CD4_GT500 do
//              sum3 := sum3 + Prop1[c];
//            sum3 := sum3 + Prop1[c1] * (DP_CD4_lower_limit[c1-1] - DP.GetMedCD4CountInit(t))
//                                        / (DP_CD4_lower_limit[c1-1] - DP_CD4_lower_limit[c1]);
//            if sum3 > 0 then
//              for c := c1-1 downto DP_CD4_GT500 do
//                Prop1[c] := Prop1[c] / sum3 * 0.5;
//            sum3 := 0;
//            for c := c1-1 downto DP_CD4_GT500 do
//              sum3 := sum3 + Prop1[c];
//            { Sum3 becomes the leftover from the upper half }
//            sum3 := 0.5 - sum3;
//            Prop1[c1] := sum3 + sum2;
//          end;
//        end { Median CD4 count entered }
//        else
//        begin
        propByMort := DP.GetNewARTPatAlloc(DP_AdvOpt_ART_ExpMort);
        propByElig := DP.GetNewARTPatAlloc(DP_AdvOpt_ART_PropElig);

        { Calculate Option 1 : allocate according to number eligible }
        { Same proportion starting treatment in each CD4 category }
        sum := 0;
        for c := DP_CD4_GT500 to DP_CD4_LT50 do
          if Eadults[c,s] > 0 then sum := sum + Eadults[c,s];
        for c := DP_CD4_GT500 to DP_CD4_LT50 do
          if Eadults[c,s] <= 0 then Prop1[c] := 0
          else Prop1[c] := Eadults[c,s] / sum;

        { Calculate Option 2 : allocation is weighted by mortality }
        { First calculate average non-ART mortality for CD category }
        for c := DP_CD4_GT500 to DP_CD4_LT50 do
          for a1 := DP_CD4_15_24 to DP_CD4_45_54 do
            EligByAge[a1, c, s] := 0;

        for c := DP_CD4_GT500 to DP_CD4_LT50 do
        begin
          for a := 15 to 80 do
          begin
            a1 := age10[a];
            EligByAge[a1, c, s] := EligByAge[a1, c, s] + pop1[t, s, a, c, DP_D_HIVPos];
          end; {a loop}
          sum := 0;
          sum1 := 0;
          for a1 := DP_CD4_15_24 to DP_CD4_45_54 do
          begin
            case a1 of
               DP_CD4_15_24:
                  sum := sum + max(mu[a1, c, s],valueMinMort[DP_A20,s]) * EligByAge[a1, c, s];
               DP_CD4_25_34:
                  sum := sum + max(mu[a1, c, s],valueMinMort[DP_A30,s]) * EligByAge[a1, c, s];
               DP_CD4_35_44:
                  sum := sum + max(mu[a1, c, s],valueMinMort[DP_A40,s]) * EligByAge[a1, c, s];
               DP_CD4_45_54:
                  sum := sum + max(mu[a1, c, s],valueMinMort[DP_A50,s]) * EligByAge[a1, c, s];
            end;
            sum1 := sum1 + EligByAge[a1, c, s];
          end;
          if sum1 = 0 then
            MortRate[c] := 0
          else
            MortRate[c] := sum / sum1;
        end; { c loop }

        RemainingNewART := NewART;
        { Second, calculate weighted average of number eligible and mortality }
        for c := DP_CD4_LT50 downto DP_CD4_GT500 do
        begin
          sum := 0;
          for c2 := c downto DP_CD4_GT500 do
            sum := sum + Eadults[c2,s] * MortRate[c2];
          if sum <= 0 then NewPatients[c,s] := 0
          else NewPatients[c,s] := RemainingNewART * Eadults[c,s] * MortRate[c] / sum;
          if NewPatients[c,s] > Eadults[c,s] then
            NewPatients[c,s] := Eadults[c,s];
          { Adjust NewART to be new ART patients not yet allocated }
          RemainingNewART := RemainingNewART - NewPatients[c,s];
        end; { c loop }

  //      sum := 0;
  //      for c := DP_CD4_LT50 downto DP_CD4_GT500 do
  //        if Eadults[c,s]>0 then sum := sum + Eadults[c,s] * MortRate[c];
  //      for c := DP_CD4_LT50 downto DP_CD4_GT500 do
  //        if Eadults[c,s]>0 then Prop2[c] := Eadults[c,s] * MortRate[c] / sum else Prop2[c] := 0;

        { Third calculate proportion starting }
        for c := DP_CD4_GT500 to DP_CD4_LT50 do
          if Eadults[c,s] <= 0 then Prop2[c] := 0
          else Prop2[c] := NewPatients[c,s] / NewART;

        { Finally, average methods one and two }
        for c := DP_CD4_GT500 to DP_CD4_LT50 do
          Prop1[c] := propByElig * Prop1[c] + propByMort * Prop2[c];
//        end; { Median CD4 count not entered }

        for c := DP_CD4_GT500 to DP_CD4_LT50 do
        begin
          NewPatients[c,s] := Prop1[c] * NewART;
          DP.SetNewARTCD4(t, c, Prop1[c]);
        end;
      end;

    end; { CalcNewARTPatients }

    { Populates NewChildART[a,s,c,d] for ages a0 <= a <= a1. Requires that
      CalcNumChildrenEligTreat has been called to populate elig[a,s,c,d].
      Declared local to Do_AIDS_Calculations so that palpha and NewChildART
      are in scope }
    procedure CalcNewChildARTPatients(t, ARTgroup : integer);
    var
      adj, retained, initByAge, denom,
      NewFLART, ARTThisYear, ARTLastYear,
      onART, totalNeed, unmetNeed, ARTdeaths : double;
      eligible : array[DP_Male..DP_Female, 0..14, DP_CD4_Per_GT30..DP_CD4_MaxCategories] of double;
      s, a, c, d, a5,
      ageMin, ageMax : integer;
    begin
      { Default: assume ARTgroup = DP_PerChildHIVRecART }
      ageMin := 0;
      ageMax := 14;
      a5 := DP_MAX_AGE; // a5 should only be used if ART inputs by 5-year age group
      case ARTgroup of
        DP_PerChildHIVRecART      : begin ageMin :=  0; ageMax := 14; end;
        DP_PerChildHIVRecART0_4   : begin ageMin :=  0; ageMax :=  4; a5 := DP_A0_4;   end;
        DP_PerChildHIVRecART5_9   : begin ageMin :=  5; ageMax :=  9; a5 := DP_A5_9;   end;
        DP_PerChildHIVRecART10_14 : begin ageMin := 10; ageMax := 14; a5 := DP_A10_14; end;
      end;

      { Cache unmet need, summed across modes of transmission }
      for s := DP_Male to DP_Female do
        for a := ageMin to ageMax do
          for c := DP_CD4_Start[a] to DP_CD4_Final[a] do
          begin
            eligible[s,a,c] := 0.0;
            for d := DP_P_Perinatal to DP_P_BF12 do
              eligible[s,a,c] := eligible[s,a,c] + elig[a,s,c,d];
          end;

      onART := 0.0;
      unmetNeed := 0.0;
      ARTdeaths := 0.0;
      for s := DP_Male to DP_Female do
        for a := ageMin to ageMax do
          for c := DP_CD4_Start[a] to DP_CD4_Final[a] do
          begin
            for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
              onART := onART + pop1[t,s,a,c,d];
            ARTdeaths := ARTdeaths + pop1[t,s,a,c,DP_D_ARTlt6m ] * (palpha[a,c,s,DP_D_ARTlt6m] + palpha[a,c,s,DP_D_ART6to12m]) / 2;
            ARTdeaths := ARTdeaths + pop1[t,s,a,c,DP_D_ARTgt12m] * palpha[a,c,s,DP_D_ARTgt12m];
            unmetNeed := unmetNeed + eligible[s,a,c];
          end;
      totalNeed := unmetNeed + onART;

      {Get children on ART in age=ARTgroup last year. If program statistics are
       provided for that age group, use them directly and convert from % to # if
       needed. If ARTgroup is 0-14 and ART was specified by 5-year age group in
       year t-1, we just sum inputs across 5-year age groups. If ARTgroup is a
       5-year age group and ART was inputted for 0-14 in t-1, then we
       approximate the 5-year age group input by multiplying the 0-14 input by
       the calculated proportion of children on ART in age=ARTgroup.}
      ARTLastYear := DP.GetChildTreatInputs(ARTgroup, t - 1);
      if ARTLastYear = DPNotAvail then
      begin
        if ARTgroup = DP_PerChildHIVRecART then
          ARTLastYear := DP.GetChildTreatInputs(DP_PerChildHIVRecART0_4,   t-1) +
            DP.GetChildTreatInputs(DP_PerChildHIVRecART5_9,   t-1) +
            DP.GetChildTreatInputs(DP_PerChildHIVRecART10_14, t-1)
        else
        begin
          denom := DP.GetChildARTCalc(DP_NumOnART, GB_Male, t-1) +
            DP.GetChildARTCalc(DP_NumOnART, GB_Female, t-1);
          if denom > 0 then
            ARTLastYear := DP.GetChildTreatInputs(DP_PerChildHIVRecART, t-1) *
              DP.GetNumChildOnARTByAge(t-1, GB_BothSexes, a5) / denom
          else
            ARTLastYear := 0.0;
        end;
      end;
//      if ARTLastYear = DPNotAvail then
//      begin
//        if ARTgroup = DP_PerChildHIVRecART
//        then ARTLastYear := DP.GetChildTreatInputs(DP_PerChildHIVRecART0_4,   t-1)
//                          + DP.GetChildTreatInputs(DP_PerChildHIVRecART5_9,   t-1)
//                          + DP.GetChildTreatInputs(DP_PerChildHIVRecART10_14, t-1)
//        else ARTLastYear := DP.GetChildTreatInputs(DP_PerChildHIVRecART, t-1)
//                          * DP.GetNumChildOnARTByAge(t-1, GB_BothSexes, a5)
//                          / (DP.GetChildARTCalc(DP_NumOnART, GB_Male, t-1) + DP.GetChildARTCalc(DP_NumOnART, GB_Female, t-1));
//      end
//      else if DP.GetChildARTByAgeGroupPerNum(ARTgroup, t - 1) = DP_Percent then
//        ARTLastYear := totalNeed * ARTLastYear / 100;

      ARTThisYear := DP.GetChildTreatInputs(ARTgroup, t);
      if DP.GetChildARTByAgeGroupPerNum(ARTgroup, t) = DP_Percent then
        ARTThisYear := totalNeed * ARTThisYear / 100;

      retained := 1 - DP.GetPercLostFollowUpChild(t) / 100;
      NewFLART := 0.5 * (ARTLastYear + ARTThisYear) - onART * retained;

      { Increase number starting ART to account for those who will die in the first year }
      NewFLART := NewFLART + ARTdeaths;
      if NewFLART < 0 then NewFLART := 0;

      { Prevent the number initiating ART from exceeding unmet need }
      if NewFLART > (unmetNeed + OnART * DP.GetPercLostFollowUpChild(t)/100)
        then NewFLART := unmetNeed + OnART * DP.GetPercLostFollowUpChild(t)/100;

      { Distribute according to IeDEA data }
      { DP.GetChildARTDist is the probability of starting ART among HIV+ not on ART}
      { First calculate using standard probabilities }
      { Need to use a+1 because the variable is defined with DP_A0 t0 DP_A14 }
      initByAge := 0.0;
      for s := DP_Male to DP_Female do
        for a := ageMin to ageMax do
          for c := DP_CD4_Start[a] to DP_CD4_Final[a] do
            initByAge := initByAge + DP.GetChildARTDist(DP_Data, a+1, t) * eligible[s,a,c];

      { Calculate adjustment factor to match the desired number of new ART patients }
      if initByAge = 0 then adj := 1 else adj := NewFLART / initByAge;

      { Repeat the calculation of those newly starting by age, adjusted to desired total }
      for s := DP_Male to DP_Female do
        for a := ageMin to ageMax do
          for c := DP_CD4_Start[a] to DP_CD4_Final[a] do
            for d := DP_P_Perinatal to DP_P_BF12 do
              NewChildART[a,s,c,d] := min(1, adj * DP.GetChildARTDist(DP_Data, a+1, t)) * elig[a,s,c,d];

    end; { CalcNewChildARTPatients }

  begin
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
    for a := 15 to 24 do age10[a] := DP_CD4_15_24;
    for a := 25 to 34 do age10[a] := DP_CD4_25_34;
    for a := 35 to 44 do age10[a] := DP_CD4_35_44;
    for a := 45 to DP_MaxSingleAges do age10[a] := DP_CD4_45_54;

    sr := 0; sum := 0; sum4 := 0;
    sum5 := 0; sum6 := 0; sum7 := 0; sum8 := 0;
    sum9 := 0;
//    sum10 := 0;
    NHB := 0; adultHIVnum := 0;
    NewHIVnum := 0; pop5 := 0; pop25 := 0; prev25 := 0;
//    temp_cd4 := 0;
    adj := 0;
    PMTCTFraction := 0;
//    v1 := 0; v2 := 0;
    {v3 := 0;} v4 := 0;
    v5 := 0; v6 := 0; AIDSdeaths := 0; CtxPer := 0; ArtPer := 0;
    BothPer := 0; NonePer := 0; TxEff := 0; InfTxEff := 0;
    ChildTxEff := 0; PrevAdj := 0; ARTpercent := 0;
    FractionProgressing := 0; progressors := 0;
    NeedForFirstLineART := 0; NeedForSecondLineART := 0;
    NeedForFirstLineARTmale := 0;
    NeedForFirstLineARTfemale := 0;
    UnmetNeedForFirstLineART := 0;
    UnmetNeedForFirstLineARTmale := 0;
    UnmetNeedForFirstLineARTfemale := 0;
    NeedForFirstLine1To4 := 0;
//    OnFirstLineArt := 0;
    OnSecondLineArt := 0;
    NewFLArt := 0;
//    OnFirstLineArtmale := 0;
//    OnFirstLineArtfemale := 0;
    OnCTX := 0; NewCTX := 0;
    NewInfFromBFGT12 := 0;
    UnmetNeedForCTX := 0;
    FYARTsurv := 0; PercentExposed := 0; IdentifiedNeed := 0;
    IdentifiedNeedFLART := 0;
    NeededART := 0; ARTLastYear := 0;
    MonthsOfExposure := 0; ExclusiveBF := 0; MixedBF := 0;
    MaleHIVChildren := 0; FemaleHIVChildren := 0;
    Excess := 0; NotOnTx := 0;
//    TotalPopAge := 0;
//    nART := 0;
    MedianDurationBF := 0;
    ArtEff := 0;
    ChildrenOnCTX := 0; CTXnum := 0;
    ProportionAverted := 0;
    SumOnART0_4 := 0; SumOnART5_9 := 0; SumOnART10_14 := 0;
    ARTsum := 0; MedCD4 := 0;
//    value := 0;
    ART_Excess := 0;
    Number_Eligible := 0; sum_start := 0;
    sum_end := 0;
//    ARTsurvivors := 0;
    Count_eligible := 0; Percent_eligible := 0;
//    SumART_15_49 := 0;
    SumART_BS := 0;
    SumART_15_49_BS := 0;
    NewART := 0; RemainingNewART := 0;
//    StartART := 0;
    totalE := 0;
//    NotOnART := 0;
//    OnART := 0;
    FeedOpt := 0;
    (*a := 0;*) (*s := 0;*) (*d := 0;*) (*h := 0;*) (*t1 := 0;*) (*a1 := 0;*) (*a2 := 0;*) a7 := 0; a5 := 0;
    ysi := 0; count := 0; (*d1 := 0;*) ag := 0; (*a9 := 0;*)
    agroup := 0; (*r := 0;*) c1 := 0; c2 := 0; (*c := 0;*) (*s1 := 0;*)
    yr  := 0;
//    InFile  := nil;
//    Outfile := nil;
    tstring := '';
    for s := DP_Male to DP_Female do
      ARTTwoYearsAgoBySex[s] := 0;
    for a := DP_MixFeed_HIVMoth to DP_RepFeed_HIVMoth do
    begin
      FeedPercent[a] := 0;
      FeedPercentIn[a] := 0;
      FeedPercentOut[a] := 0;
    end;

    temps := ''; temps1 := ''; temps2 := ''; temps3 := '';
    for s := DP_Male to DP_Female do
    begin
      for c := DP_CD4_GT500 to DP_CD4_LT50 do
      begin
        EligibleSpecialPop[c,s] := 0;
      end;
    end;

    for c := DP_CD4_GT500 to DP_CD4_LT50 do
      Prop3[c] := 0;

    for a := 1 to 4 do
      for c := DP_CD4_GT500 to DP_CD4_LT50 do
        for s := DP_Male to DP_Female do
        begin
          lambda[a,c,s] := 0;
        end;

    for a := DP_CD4_0_2 to DP_CD4_5_14 do
      for h := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
        for s := DP_Male to DP_Female do
          for r := DP_P_Perinatal to DP_P_BF12 do
            plambda[a,h,s,r] := 0;

    for a := 1 to 4 do
      for c := DP_CD4_GT500 to DP_CD4_LT50 do
        for s := DP_Male to DP_Female do
          mu[a,c,s] := 0;
    for a := DP_CD4_0_2 to DP_CD4_5_14 do
      for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
        for s := DP_Male to DP_Male do
          for r := DP_P_Perinatal to DP_P_BF12 do
            pmu[a,c,s,r] := 0;

    for a := 1 to 4 do
      for c := DP_CD4_GT500 to DP_CD4_LT50 do
        for s := DP_Male to DP_Female do
          for h := DP_D_ARTlt6m to DP_D_ARTgt12m do
            alpha[a,c,s,h] := 0;

    for a := 1 to 4 do
      for c := DP_CD4_GT500 to DP_CD4_LT50 do
        for s := DP_Male to DP_Female do
        begin
          beta1[a,c,s] := 0;
          beta2[a,c,s] := 0;
        end;

    for a := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
      CD_perc_at_infection[a] := 0;

    for a := DP_CD4_GT500-1 to DP_CD4_LT50 do
      DP_CD4_lower_limit[a] := 0;

    for c := DP_CD4_GT500 to DP_CD4_LT50 do
      for r := DP_D_HIVPos to DP_D_Max do
      begin
        entrants[c,r] := 0;
        exits[c,r] := 0;
      end;

    for c := 0 to DP_CD4_LT50 do
      for s := DP_Male to DP_Female do
      begin
        NewPatients[c,s] := 0;
        Eadults[c,s] :=0;
      end;

    for a := 0 to 14 do
      for s := DP_BothSexes to DP_Female do
        for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
          for h := DP_D_All to DP_P_BF12 do
            NewChildART[a,s,c,h] := 0;

    for c := DP_CD4_GT500 to DP_CD4_LT50 do
      Age0to14CD4Distr[c] := 0;
    for a := 15 to DP_MaxSingleAges do
      for s := DP_Male to DP_Female do
        valueMinMort[a,s] := 0;

    for a := 15 to DP_MaxSingleAges do
      for s := DP_Male to DP_Female do
        for d := DP_NoTreat to DP_OnART do
          AdultAIDSDeaths[a,s,d] := 0.0;

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////
    for s := DP_BothSexes to DP_Female do
      SetLength(EligibleAdults[s], GB.GetNumDynYears(p));

    { Transfer data from editor arrays to equation parameters }
    for a := DP_CD4_15_24 to DP_CD4_45_54 do
    begin
      for c := DP_CD4_GT500 to DP_CD4_LT50 do
      begin
        for s := DP_Male to DP_Female do
        begin
          if c = DP_CD4_LT50 then
            lambda[a, c, s] := 0.0000001
          else
  {          lambda[a, c, s] := 1 - exp(-1 * 1 / DP.GetAdultAvgDurByCD4(DP_Data, s, a, c)
                / DP_TimeStep);}
            // lambda[a,c,s] := DP.GetAdultAnnRateProgressLowerCD4(DP_Data, s, a, c) / DP_TimeStep;
             lambda[a, c, s] := 1
               - exp(-1 * DP.GetAdultAnnRateProgressLowerCD4(DP_Data, s, a, c)
                         / DP_TimeStep);


          if lambda[a, c, s] <= 0 then
          begin
            lambda[a, c, s] := 1;
            if lambda[a, c, s] < 0 then
              DP_AddLogEntry(GB, DP_DPAIM, 'DPPROJA', '~1-20', 'lambda[a, c, s]',
                DP_LTEZ, lambda[a, c, s], 'a = ' + IntToStr(a)
                  + ',c = ' + IntToStr(c) + ',s = ' + IntToStr(s));
          end;
          mu[a, c, s] := DP.GetAdultMortByCD4NoART(DP_Data, s, a, c);
          alpha[a, c, s, DP_D_ARTlt6m] := DP.GetAdultMortByCD4WithART0to6(DP_Data, s, a, c);
          alpha[a, c, s, DP_D_ART6to12m] := DP.GetAdultMortByCD4WithART7to12(DP_Data, s, a, c);
          alpha[a, c, s, DP_D_ARTgt12m] := DP.GetAdultMortByCD4WithARTGT12(DP_Data, s, a, c);
          beta1[a, c, s] := DP.GetIncCD4FirstYrWithART(DP_Data);
          beta2[a, c, s] := DP.GetIncCD4SubYrWithART(DP_Data);
        end; {s loop}
      end; {c loop}
    end; {a loop}

  (* Transfer pediatric model parameters *)
    for c := DP_CD4_Per_GT30 to DP_CD4_Per_5_10 do
      for s := DP_Male to DP_Female do
        for a := DP_CD4_0_2 to DP_CD4_5_14 do
          for d := DP_P_Perinatal to DP_P_BF12 do
            plambda[a, c, s, d] := DP.GetChildAnnRateProgressLowerCD4(DP_Data,s,a+1,c);

    for a := DP_CD4_0_2 to DP_CD4_5_14 do
      for d := DP_P_Perinatal to DP_P_BF12 do
        for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
          for s := DP_Male to DP_Female do
            pmu[a,c,s,d] := DP.GetChildMortByCD4NoART(DP_Data,a+1,d,c);

    a1 := DP_ART_0_1;
    a2 := DP_CD4_0t4;
    for a := 0 to 14 do
    begin
      case a of
        0      : begin a1 := DP_ART_0_1;   a2 := DP_CD4_0t4;  end;
        1..2   : begin a1 := DP_ART_1_2;   a2 := DP_CD4_0t4;  end;
        3..4   : begin a1 := DP_ART_2_5;   a2 := DP_CD4_0t4;  end;
        5..9   : begin a1 := DP_ART_5_10;  a2 := DP_CD4_5t14; end;
        10..14 : begin a1 := DP_ART_10_15; a2 := DP_CD4_5t14; end;
      end;
      v1 := DP.GetChildMortalityRatesMultiplier;
      v2 := DP.GetChildMortalityRates(DP_Data, a2, DP_MortRates_LT12Mo, t);
      v3 := DP.GetChildMortalityRates(DP_Data, a2, DP_MortRates_GT12Mo, t);
      for s := DP_Male to DP_Female do
        for c := DP_CD4_Ped_Top to DP_CD4_Ped_LT200 do
        begin
          palpha[a,c,s,DP_D_ARTlt6m  ] := DP.GetChildMortByCD4WithART0to6( DP_Data,s,a1,c) * v1 * v2;
          palpha[a,c,s,DP_D_ART6to12m] := DP.GetChildMortByCD4WithART7to12(DP_Data,s,a1,c) * v1 * v2;
          palpha[a,c,s,DP_D_ARTgt12m ] := DP.GetChildMortByCD4WithARTGT12( DP_Data,s,a1,c) * v1 * v3;
        end;
    end;

    for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
      CD_perc_at_infection[c] := DP.GetChildDistNewInfectionsCD4(DP_Data, c);

    for s := DP_Male to DP_Female do
      for a := 15 to DP_MaxSingleAges do
        valueMinMort[a,s] := 0;

    DP_CD4_lower_limit[DP_CD4_GT500-1] := 1000;
    DP_CD4_lower_limit[DP_CD4_GT500] := 500;
    DP_CD4_lower_limit[DP_CD4_350_500] := 350;
    DP_CD4_lower_limit[DP_CD4_250_349] := 250;
    DP_CD4_lower_limit[DP_CD4_200_249] := 200;
    DP_CD4_lower_limit[DP_CD4_100_199] := 100;
    DP_CD4_lower_limit[DP_CD4_50_99] := 50;
    DP_CD4_lower_limit[DP_CD4_LT50] := 0;

    MedianDurationBF := 0;

    (* Calculate number initiating treatment. Input may be % or #
       For population, t refers to July 1
       For ART inputs, t refers to Dec 31 of same year
       Increase in number on ART from t-1 to t will be
       based on ART growth from t-2 to t-1 for July 1 (t-1) to Dec 31
       and from t-1 to t for Dec 31 to July 1 (t) *)

    if t > 1 then
    begin

      CD4_Threshold := DP.GetCD4ThreshHoldAdults(t);
      if CD4_Threshold = 200 then CD4_Threshold := 250;

      { Start adult calculations }

      CalcNeedForART(pop1);

      if GB.ModuleActive(p, GB_HA) then
      begin
        HA_BehaveChangeFactorCalc(p,t);
        AdjustMortForViralSuppression(p,t);
      end;

      (* if Goals is active, Goals PreInnerLoopCalcs *)
      if (GB.ModuleActive(p, GB_HV)) then
      begin

        {initialize variables used in mapping between AIM and Goals}
        CalcInitGoalsPopVars(t);

        {set AIM population sizes. Used in aging and transition rates calcs}
        CalcGoalsPopVars(t);

        {set AIM 15 year olds new entrants to Goals}
        CalcGoalsNewEntrants(t);

        {set AIM 49 year olds. Used as exit rate for Goals}
        CalcGoalsAgingExits(t);

        { Calculate reduced mortality by ART by 50% as viral suppression increases to 95% }
        { Input is percent not virally suppressed, initial default is 0.25 }
        { Fast-Track target is 0.05 }
        alpha_mult := 1;
        if DP_GetHVInfectMultiplierOnART(p,1) > 0.05 then
          alpha_mult := 1 - 0.5
                        * (DP_GetHVInfectMultiplierOnART(p,1) - DP_GetHVInfectMultiplierOnART(p,t))
                        / (DP_GetHVInfectMultiplierOnART(p,1) - 0.05);
        if alpha_mult > 1 then alpha_mult := 1;

        for a := DP_CD4_15_24 to DP_CD4_45_54 do
          for c := DP_CD4_GT500 to DP_CD4_LT50 do
            for s := DP_Male to DP_Female do
            begin
              alpha[a, c, s, DP_D_ARTlt6m]   := alpha_mult * DP.GetAdultMortByCD4WithART0to6(DP_Data, s, a, c);
              alpha[a, c, s, DP_D_ART6to12m] := alpha_mult * DP.GetAdultMortByCD4WithART7to12(DP_Data, s, a, c);
              alpha[a, c, s, DP_D_ARTgt12m]  := alpha_mult * DP.GetAdultMortByCD4WithARTGT12(DP_Data, s, a, c);
            end;

        HVObj.PreInnerLoopCalcs(t);

      end;{if (GB.ModuleActive(GB_HV))}

//      CD4_Threshold := DP.GetCD4ThreshHoldAdults(t);
//      if CD4_Threshold = 200 then CD4_Threshold := 250;

      ltfu := -Ln(1.0 - DP.GetPercLostFollowUp(t) / 100); // convert the LTFU % to a per-capita annual rate
      for a := 15 to DP_MaxSingleAges do
      begin
        a10 := age10[a];
        for s := DP_Male to DP_Female do
          for c := DP_CD4_GT500 to DP_CD4_LT50 do
            for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
            begin
              if d = DP_D_ARTgt12m then
                d1 := DP_MortRates_GT12Mo
              else
                d1 := DP_MortRates_LT12Mo;
              mort := max(valueMinMort[a,s],
              alpha[a10,c,s,d] * (DP.GetMortalityRates(DP_Data,d1,t) * DP.GetMortalityRatesMultiplier));
              ARTPropMortOrLTFU[s,a,c,d] := mort + ltfu;
            end;
      end;

      for s := DP_Male to DP_Female do
      begin
        { Calculate number eligible for treatment }
        EligibleAdults[s,t] := 0;
        EligibleAdults[s,t - 1] := 0;
        if t > 2 then
          EligibleAdults[s, t - 2] := 0;
        ARTTwoYearsAgoBySex[s] := 0;

        for a := 15 to DP_MaxSingleAges do
        begin
          for c := DP_CD4_GT500 to DP_CD4_LT50 do
          begin
            if DP_CD4_lower_limit[c] < CD4_Threshold then {DP.GetCD4ThreshHoldAdults(t) then}
            begin
              EligibleAdults[s,t] := EligibleAdults[s,t] + pop1[t, s, a, c, DP_D_HIVPos];
              EligibleAdults[s,t - 1] := EligibleAdults[s,t - 1]
                                         + pop1[t - 1, s, a, c, DP_D_HIVPos];
              if t > 2 then
                EligibleAdults[s,t - 2] := EligibleAdults[s,t - 2]
                                            + pop1[t - 2, s, a, c, DP_D_HIVPos];
            end
            else {Add any special populations }
            begin
              { TH/HIV, SDC, SW, MSM, IDU, Other }
              for r := DP_EligTreatTB_HIV to DP_EligTreatOtherPop do
              begin
                if (r = DP_EligTreatTB_HIV) or
                   (r = DP_EligTreatDiscordantCouples) or
                   ((r = DP_EligTreatSexWorkers) and (s = DP_Female)) or
                   ((r = DP_EligTreatMSM) and (s = DP_Male)) or
                   (r = DP_EligTreatIDU) or
                   (r = DP_EligTreatOtherPop)
                then
                begin
                  if DP.GetPopsEligTreat(r).Eligible
                    and (DP.GetPopsEligTreat(r).Year <= (DP.GetFirstYear + t - 1))
                  then
                  begin
                    EligibleAdults[s,t] := EligibleAdults[s,t] +
                               DP.GetPopsEligTreat(r).PercentHIV / 100
                               * pop1[t, s, a, c, DP_D_HIVPos];

                    EligibleAdults[s,t - 1] := EligibleAdults[s,t - 1] +
                               DP.GetPopsEligTreat(r).PercentHIV / 100
                               * pop1[t - 1, s, a, c, DP_D_HIVPos];
                    if t > 2 then

                      EligibleAdults[s,t - 2] := EligibleAdults[s,t - 2] +
                                DP.GetPopsEligTreat(r).PercentHIV / 100
                                * pop1[t - 2, s, a, c, DP_D_HIVPos];
                  end;
                end;
              end; {r loop }
              if (s = DP_Female) and (a>=15) and (a<=49) and
                (DP.GetPopsEligTreat(DP_EligTreatPregnantWomen).Eligible)
                 and (DP.GetPopsEligTreat(DP_EligTreatPregnantWomen).Year
                      <= (DP.GetFirstYear + t - 1)) then
              begin
                EligibleAdults[s,t] := EligibleAdults[s,t] + EligiblePregnantWomen(a,t);
                EligibleAdults[s,t-1] := EligibleAdults[s,t-1] + EligiblePregnantWomen(a,t-1);
                if t>2 then
                  EligibleAdults[s,t-2] := EligibleAdults[s,t-2] + EligiblePregnantWomen(a,t-2);
              end;
            end; { Special pops }
          end; {c loop }
        end; { a loop }

        { Calculate number on treatment }
        sum1 := 0;
        for a := 15 to 80 do
        begin
          for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
          begin
            for c := DP_CD4_GT500 to DP_CD4_LT50 do
              sum1 := sum1 + pop1[t, s, a, c, d];
          end;
        end;
        DP.SetARTCalc_15Plus(DP_NumRecART, s, t, sum1);

        if t > 2 then
        begin
          if DP.GetHAARTBySexPerNum(t - 2, s) = DP_Percent then
//            ARTTwoYearsAgoBySex[s] := ConvertHAARTBySexPerToNum(DP, DP.GetHAARTBySex(t - 2, s), t, s)
            ARTTwoYearsAgoBySex[s] := DP.GetARTCalc_15Plus(DP_TotNeedART, s, t) *
                                      DP.GetHAARTBySex(t - 2, s) / 100
          else
            ARTTwoYearsAgoBySex[s] := DP.GetHAARTbySex(t - 2, s);

          if (GB.GetTempProj > 1) and (not GB.GetInCalcStateMode(p)) then
            if t = GBTempProjIdx(GB, GB_CalcYrIdxOfPrimProj) then
              ARTTwoYearsAgoBySexCalcState[s] := ARTTwoYearsAgoBySex[s];
        end
        else  (* t = 2 *)
          if GB.GetInCalcStateMode(p) then
            ARTTwoYearsAgoBySex[s] := ARTTwoYearsAgoBySexCalcState[s]
          else
            ARTTwoYearsAgoBySex[s] := 0;

        if DP.GetHAARTBySexPerNum(t - 1, s) = DP_Percent then
//            ARTLastYearBySex[s] := ConvertHAARTBySexPerToNum(DP, DP.GetHAARTBySex(t - 1, s), t, s)
          ARTLastYearBySex[s] := DP.GetARTCalc_15Plus(DP_TotNeedART, s, t)
            * DP.GetHAARTbySex(t - 1, s) / 100
        else
          ARTLastYearBySex[s] := DP.GetHAARTbySex(t - 1, s);
        if DP.GetHAARTbySexPerNum(t, s) = DP_Percent then
//          ARTThisYearBySex[s] := ConvertHAARTBySexPerToNum(DP, DP.GetHAARTBySex(t, s), t, s)
          ARTThisYearBySex[s] := DP.GetARTCalc_15Plus(DP_TotNeedART, s, t)
            * DP.GetHAARTbySex(t, s) / 100
        else
          ARTThisYearBySex[s] := DP.GetHAARTbySex(t, s);

        for a5 := DP_A0_4 to DP_A80_up do
          DP.SetNewlyNeedART(t, s, a5, 0);

        for a5 := DP_A0_4 to DP_A80_up do
          DP.SetAdultsChildrenStartingART(t, a5, s, 0);

       { Progress current infections in adults }
        for a5 := DP_A15_19 to DP_A80_up do
        begin
          DP.SetAIDSDeaths(t, a5, s, 0);
          DP.SetAIDSDeathsNoART(t, a5, s, 0);
          DP.SetAIDSDeathsART(t, a5, s, 0);
        end;

      end;{s loop}

      for t1 := 1 to DP_TimeStep do
      begin
        { Cache the number of adults not on ART by sex and CD4 count }
        for s := DP_Male to DP_Female do
          for c := DP_CD4_GT500 to DP_CD4_LT50 do
          begin
            AdultsNotOnART[s,c] := 0.0;
            for a := 15 to DP_MaxSingleAges do
              AdultsNotOnART[s,c] := AdultsNotOnART[s,c] + pop1[t,s,a,c,DP_D_HIVPos];
          end;

        for s := DP_Male to DP_Female do
        begin
          { Calculate number eligible for treatment }
          EligibleAdults[s,t] := 0;
          for a := 15 to DP_MaxSingleAges do
            for c := DP_CD4_GT500 to DP_CD4_LT50 do
              if DP_CD4_lower_limit[c] < DP.GetCD4ThreshHoldAdults(t) then
                EligibleAdults[s,t] := EligibleAdults[s,t] + pop1[t, s, a, c, DP_D_HIVPos]
              else
               for r := DP_EligTreatTB_HIV to DP_EligTreatOtherPop do
                if DP.GetPopsEligTreat(r).Eligible and (DP.GetPopsEligTreat(r).Year <= (DP.GetFirstYear + t - 1))then
                begin
                 if not(((r = DP_EligTreatMSM) and (s = DP_Female)) or ((r = DP_EligTreatSexWorkers) and (s = DP_male))) then
                 begin
                  EligibleAdults[s,t] := EligibleAdults[s,t] + DP.GetPopsEligTreat(r).PercentHIV / 100 * pop1[t, s, a, c, DP_D_HIVPos];
                 end;
                end;

          if DP.GetARTCoverageSelection = DP_NumNewARTPats then
          begin
            NewART := DP.GetNumNewARTPats(s,t) / DP_TimeStep;
          end
          else
          begin
            { Number of new ART patients is difference between need and those already on ART plus dropouts }
            ARTsurvivors := 0;
            for a := 15 to DP_MaxSingleAges do
              for c := DP_CD4_GT500 to DP_CD4_LT50 do
                for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
                begin
                  if (GB.GetYear(DP.proj, t) = 2020) and (DP.GetAdultPercentLTFUByMonth(0, t1) > 0) then
                    ARTPropMortOrLTFU[s,a,c,d] := ARTPropMortOrLTFU[s,a,c,d] - ltfu + DP.GetAdultPercentLTFUByMonth(0, t1) / 100;
                  ARTsurvivors := ARTsurvivors + pop1[t,s,a,c,d] * (1.0 - ARTPropMortOrLTFU[s,a,c,d] / DP_TimeStep);
                end;

            { Calculate number that should be on treatment }
            if DP.GetHAARTbySexPerNum(t, s) = DP_Percent then
              NeededART := ARTsurvivors + (ARTThisYearBySex[s] - ARTsurvivors)
                / DP_TimeStep * t1
            else
            begin

               { Adjust Dec 31 inputs if monthly data entered for 2020 }
               { For Spectrum 2020 calculations we only use the ART by month from Jan to June }
               if ({(t + DP.GetFirstYear - 1)}GB.GetYear(DP.proj, t) = 2020) then
               begin
                 sum := 0;
                 for c := GB_January to GB_June do
                   sum := sum + DP.GetAdultARTByMonth(s,0,c); { For second index, 0=2020, 1=2021 }
                 sum := sum * 5 / 6; {Convert from sum over 6 months to 5 time steps}
                 if (sum>0) then ARTThisYearBySex[s] := (sum - 3.5*ARTLastYearBySex[s])/1.5;
               end;
               { Adjust Dec 31 inputs if monthly data entered for 2020 }
               { For Spectrum 2021 calculations we use the ART by month from July to Dec }
               if ({(t + DP.GetFirstYear - 1)}GB.GetYear(DP.proj, t) = 2021) then
               begin
                 sum := 0;
                 for c := GB_July to GB_December do
                   sum := sum + DP.GetAdultARTByMonth(s,0,c); { For second index, 0=2020, 1=2021 }
                 sum := sum * 5 / 6; {Convert from sum over 6 months to 5 time steps}
                 if (sum>0) then ARTLastYearBySex[s] := (sum - ARTTwoYearsAgoBySex[s])/4 ;
               end;

              if t1 < (DP_TimeStep div 2) then
                NeededART := ARTTwoYearsAgoBySex[s] + (ARTLastYearBySex[s] - ARTTwoYearsAgoBySex[s])
                  / DP_TimeStep * (t1 + (DP_TimeStep / 2))
              else
                NeededART := ARTLastYearBySex[s] + (ARTThisYearBySex[s] - ARTLastYearBySex[s]) / DP_TimeStep *
                  (t1 - (DP_TimeStep / 2));
            end;

            NewART := NeededART - ARTsurvivors;
          end;

          if NewART < 0 then
          begin
            NewART := 0;
            DP_AddLogEntry(GB, DP_DPAIM, 'DPPROJA', '~1-30', 'NewART', DP_LTZ, NewART,' ');
          end;
          { If number of new ART patients is greater than the number eligible
            set it equal to eligibility }
          if NewART > EligibleAdults[s,t] then
          begin
            NewART := EligibleAdults[s,t];
            DP_AddLogEntry(GB, DP_DPAIM, 'DPPROJA', '~1-40', 'NewART',
              'Num new ART patients greater than num elibigle', NewART,
              't = ' + IntToStr(t));
          end;

          if NewART > 0 then CalcNewARTPatients(s);

          for a := 15 to DP_MaxSingleAges do
          begin
            { Calculate index of 5 year age groups }
            a5 := trunc(a / 5) + 1;
            { Calculate index for 10 year adult age groups }
            a10 := age10[a];

            for c := DP_CD4_GT500 to DP_CD4_LT50 do
            begin
              { For those not on ART }
              d := DP_D_HIVPos;

              { Determine the number newly starting ART }
              { Allocate new patients in CD4 category by age according the
                distribution of the HIV+ population not on ART }
              if (AdultsNotOnART[s,c] > 0) then
              begin
                if (DP.GetARTCovByAge(s, a5, t)>0) then
                  {Calculate StartART as the number needing ART (pop x cov) minus those on ART divided
                   by DP_TimeStep-t1+1 to adjust to time within the t1 loop }
                  StartART := (pop1[t, s, a, c, DP_D_All] * DP.GetARTCovByAge(s, a5, t)/100 -
                              (pop1[t,s,a,c,DP_D_ARTlt6m]+pop1[t,s,a,c,DP_D_ART6to12m]+pop1[t,s,a,c,DP_D_ARTgt12m])) /
                              (DP_TimeStep - t1 + 1)
                else StartART := NewPatients[c,s] * pop1[t, s, a, c, d] / AdultsNotOnART[s,c]
              end
              else
                StartART := 0;
{ New code }
              if DP.GetARTCoverageSelection = DP_CD4Number then
              begin
                sum10 := AdultsNotOnART[DP_Male,c] + AdultsNotOnART[DP_Female,c];
                if (sum10 = 0) or (DP.GetCD4Coverage(DP_CD4Number,t,c) = 0) then
                  StartART := 0
                else
                begin
                  { Calculate number already on ART for this age, CD4 category and sex}
                  sum := pop1[t,s,a,c,DP_D_ARTlt6m] + pop1[t,s,a,c,DP_D_ART6to12m] + pop1[t,s,a,c,DP_D_ARTgt12m];
                  StartART := (DP.GetCD4Coverage(DP_CD4Number,t,c)
                              * pop1[t, s, a, c, d] / sum10
                              - sum ) / (DP_TimeStep - t1 + 1);
                end;
              end;
              if StartART < 0 then
                StartART := 0;

              if DP.GetARTCoverageSelection = DP_CD4Percent then
                StartART := 0;

              { In CD4 > 500 new entrants are new infections
                in all other CD4 categories new entrants are those
                progressing from next highest CD4 category
                New entrants also include those dropping out of ART }
//              if c = DP_CD4_GT500 then
//                { Add ART dropouts, new infections are calculated at the end }
//                entrants[c, d] := (pop1[t, s, a, c, DP_D_ARTlt6m]
//                  + pop1[t, s, a, c, DP_D_ART6to12m] + pop1[t, s, a, c, DP_D_ARTgt12m])
//                  * DP.GetPercLostFollowUp(t)/100 / DP_TimeStep
//              else
//                { Add progressors from higher CD4 category and ART dropouts }
//                entrants[c, d] := pop1[t, s, a, c - 1, d] * lambda[a10, c - 1, s]
//                  + (pop1[t, s, a, c, DP_D_ARTlt6m] + pop1[t, s, a, c, DP_D_ART6to12m]
//                  + pop1[t, s, a, c, DP_D_ARTgt12m]) * DP.GetPercLostFollowUp(t)/100
//                  / DP_TimeStep;

              if c = DP_CD4_GT500 then
                { Add ART dropouts, new infections are calculated at the end }
                { CD4>500 receives ART drop outs from all durations for those starting at CD4>500
                  and ART drop outs from those starting CD4_350-500 on ART for more than one year }
                entrants[c, d] := (pop1[t, s, a, c, DP_D_ARTlt6m]
                  + pop1[t, s, a, c, DP_D_ART6to12m] + pop1[t, s, a, c, DP_D_ARTgt12m]
                  + pop1[t, s, a, c+1, DP_D_ARTgt12m])
                  * ltfu / DP_TimeStep
              else if c = DP_CD4_LT50 then
                { Add ART dropouts from those starting at CD4<50 only for those
                  on ART < 1 year. Those on ART>1 year go to CD4_50_99 group}
                entrants[c, d] := pop1[t, s, a, c - 1, d] * lambda[a10, c - 1, s]
                  + (pop1[t, s, a, c, DP_D_ARTlt6m] + pop1[t, s, a, c, DP_D_ART6to12m])
                  * ltfu / DP_TimeStep
              else
                { Add progressors from higher CD4 category and ART dropouts }
                { ART drop outs come from those starting ART at one lower CD4 category
                  for those on ART more than 1 year }
                entrants[c, d] := pop1[t, s, a, c - 1, d] * lambda[a10, c - 1, s]
                  + (pop1[t, s, a, c, DP_D_ARTlt6m] + pop1[t, s, a, c, DP_D_ART6to12m]
                  + pop1[t, s, a, c + 1, DP_D_ARTgt12m])
                  * ltfu / DP_TimeStep;

              { Exits are those who leave because they progress to the
                next lower CD4 category, die or start on ART }
              v1 := pop1[t,s,a,c,DP_D_ARTlt6m] + pop1[t,s,a,c,DP_D_ART6to12m] + pop1[t,s,a,c,DP_D_ARTgt12m];
              v2 := v1 + pop1[t,s,a,c,DP_D_HIVPos];
              if (v2>0) then v3 := 1 - v1 / v2 else v3 := 1;
              exits[c, d] := pop1[t, s, a, c, d] *
                (lambda[a10, c, s] + max(mu[a10, c, s] * v3, valueMinMort[a,s]) / DP_TimeStep)
                + min(StartART, pop1[t,s,a,c,d]);
              {exits[c, d] := pop1[t, s, a, c, d] *
                (lambda[a10, c, s] + max(mu[a10, c, s],valueMinMort[a,s]) / DP_TimeStep)
                + min(StartART, pop1[t,s,a,c,d]);}

              { Newly needing ART }
              if DP_CD4_lower_limit[c] = DP.GetCD4ThreshHoldAdults(t) then
                DP.SetNewlyNeedART(t, s, a5, DP.GetNewlyNeedART(t, s, a5)
                    + pop1[t, s, a, c, d] * lambda[a10, c, s]);

              { Save AIDS deaths }
              temp := max(mu[a10, c, s] * v3, valueMinMort[a,s]) * pop1[t, s, a, c, d] / DP_TimeStep;
              AdultAIDSDeaths[a,s,DP_NoTreat] := AdultAIDSDeaths[a,s,DP_NoTreat] + temp;

              if (GB.ModuleActive(p, GB_HV)) then
              begin
                temp_cd4 := max(mu[a10, c, s],valueMinMort[a,s]) * pop1[t, s, a, c, d]/DP_TimeStep;//need a t1 rate
                DP.SetAIDSCd4MortNoART(t, t1, a5, s, c, DP.GetAIDSCd4MortNoART(t, t1, a5, s, c)+temp_cd4);

                temp_cd4:= pop1[t, s, a, c, d] * lambda[a10, c, s];//lambda already a t1 rate
                DP.SetAIDSCd4Trans(t, t1, a5, s, c,  DP.GetAIDSCd4Trans(t, t1, a5, s, c)+temp_cd4);
              end;

              { For those on ART }

        (*EPW*)
              DP.SetNewlyOnART(t,t1,s,a5,c,DP.GetNewlyOnART(t,t1,s,a5,c)+StartART);
              DP.SetAdultsChildrenStartingART(t, a5, s, DP.GetAdultsChildrenStartingART(t, a5, s) + StartART);

              d := DP_D_ARTlt6m;
              entrants[c,d] := min(StartART,pop1[t,s,a,c,DP_D_HIVPos]);
              exits[c,d] := pop1[t,s,a,c,d] * (ARTPropMortOrLTFU[s,a,c,d] + 12 / 6) / DP_TimeStep;

              d := DP_D_ART6to12m;
              entrants[c,d] := pop1[t,s,a,c,DP_D_ARTlt6m] * (12 / 6) / DP_TimeStep;
              exits[c,d] := pop1[t,s,a,c,d] * (ARTPropMortOrLTFU[s,a,c,d] + 12 / 6) / DP_TimeStep;

              d := DP_D_ARTgt12m;
              entrants[c,d] := pop1[t,s,a,c,DP_D_ART6to12m] * (12 / 6) / DP_TimeStep;
              exits[c,d] := pop1[t,s,a,c,d] * ARTPropMortOrLTFU[s,a,c,d] / DP_TimeStep;

              if c < DP_CD4_LT50 then
              begin
                d := DP_D_ARTlt6m;
                entrants[c,d] := entrants[c,d] + pop1[t,s,a,c+1,d] * beta1[a10,c+1,s] / DP_TimeStep;

                d := DP_D_ART6to12m;
                entrants[c,d] := entrants[c,d] + pop1[t,s,a,c+1,d] * beta2[a10,c+1,s] / DP_TimeStep;

                d := DP_D_ARTgt12m;
                entrants[c,d] := entrants[c,d] + pop1[t,s,a,c+1,d] * beta2[a10,c+1,s] / DP_TimeStep;
              end;

              if c > DP_CD4_GT500 then
              begin
                d := DP_D_ARTlt6m;
                exits[c,d] := exits[c,d] + pop1[t,s,a,c,d] * beta1[a10,c,s] / DP_TimeStep;

                d := DP_D_ART6to12m;
                exits[c,d] := exits[c,d] + pop1[t,s,a,c,d] * beta2[a10,c,s] / DP_TimeStep;

                d := DP_D_ARTgt12m;
                exits[c,d] := exits[c,d] + pop1[t,s,a,c,d] * beta2[a10,c,s] / DP_TimeStep;
              end;

              for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
              begin
                { Save AIDS deaths }
                if d = DP_D_ARTgt12m then
                  d1 := DP_MortRates_GT12Mo
                else
                  d1 := DP_MortRates_LT12Mo;
                temp := max(valueMinMort[a,s], alpha[a10, c, s, d] *
                  (DP.GetMortalityRates(DP_Data,d1,t) * DP.GetMortalityRatesMultiplier))
                   * pop1[t, s, a, c, d] / DP_TimeStep;
                AdultAIDSDeaths[a,s,DP_OnART] := AdultAIDSDeaths[a,s,DP_OnART] + temp;

                if (GB.ModuleActive(p, GB_HV)) then
                begin
                  temp_cd4 := max(valueMinMort[a,s], alpha[a10, c, s, d] *
                    (DP.GetMortalityRates(DP_Data,d1,t) * DP.GetMortalityRatesMultiplier)) *
                    pop1[t, s, a, c, d] / DP_TimeStep;//need a t1 rate
                  DP.SetAIDSCd4MortART(t, t1, a5, s, c, DP.GetAIDSCd4MortART(t, t1, a5, s, c)+temp_cd4);
                end;

              end; { d loop }
            end; { c loop }

            for c := DP_CD4_GT500 to DP_CD4_LT50 do
            begin
              { For those not on ART }
              d := DP_D_HIVPos;
              pop1[t, s, a, c, d] := pop1[t, s, a, c, d] + entrants[c, d] - exits[c, d];
              { For those on ART }
              for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
              begin
                pop1[t, s, a, c, d] := pop1[t, s, a, c, d]
                            + entrants[c, d] - exits[c, d];
              end; { d loop }
            end; { c loop }
          end; { a loop }

          { Add new ART if ART coverage by CD4 count }
          if DP.GetARTCoverageSelection = DP_CD4Percent then
          begin
            d := DP_D_HIVPos; {those who are not on ART }
            for c := DP_CD4_GT500 to DP_CD4_LT50 do
            begin
              StartART := 0;
              TotalPopAge := 0;
              for a := 15 to DP_MaxSingleAges do
              begin
                { Calculate index of 5 year age groups }
                a5 := trunc(a / 5) + 1;
                { Calculate index for 10 year adult age groups }
                a10 := age10[a];

                if (GB.ModuleActive(p, GB_HV)) then
                StartART := StartART +
                           HVObj.GetDP_RiskGrpARTPerc(s,c,t)
                           * (pop1[t, s, a, c, d]
                              + pop1[t, s, a, c, DP_D_ARTlt6m]
                              + pop1[t, s, a, c, DP_D_ART6to12m]
                              + pop1[t, s, a, c, DP_D_ARTgt12m])
                            - (pop1[t, s, a, c, DP_D_ARTlt6m]
                              + pop1[t, s, a, c, DP_D_ART6to12m]
                              + pop1[t, s, a, c, DP_D_ARTgt12m])
                            + (pop1[t, s, a, c, DP_D_ARTlt6m]
                                + pop1[t, s, a, c, DP_D_ART6to12m]
                                + pop1[t, s, a, c, DP_D_ARTgt12m])
                              * ltfu / DP_TimeStep
                else
                  StartART := StartART +
                             DP.GetCD4Coverage(DP_CD4Percent,t,c) / 100
                             * (pop1[t, s, a, c, d]
                                + pop1[t, s, a, c, DP_D_ARTlt6m]
                                + pop1[t, s, a, c, DP_D_ART6to12m]
                                + pop1[t, s, a, c, DP_D_ARTgt12m])
                             - (pop1[t, s, a, c, DP_D_ARTlt6m]
                                + pop1[t, s, a, c, DP_D_ART6to12m]
                                + pop1[t, s, a, c, DP_D_ARTgt12m])
                              + (pop1[t, s, a, c, DP_D_ARTlt6m]
                                 + pop1[t, s, a, c, DP_D_ART6to12m]
                                 + pop1[t, s, a, c, DP_D_ARTgt12m])
                                * ltfu / DP_TimeStep;
                TotalPopAge := TotalPopAge + pop1[t, s, a, c, d];


              end; {a loop}
              if StartART < 0 then StartART := 0;
              for a := 15 to DP_MaxSingleAges do
              begin
                if TotalPopAge > 0 then
                begin
                //modify nART:
                  nART := pop1[t, s, a, c, d] * StartART / TotalPopAge;
                 //nART := min(pop1[a, t, s, c, d] , pop1[a, t, s, c, d] * StartART / TotalPopAge);
                  pop1[t, s, a, c, d] := pop1[t, s, a, c, d] - nART;
                  pop1[t, s, a, c, DP_D_ARTlt6M] := pop1[t, s, a, c, DP_D_ARTlt6M] + nART;
                end;
              end; {a loop}
            end; {c loop}
          end; {If coverage by CD4 count}

      (*EPW Start*)
          //The Null scenario - no ART
          if GB.GetSpectrumMode=GB_ChoiceMode then
          begin


            if (GetRNBaseScenHVType = DP_SG_NullScenario)
            or (GetRNScaleScenHVType = DP_SG_NullScenario) then
            begin

              if(DP.GetHAARTBySex(t, s)=0) and (t>=GetRNBaseYearSGIndex) then
              begin
                for a := 0 to 14 do
                begin
                  for c := DP_CD4_Ped_GT1000 to DP_CD4_Ped_LT200  do
                  begin
                    d := DP_D_HIVPos; {those who are not on ART }
                    //move from on ART to no on ART
                    pop1[t, s, a, c, d] := pop1[t, s, a, c, d]+
                                           pop1[t, s, a, c, DP_D_ARTlt6m]+
                                           pop1[t, s, a, c, DP_D_ART6to12m]+
                                           pop1[t, s, a, c, DP_D_ARTgt12m];

                  pop1[t, s, a, c, DP_D_ART6to12m]:=0;
                  pop1[t, s, a, c, DP_D_ARTgt12m]:=0;
                  end;//for c

                  for c := DP_CD4_GT500 to DP_CD4_LT50 do
                  begin
                    d := DP_D_HIVPos; {those who are not on ART }
                    //move from on ART to no on ART
                    pop1[t, s, a, c, d] := pop1[t, s, a, c, d]+
                                           pop1[t, s, a, c, DP_D_ARTlt6m]+
                                           pop1[t, s, a, c, DP_D_ART6to12m]+
                                           pop1[t, s, a, c, DP_D_ARTgt12m];

                    pop1[t, s, a, c, DP_D_ARTlt6m]:=0;
                    pop1[t, s, a, c, DP_D_ART6to12m]:=0;
                    pop1[t, s, a, c, DP_D_ARTgt12m]:=0;
                  end;//for c

                end;//for a

                for a := 15 to DP_MaxSingleAges do
                  for c := DP_CD4_GT500 to DP_CD4_LT50 do
                  begin
                    d := DP_D_HIVPos; {those who are not on ART }
                    //move from on ART to no on ART
                    pop1[t, s, a, c, d] := pop1[t, s, a, c, d]+
                                           pop1[t, s, a, c, DP_D_ARTlt6m]+
                                           pop1[t, s, a, c, DP_D_ART6to12m]+
                                           pop1[t, s, a, c, DP_D_ARTgt12m];

                    pop1[t, s, a, c, DP_D_ARTlt6m]:=0;
                    pop1[t, s, a, c, DP_D_ART6to12m]:=0;
                    pop1[t, s, a, c, DP_D_ARTgt12m]:=0;
                  end;//for c/a
              end;//if(DP.GetHAARTBySex(t, s)=0) and (t>=GetRNBaseYearSGIndex) then
            end;//if null scenario
          end;//if in choice mode
      (*EPW End*)
        end; {s loop}

        (*if Goals is active, run Goals InnerLoopCalcs*)
        if (GB.ModuleActive(p, GB_HV)) then
        begin
            (*set number reveiving ART for Goals*)
            CalcGoalsReceiveART(t);
            HVObj.InnerLoopCalcs(t,t1);
            CalcNewHIVGoals(t1);// add new infections from Goals to AIM
        end
        else if (GB.ModuleActive(p, GB_HA)) then
        begin
          //if ((t1 = 1) and (t = 2020-1970+1)) then
          if (t1 = 1) then
          begin
            HAObj.CalcHIVTransmitters(t, pop1);
          end;

          if(GetHAData(p).GetWriteTransmittersFlag=TRUE) then
          if ((t1 = 1) and (t = 2020-GetGBFirstYear(p)+1)) then
            HAObj.WriteHIVTransmitters(t, p);

          CalcInfectionsByAge(p, t, t1);

        end;


      end; { t1 loop }

      { Bookkeeping: transfer cached adult AIDS deaths to DPData arrays }
      for a := 15 to DP_MaxSingleAges do
      begin
        a5 := a Div 5 + 1;
        for s := DP_Male to DP_Female do
        begin
          DP.SetAIDSDeathsNoART(t, a5, s, DP.GetAIDSDeathsNoART(t, a5, s) + AdultAIDSDeaths[a,s,DP_NoTreat]);
          DP.SetAIDSDeathsART(t, a5, s, DP.GetAIDSDeathsART(t, a5, s) + AdultAIDSDeaths[a,s,DP_OnART]);
          DP.SetAIDSDeathsARTSingleAge(t, s, a + 1, AdultAIDSDeaths[a, s, DP_OnART]);
          AIDSDeathsByAge[a,s,t] := AdultAIDSDeaths[a,s,DP_NoTreat] + AdultAIDSDeaths[a,s,DP_OnART];
        end;
      end;
      for a5 := DP_A15_19 to DP_A80_up do
        for s := DP_Male to DP_Female do
          DP.SetAIDSDeaths(t, a5, s, DP.GetAIDSDeathsNoART(t, a5, s) + DP.GetAIDSDeathsART(t, a5, s));

      SumART_BS := 0;
      SumART_15_49_BS := 0;
      for s := DP_Male to DP_Female do
      begin
        sumART := 0;
        for a := 15 to 80 do
          for c := DP_CD4_GT500 to DP_CD4_LT50 do
            for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
              sumART := sumART + pop1[t, s, a, c, d];
        DP.SetARTCalc_15Plus(DP_NumRecART, s, t, sumART);
        SumART_BS := SumART_BS + SumART;

        SumART_15_49 := 0;
        for a := 15 to 49 do
          for c := DP_CD4_GT500 to DP_CD4_LT50 do
            for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
              SumART_15_49 := SumART_15_49 + pop1[t, s, a, c, d];
        SumART_15_49_BS := SumART_15_49_BS + SumART_15_49;

        if SumART = 0 then
          DP.SetPercART50Plus(s, t, 0)
        else
          DP.SetPercART50Plus(s, t, (SumART - SumART_15_49) / SumART);
      end;
      if SumART_BS = 0 then
        DP.SetPercART50Plus(DP_BothSexes, t, 0)
      else
        DP.SetPercART50Plus(DP_BothSexes, t, (SumART_BS - SumART_15_49_BS) / SumART_BS);

      DP.SetARTCalc_15Plus(DP_NumRecART, DP_BothSexes, t,
              DP.GetARTCalc_15Plus(DP_NumRecART, DP_Male, t)
               + DP.GetARTCalc_15Plus(DP_NumRecART, DP_Female, t));

      (* if Goals is active, run Goals PostInnerLoopCalcs*)
      if (GB.ModuleActive(p, GB_HV)) then
        HVObj.PostInnerLoopCalcs(t);

      { Sum for all status and duration states }
      SumPopulationByHIVStatus(t);

      { Recalculate Need for ART if number on ART is greater than need }
      sum1 := 0;
      for s := DP_Male to DP_Female do
        for a := 15 to 80 do
          for c := DP_CD4_GT500 to DP_CD4_LT50 do
            sum1 := sum1 + pop1[t,s,a,c,DP_D_All] - pop1[t,s,a,c,DP_D_HIVPos];
      if sum1 > DP.GetARTCalc_15Plus(DP_TotNeedART, DP_BothSexes, t) then
        CalcNeedForART(pop1);

      { Include calculation of new adult infections }
      {If Goals is active, then new infections added as part of InnerLoopCalcs}
      if ((not GB.ModuleActive(p, GB_HV)) and (not GB.ModuleActive(p, GB_HA))) then
        CalcNewAdultInfections(t);

      { Count those newly in need of ART under treat-all }
      if DP.GetCD4ThreshHoldAdults(t) >= DP_CD4_lower_limit[DP_CD4_GT500-1]-1 then
        for s := GB_Male to GB_Female do
          for a5 := DP_A15_19 to DP_A80_up do
            DP.SetNewlyNeedART(t, s, a5, DP.GetNewInfections(t, a5, s));

      { Start child calculations }
      { Include calculation of HIV+ births }
      CalcHIVBirths(t);

      DP.SetAIDSBirths(t, DP_BothSexes, HIVbirths);
      DP.SetAIDSBirths(t, DP_Male, HIVbirths * GetDP_SexBirthRatioPercent(DP, t, DP_Male));
      DP.SetAIDSBirths(t, DP_Female, HIVbirths * GetDP_SexBirthRatioPercent(DP, t, DP_Female));

      { Distribute new infections by CD4 percent }
      for s := DP_Male to DP_Female do
      begin
        pop1[t, s, 0, DP_H_All, DP_P_Perinatal] := 0;
        pop1[t, s, 0, DP_H_All, DP_P_BF0] := 0;
        pop1[t, s, 0, DP_H_All, DP_P_BF7] := 0;
        pop1[t, s, 0, DP_H_All, DP_P_BF12] := 0;
      end;
      for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
      begin
        for s := DP_Male to DP_Female do
        begin
          { Perinatal infections }
          pop1[t, s, 0, c, DP_P_Perinatal] := DP.GetAIDSBirths(t, s)
                                         * CD_perc_at_infection[c] / 100;
          pop1[t, s, 0, DP_H_All, DP_P_Perinatal] := pop1[t, s, 0, DP_H_All, DP_P_Perinatal]
                                                + pop1[t, s, 0, c, DP_P_Perinatal];
          { BF infections at 0-6 months plus incident infections }
          pop1[t, s, 0, c, DP_P_BF0] := (NewInfFromBFLT6 + IncidentInfectionsBF)
                                               * GetDP_SexBirthRatioPercent(DP, t, s)
                                               * CD_perc_at_infection[c] / 100;
          pop1[t, s, 0, DP_H_All, DP_P_BF0] := pop1[t, s, 0, DP_H_All, DP_P_BF0]
                                                + pop1[t, s, 0, c, DP_P_BF0];

          { BF infections at 7-12 months }
          pop1[t, s, 0, c, DP_P_BF7] := NewInfFromBF6to12
                                               * GetDP_SexBirthRatioPercent(DP, t, s)
                                               * CD_perc_at_infection[c] / 100;
          pop1[t, s, 0, DP_H_All, DP_P_BF7] := pop1[t, s, 0, DP_H_All, DP_P_BF7]
                                                + pop1[t, s, 0, c, DP_P_BF7];
        end;

      end;

      pop1[t,DP_BothSexes,0,DP_H_All,DP_P_Perinatal] := pop1[t,DP_Male,0,DP_H_All,DP_P_Perinatal]
                                                      + pop1[t,DP_Female,0,DP_H_All,DP_P_Perinatal];
      pop1[t,DP_BothSexes,0,DP_H_All,DP_P_BF0] := pop1[t,DP_Male,0,DP_H_All,DP_P_BF0]
                                                + pop1[t,DP_Female,0,DP_H_All,DP_P_BF0];
      pop1[t,DP_BothSexes,0,DP_H_All,DP_P_BF7] := pop1[t,DP_Male,0,DP_H_All,DP_P_BF7]
                                                + pop1[t,DP_Female,0,DP_H_All,DP_P_BF7];


      for s := DP_BothSexes to DP_Female do
      begin
  {      pop1[0, t, s, DP_H_No, DP_H_No] := pop1[0, t, s, DP_H_No, DP_H_No]
                                         - pop1[0, t, s, DP_P_Perinatal, DP_H_No]
                                         - pop1[0, t, s, DP_P_BF0, DP_H_No]
                                         - pop1[0, t, s, DP_P_BF7, DP_H_No];
  }
         pop1[t, s, 0, DP_H_No, DP_D_HIVNeg] := pop1[t, s, 0, DP_H_No, DP_D_HIVNeg]
                                         - pop1[t, s, 0, DP_H_All, DP_P_Perinatal]
                                         - pop1[t, s, 0, DP_H_All, DP_P_BF0]
                                         - pop1[t, s, 0, DP_H_All, DP_P_BF7];
        if pop1[t, s, 0, DP_H_No, DP_D_HIVNeg] < 0 then
          pop1[t, s, 0, DP_H_No, DP_D_HIVNeg] := 0;
        DP.SetNewInfections(t, DP_A0_4, s, DP.GetNewInfections(t, DP_A0_4, s)
                                         + DP.GetAIDSBirths(t, s)
                                         + (NewInfFromBFLT6 + NewInfFromBF6to12
                                         + IncidentInfectionsBF)
                                         * GetDP_SexBirthRatioPercent(DP,t,s));
        DP.SetNewInfectionsBySingleAge(t, s, DP_A0,
                                       DP.GetAIDSBirths(t, s)
                                       + (NewInfFromBFLT6 + NewInfFromBF6to12
                                       + IncidentInfectionsBF)
                                       * GetDP_SexBirthRatioPercent(DP,t,s));
      end;
      DP.SetNewInfections(t, DP_A0_4, DP_BothSexes,
        DP.GetNewInfections(t, DP_A0_4, DP_Male) + DP.GetNewInfections(t, DP_A0_4, DP_Female));

      { Calculate new infections through breastfeeding for those over the age of one }
//      if DP.GetBirths(t, DP_AllAges, DP_BothSexes) = 0 then
//        PercentExposed := 0
//      else
//        PercentExposed := (births - HIVbirthsToExistingInfections - CumNewInfFromBF) / DP.GetBirths
//          (t, DP_AllAges, DP_BothSexes);
      CalcBFtransmission(p, t, 7, 12);
      a := 1;
      for s := DP_Male to DP_Female do
      begin
        if (pop1[t, DP_Male, a, DP_H_No, DP_D_HIVNeg] + pop1[t, DP_Female, a, DP_H_No, DP_D_HIVNeg]) = 0
          then NewInfFromBFGT12 := 0
          else NewInfFromBFGT12 := births * BFTR
                                   * pop1[t, s, a, DP_H_No, DP_D_HIVNeg]
                                   / (pop1[t, DP_Male, a, DP_H_No, DP_D_HIVNeg] + pop1[t, DP_Female, a, DP_H_No, DP_D_HIVNeg]);
        FirstYearNewInfFromBF := NewInfFromBFGT12;
        for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
          pop1[t, s, a, c, DP_P_BF12] := NewInfFromBFGT12 *  CD_perc_at_infection[c] / 100;
        pop1[t, s, a, DP_H_No, DP_D_HIVNeg] := pop1[t, s, a, DP_H_No, DP_D_HIVNeg] - NewInfFromBFGT12;
        DP.SetNewInfections(t, DP_A0_4, s, DP.GetNewInfections(t, DP_A0_4, s)
                                             + NewInfFromBFGT12);
  {                                           + NewInfFromBFGT12 * DP.GetSexBirthRatioPercent(t,s));}
        DP.SetNewInfectionsBySingleAge(t, s, DP_A1, NewInfFromBFGT12);
        CumNewInfFromBF := CumNewInfFromBF + NewInfFromBFGT12;
      end;

//      PercentExposed := PercentExposed * (1 - BFTR);
      CalcBFtransmission(p, t, 13, 18);
      a := 2;
      for s := DP_Male to DP_Female do
      begin
        if (pop1[t, DP_Male, a, DP_H_No, DP_D_HIVNeg] + pop1[t, DP_Female, a, DP_H_No, DP_D_HIVNeg]) = 0
          then NewInfFromBFGT12 := 0
          else NewInfFromBFGT12 := births * BFTR
                                   * pop1[t, s, a, DP_H_No, DP_D_HIVNeg]
                                   / (pop1[t, DP_Male, a, DP_H_No, DP_D_HIVNeg] + pop1[t, DP_Female, a, DP_H_No, DP_D_HIVNeg]);
//        NewInfFromBFGT12 := pop1[t, s, a, DP_H_No, DP_D_HIVNeg] * PercentExposed * BFTR;
        SecondYearNewInfFromBF := NewInfFromBFGT12;
        for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
          pop1[t, s, a, c, DP_P_BF12] := pop1[t, s, a, c, DP_P_BF12]
                                         + NewInfFromBFGT12 *  CD_perc_at_infection[c] / 100;
        pop1[t, s, a, DP_H_No, DP_D_HIVNeg] := pop1[t, s, a, DP_H_No, DP_D_HIVNeg] - NewInfFromBFGT12;
        DP.SetNewInfections(t, DP_A0_4, s, DP.GetNewInfections(t, DP_A0_4, s)
                                             + NewInfFromBFGT12);
  {                                           + NewInfFromBFGT12 * DP.GetSexBirthRatioPercent(t,s));}
        DP.SetNewInfectionsBySingleAge(t, s, DP_A2, NewInfFromBFGT12);
        CumNewInfFromBF := CumNewInfFromBF + NewInfFromBFGT12;
      end;

      if HIVbirths > 0 then
        DP.SetPerinatalTransmission(t, PTR1 * (1 + CumNewInfFromBF/ HIVbirths));

      DP.SetNewInfections(t, DP_A0_4, DP_BothSexes,
                         DP.GetNewInfections(t, DP_A0_4, DP_Male)
                        + DP.GetNewInfections(t, DP_A0_4, DP_Female));

  {
  if t = 6 then
  begin
  AssignFile(OutFile,'Children.txt');
  Rewrite(OutFile);
  writeln(Outfile, 's  a  c  Perinatal, BF0, BF7, BF12');
  for s := DP_Male to DP_Female do
  for a := 0 to 4 do
  for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
  writeln(Outfile, s:2, a:3, c: 2, pop1[a,t,s,c,DP_P_Perinatal]:8:3,pop1[a, t, s, c, DP_P_BF0]:8:3,
                                   pop1[a,t,s,c,DP_P_BF7]:8:3, pop1[a,t,s,c,DP_P_BF12]:8:3);
  writeln(OutFile, 'births ',HIVbirths:4:3, ' BF0-6 ', NewInfFromBFLT6:4:3,' BF6-12 ',NewInfFromBF6to12:4:3,' FY BF12+ ',
                    FirstYearNewInfFromBF:4:3,' SY BF12+ ',SecondYearNewInfFromBF:4:3,' New Inf ',DP.GetNewInfections(t, DP_A0_4, DP_Male):8:3,
                    DP.GetNewInfections(t, DP_A0_4, DP_Female):8:3);
  CloseFile(OutFile);
  end;
  }
      { Progress HIV+ children not on ART }
      for a := 0 to 14 do
      begin
        a5 := trunc(a/5) + 1;
        a1 := DP_CD4_0_2;
        case a of
          0 .. 2:  a1 := DP_CD4_0_2;
          3 .. 4:  a1 := DP_CD4_3_4;
          5 .. 14: a1 := DP_CD4_5_14;
        end;









        { For those not on ART }
        ART_Excess := 0;
        for s := DP_Male to DP_Female do
        begin
          v2 := 0;
          for d := DP_P_Perinatal to DP_P_BF12 do
          begin
            if a < 5 then
            begin
              for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
              begin
                { entrants[c,d] are progressors in from higher CD4 category }
                if c <> DP_CD4_Per_GT30 then
                  entrants[c,d] := v2
                else
                  entrants[c,d] := 0;
                { v1 is AIDS deaths }
                v1 := pop1[t,s,a,c,d] * pmu[a1,c,s,d];
                { Adjust mortality for effects of cotrimoxazole }
                v1 := v1 * (1 - CTXeff * CTXcoverage);
                { v2 is progression out to lower CD4 category }
                v2 := (pop1[t,s,a,c,d] - v1/2) * plambda[a1,c,s,d];
                exits[c,d] := v1 + v2; { + NewChildART[a,s,c,d] + ART_Excess;}
                DP.SetAIDSDeaths(t,a5,s, DP.GetAIDSDeaths(t,a5,s) + v1);
                AidsDeathsByAge[a,s,t] := AidsDeathsByAge[a,s,t] + v1;
                pop1[t,s,a,c,d] := pop1[t,s,a,c,d] + entrants[c,d] - exits[c,d];
                if pop1[t,s,a,c,d] < 0 then
                begin
                  ART_Excess := ART_Excess - pop1[t,s,a,c,d];
                  pop1[t,s,a,c,d] := 0;
                end;
              end;
            end
            else {a > 4}
            begin
              for c := DP_CD4_Ped_GT1000 to DP_CD4_Ped_LT200 do
              begin
                { entrants[c,d] are progressors in from higher CD4 category }
                if c <> DP_CD4_Ped_GT1000 then
                  entrants[c,d] := v2
                else
                  entrants[c,d] := 0;
                { v1 is AIDS deaths }
                v1 := pop1[t,s,a,c,d] * pmu[a1,c,s,d];
                { Adjust mortality for effects of cotrimoxazole }
                v1 := v1 * (1 - CTXeff * CTXcoverage);
                { v2 is progression out to lower CD4 category }
                v2 := (pop1[t,s,a,c,d] - v1/2) * plambda[a1,c,s,d];
                exits[c,d] := v1 + v2 {+ NewChildART[a,s,c,d]};
                DP.SetAIDSDeaths(t,a5,s, DP.GetAIDSDeaths(t,a5,s) + v1);
                AidsDeathsByAge[a,s,t] := AidsDeathsByAge[a,s,t] + v1;
                pop1[t,s,a,c,d] := pop1[t,s,a,c,d] + entrants[c,d] - exits[c,d];
                if pop1[t,s,a,c,d] < 0 then pop1[t,s,a,c,d] := 0;
              end;
            end; {a > 4 }
          end;
        end;
      end;

      CalcNumChildrenEligTreat(t);

      if DP.GetChildTreatInputs(DP_PerChildHIVRecART, t) <> DPNotAvail then
        CalcNewChildARTPatients(t, DP_PerChildHIVRecART)
      else
      begin
        CalcNewChildARTPatients(t, DP_PerChildHIVRecART0_4);
        CalcNewChildARTPatients(t, DP_PerChildHIVRecART5_9);
        CalcNewChildARTPatients(t, DP_PerChildHIVRecART10_14);
      end;

      OnFirstLineArt := 0;
      OnFirstLineArtmale := 0;
      OnFirstLineArtfemale := 0;
      for s := DP_Male to DP_Female do
        for a := 0 to 14 do
          for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
            for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
            begin
              OnFirstLineArt := OnFirstLineArt + pop1[t, s, a, c, d];
              if s = DP_Male then
                 OnFirstLineArtmale := OnFirstLineArtmale
                                     + pop1[t, s, a, c, d]
              else
                OnFirstLineArtfemale := OnFirstLineArtfemale
                                      + pop1[t, s, a, c, d];
      end;
      NeedForFirstLineART := UnmetNeed + OnFirstLineArt;

      { Calculate need for ART as the sum of UnmetNeed and those on ART }
      DP.SetChildARTCalc(DP_NumNeedART, DP_Male, t, UnmetNeedMale + OnFirstLineARTMale);
      DP.SetChildARTCalc(DP_NumNeedART, DP_Female, t, UnmetNeedFemale + OnFirstLineARTFemale);
      DP.SetChildARTCalc(DP_NumNeedART, DP_BothSexes, t, UnmetNeed + OnFirstLineARTMale + OnFirstLineARTFemale);

      (* Calculations for child CTX need and use *)
      CalcNeedForCTX(t);

      { Progress HIV+ children }
      for a := 0 to 14 do
      begin
        a5 := trunc(a/5) + 1;

        { For those not on ART progression already done
          just need to adjust for ART initiation and LTFU}
        for s := DP_Male to DP_Female do
        begin

          { Children who stop ART need to be assigned a transmission timing. We
            assign timings proportional to the distribution in children off ART.
            If everyone is on ART, we split dropouts evenly across timing strata }
          for c := DP_CD4_Start[a] to DP_CD4_Final[a] do
          begin
            pedsum := 0.0;
            for d := DP_P_Perinatal to DP_P_BF12 do
            begin
              ChildTiming[c,d] := pop1[t,s,a,c,d];
              pedsum := pedsum + ChildTiming[c,d];
            end;

            if pedsum > 0.0 then
            begin
              for d := DP_P_Perinatal to DP_P_BF12 do
                ChildTiming[c,d] := ChildTiming[c,d] / pedsum;
            end
            else
            begin
              for d := DP_P_Perinatal to DP_P_BF12 do
                ChildTiming[c,d] := 0.25;
            end;
          end;

          for d := DP_P_Perinatal to DP_P_BF12 do
          begin
            if a < 5 then
            begin
              for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
              begin
                peddrop := DP.GetPercLostFollowUpChild(t) / 100 * (pop1[t,s,a,c,DP_D_ARTlt6m] + pop1[t,s,a,c,DP_D_ARTgt12m]);
                peddrop := peddrop * ChildTiming[c,d];

                pop1[t,s,a,c,d] := pop1[t,s,a,c,d] - NewChildART[a,s,c,d] + peddrop;
                if pop1[t,s,a,c,d] < 0 then pop1[t,s,a,c,d] := 0;
              end;
            end
            else {a > 4}
            begin
              for c := DP_CD4_Ped_GT1000 to DP_CD4_Ped_LT200 do
              begin
                peddrop := DP.GetPercLostFollowUpChild(t) / 100 * (pop1[t,s,a,c,DP_D_ARTlt6m] + pop1[t,s,a,c,DP_D_ARTgt12m]);
                peddrop := peddrop * ChildTiming[c,d];

                pop1[t,s,a,c,d] := pop1[t,s,a,c,d] - NewChildART[a,s,c,d] + peddrop;
                if pop1[t,s,a,c,d] < 0 then pop1[t,s,a,c,d] := 0;
              end;
            end; {a > 4 }
          end;
        end;
        { For those on ART }
        for s := DP_Male to DP_Female do
        begin
          ChildAIDSDeathsOnART[a,s] := 0;
          if a < 5 then
          begin
            for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
            begin
              { Entrants to LT6M category are new ART patients }
              entrants[c,DP_D_ARTlt6m] := NewChildART[a,s,c,DP_P_Perinatal]
                                        + NewChildART[a,s,c,DP_P_BF0]
                                        + NewChildART[a,s,c,DP_P_BF7]
                                        + NewChildART[a,s,c,DP_P_BF12];
              { Entrants are subject to mortality, adjusted for CTX effect }
              v2 := entrants[c,DP_D_ARTlt6m] * (palpha[a,c,s,DP_D_ARTlt6m] + palpha[a,c,s,DP_D_ART6to12M]) / 2;
              v2 := v2 * ( 1 - CTXcoverage * DP.GetEffectTreatChild(DP_Data, DP_ChildEffWithART, 1));
              entrants[c,DP_D_ARTlt6m] := entrants[c,DP_D_ARTlt6m] - v2;
              { Estimate mortality of those already on ART LT6M }
              v1 := pop1[t,s,a,c,DP_D_ARTlt6m] * (palpha[a,c,s,DP_D_ARTlt6m] + palpha[a,c,s,DP_D_ART6to12M]) / 2;
              v1 := v1 * ( 1 - CTXcoverage * DP.GetEffectTreatChild(DP_Data, DP_ChildEffWithART, 1));
              { Entrants to GT12M category are those surviving from LT6M last year }
              entrants[c,DP_D_ARTgt12m] := pop1[t,s,a,c,DP_D_ARTlt6m] - v1;
              { Those on ART LT6 are just the new patients }
              pop1[t,s,a,c,DP_D_ARTlt6m] := entrants[c,DP_D_ARTlt6m];
              { Exits include those dying on ART, adjusted for CTX effects }
              v3 := pop1[t,s,a,c,DP_D_ARTgt12m] * palpha[a,c,s,DP_D_ARTgt12m];
              v3 := v3 * (1 - CTXARTeff * CTXcoverage);
              { Exits also include children who interrupt ART }
              v4 := (pop1[t,s,a,c,DP_D_ARTgt12m] + pop1[t,s,a,c,DP_D_ARTlt6m]) * DP.GetPercLostFollowUpChild(t)/100;
              exits[c,DP_D_ARTgt12m] := v3 + v4;
              pop1[t,s,a,c,DP_D_ARTgt12m] := pop1[t,s,a,c,DP_D_ARTgt12m] + entrants[c,DP_D_ARTgt12m] - exits[c,DP_D_ARTgt12m];

              ChildAIDSDeathsOnART[a,s] := ChildAIDSDeathsOnART[a,s] + v1 + v2 + v3;
            end; {c loop}
          end
          else { a> 5 }
          begin
            for c := DP_CD4_Ped_GT1000 to DP_CD4_Ped_LT200 do
            begin
              entrants[c,DP_D_ARTlt6m] := NewChildART[a,s,c,DP_P_Perinatal]
                                        + NewChildART[a,s,c,DP_P_BF0]
                                        + NewChildART[a,s,c,DP_P_BF7]
                                        + NewChildART[a,s,c,DP_P_BF12];
              v2 := entrants[c,DP_D_ARTlt6m] * (palpha[a,c,s,DP_D_ARTlt6m] + palpha[a,c,s,DP_D_ART6to12M]) / 2;
              v2 := v2 * ( 1 - CTXcoverage * DP.GetEffectTreatChild(DP_Data, DP_ChildEffWithART, 1));
              entrants[c,DP_D_ARTlt6m] := entrants[c,DP_D_ARTlt6m] - v2;

              v1 := pop1[t,s,a,c,DP_D_ARTlt6m] * (palpha[a,c,s,DP_D_ARTlt6m] + palpha[a,c,s,DP_D_ART6to12M]) / 2;
              v1 := v1 * ( 1 - CTXcoverage * DP.GetEffectTreatChild(DP_Data, DP_ChildEffWithART, 1));
              entrants[c,DP_D_ARTgt12m] := pop1[t,s,a,c,DP_D_ARTlt6m] - v1;
              pop1[t,s,a,c,DP_D_ARTlt6m] := entrants[c,DP_D_ARTlt6m];

              v3 := pop1[t,s,a,c,DP_D_ARTgt12m] * palpha[a,c,s,DP_D_ARTgt12m];
              v3 := v3 * (1 - CTXARTeff * CTXcoverage);
              v4 := (pop1[t,s,a,c,DP_D_ARTgt12m] + pop1[t,s,a,c,DP_D_ARTlt6m]) * DP.GetPercLostFollowUpChild(t)/100;
              exits[c,DP_D_ARTgt12m] := v3 + v4;
              pop1[t,s,a,c,DP_D_ARTgt12m] := pop1[t,s,a,c,DP_D_ARTgt12m] + entrants[c,DP_D_ARTgt12m] - exits[c,DP_D_ARTgt12m];

              ChildAIDSDeathsOnART[a,s] := ChildAIDSDeathsOnART[a,s] + v1 + v2 + v3;
            end; {c loop}
          end;
          { Add pediatric deaths on ART to bookkeeping variables }
          DP.SetAIDSDeaths(t, a5, s, DP.GetAIDSDeaths(t,a5,s) + ChildAIDSDeathsOnART[a,s]);
          DP.SetAIDSDeathsARTSingleAge(t, s, a + 1, ChildAIDSDeathsOnART[a, s]);
          AidsDeathsByAge[a,s,t] := AidsDeathsByAge[a,s,t] + ChildAIDSDeathsOnART[a,s];
        end; {s loop}
      end; {a loop}

      DP.SetChildARTCalc(DP_NumOnART, DP_Male, t, 0);
      DP.SetChildARTCalc(DP_NumOnART, DP_Female, t, 0);
      DP.SetChildARTCalc(DP_NumOnART, DP_BothSexes, t, 0);
      for a := 0 to 14 do
        for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
          for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
          begin
            DP.SetChildARTCalc(DP_NumOnART, DP_Male, t,
                        DP.GetChildARTCalc(DP_NumOnART, DP_Male, t)
                        + pop1[t, DP_Male, a, c, d]);
            DP.SetChildARTCalc(DP_NumOnART, DP_Female, t,
                        DP.GetChildARTCalc(DP_NumOnART, DP_Female, t)
                        + pop1[t, DP_Female, a, c, d]);
            DP.SetChildARTCalc(DP_NumOnART, DP_BothSexes, t,
                        DP.GetChildARTCalc(DP_NumOnART, DP_BothSexes, t)
                        + pop1[t, DP_Male, a, c, d]
                        + pop1[t, DP_Female, a, c, d]);
          end;

      (* determine number of children on ART by sex and 5 year age group *)
      for s := DP_Male to DP_Female do
      begin
        SumOnART0_4 := 0;
        SumOnART5_9 := 0;
        SumOnART10_14 := 0;
        for a := DP_A0 to DP_A14 do
        begin
          for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
            for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
            begin
              case a of
                DP_A0..DP_A4   :  SumOnART0_4 := SumOnART0_4 + pop1[t, s, a-1, c, d];
                DP_A5..DP_A9   :  SumOnART5_9 := SumOnART5_9 + pop1[t, s, a-1, c, d];
                DP_A10..DP_A14 :  SumOnART10_14 := SumOnART10_14 + pop1[t, s, a-1, c, d];
              end;
            end;
        end;
        DP.SetNumChildOnARTByAge(t, s, DP_A0_4, SumOnART0_4);
        DP.SetNumChildOnARTByAge(t, s, DP_A5_9, SumOnART5_9);
        DP.SetNumChildOnARTByAge(t, s, DP_A10_14, SumOnART10_14);
      end;
      for a := DP_A0_4 to DP_A10_14 do
        DP.SetNumChildOnARTByAge(t, DP_BothSexes, a, DP.GetNumChildOnARTByAge(t, DP_Male, a) +
                                                       DP.GetNumChildOnARTByAge(t, DP_Female, a));
//DOUBLE CHECK BELOW  -NEW POP1
      (* Determine number of children age 0 on ART by sex *)
      for s := DP_Male to DP_Female do
      begin
        SumOnART_Age0 := 0;
        for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
          for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
          SumOnART_Age0 := SumOnART_Age0 + pop1[t, s, DP_A0-1, c, d];

        DP.SetNumChildOnARTUnder1(t, s, SumOnART_Age0);
      end;
      DP.SetNumChildOnARTUnder1(t, DP_BothSexes, DP.GetNumChildOnARTUnder1(t, DP_Male) +
                                                 DP.GetNumChildOnARTUnder1(t, DP_Female));

      for s := DP_Male to DP_Female do
        for a := 0 to 14 do
        begin
          a5 := a div 5 + 1;
          for c := DP_CD4_Start[a] to DP_CD4_Final[a] do
            for d := DP_P_Perinatal to DP_P_BF12 do
              DP.SetAdultsChildrenStartingART(t, a5, s,
                DP.GetAdultsChildrenStartingART(t, a5, s) + NewChildART[a,s,c,d]);
        end;

      for a5 := DP_A0_4 to DP_A80_up do
        DP.SetAdultsChildrenStartingART(t, a5, DP_BothSexes,
          DP.GetAdultsChildrenStartingART(t, a5, DP_Male)
          + DP.GetAdultsChildrenStartingART(t, a5, DP_Female));

  {
  if t = 6 then
  begin
  AssignFile(OutFile,'Children.txt');
  Rewrite(OutFile);
  writeln(Outfile, 's  a  c  Perinatal, BF0, BF7, BF12');
  for s := DP_Male to DP_Female do
  for a := 0 to 4 do
  for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
  writeln(Outfile, s:2, a:3, c: 2, pop1[a,t,s,c,DP_P_Perinatal]:8:0,pop1[0, t, s, c, DP_P_BF0]:8:0,
                                   pop1[a,t,s,c,DP_P_BF7]:8:0, pop1[a,t,s,c,DP_P_BF12]:8:0);
  writeln(OutFile, 'births ',HIVbirths:4:1, ' BF0-6 ', NewInfFromBFLT6:4:1,' BF6-12 ',NewInfFromBF6to12:4:1,' BF12+ ',
                    NewInfFromBFGT12:4:1,' New Inf ',DP.GetNewInfections(t, DP_A0_4, DP_Male):8:1,
                    DP.GetNewInfections(t, DP_A0_4, DP_Female):8:1);
  CloseFile(OutFile);
  end;
  }
      { Put single year pop data into five year arrays }
      sum1 := 0;
      for s := DP_Male to DP_Female do
      begin
        for a5 := DP_A0_4 to DP_A80_up do
        begin
          DP.SetHIV(t, s, a5, 0);
          DP.SetOnART(t, s, a5,0);
          DP.SetOnSL(t, s, a5, 0);
        end;

        for a := 0 to 80 do
        begin
          a5 := trunc(a / 5) + 1;
          case a of
            0..4:
              begin
                for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
                begin
                  for d := DP_P_Perinatal to DP_P_BF12 do
                  begin
                    DP.SetHIV(t, s, a5, DP.GetHIV(t, s, a5) + pop1[t, s, a, c, d]);
                  end;
                  for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
                  begin
                    DP.SetHIV(t, s, a5, DP.GetHIV(t, s, a5) + pop1[t, s, a, c, d]);
                    DP.SetOnART(t, s, a5, DP.GetOnART(t, s, a5) + pop1[t, s, a, c, d]);
                  end;
                end;
              end;
            5..14:
              begin
                for c := DP_CD4_Ped_GT1000 to DP_CD4_Ped_LT200 do
                begin
                  for d := DP_P_Perinatal to DP_P_BF12 do
                  begin
                    DP.SetHIV(t, s, a5, DP.GetHIV(t, s, a5) + pop1[t, s, a, c, d]);
                  end;
                  for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
                  begin
                    DP.SetHIV(t, s, a5, DP.GetHIV(t, s, a5) + pop1[t, s, a, c, d]);
                    DP.SetOnART(t, s, a5, DP.GetOnART(t, s, a5) + pop1[t, s, a, c, d]);
                  end;
                end;
              end;
            15..80:
              begin
                for c := DP_CD4_GT500 to DP_CD4_LT50 do
                begin
                  DP.SetHIV(t, s, a5, DP.GetHIV(t, s, a5) + pop1[t, s, a, c, DP_D_HIVPos]);
                  for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
                  begin
                    sum1 := sum1 + pop1[t, s, a, c, d];
                    DP.SetHIV(t, s, a5, DP.GetHIV(t, s, a5) + pop1[t, s, a, c, d]);
                    DP.SetOnART(t, s, a5, DP.GetOnART(t, s, a5) + pop1[t, s, a, c, d]);
                  end; {d loop}
                end; {c loop}
              end;
          end;
        end; {a loop}

        { Adjust the 2020 mid-year number on ART for display if monthly figures entered }
        if (GB.GetYear(DP.proj, t) = 2020) and (DP.GetAdultARTByMonth(s, DP_ART_2020, GB_July) > 0) then
        begin
          sum2 := 0;
          for a5 := DP_A15_19 to DP_A80_up do
            sum2 := sum2 + DP.GetOnART(t, s, a5);
          for a5 := DP_A15_19 to DP_A80_up do
            DP.SetOnART(t,s,a5,DP.GetAdultARTByMonth(s, DP_ART_2020, GB_July) * DP.GetOnART(t, s, a5)/sum2);
        end;
      end; {s loop}

      DP.SetARTCalc_15Plus(DP_NumRecART, DP_BothSexes, t, sum1);
      DP.SetARTCalc_15Plus(DP_UnmetNeedART, DP_BothSexes, t,
                                   EligibleAdults[DP_Male,t]
                                   + EligibleAdults[DP_Female,t]);

      for a5 := DP_A0_4 to DP_A80_up do
        for s := DP_Male to DP_Female do
        begin
          DP.SetUnmetNeedART(t, s, a5, DP.GetNeedART(t, s, a5) - DP.GetOnART
              (t, s, a5));
          if GB.ModuleActive(p,GB_RN) then
            if (DP.GetOnART(t-1, s, a5) > 0) then
              DP.SetOnSL(t, s, a5, DP.GetOnSL(t - 1, s, a5) +
                                    (DP.GetOnART(t-1, s, a5) - DP.GetOnSL(t - 1, s, a5)) *
                                     DP_GetRNMigFirstToSecondLine(p, t) / 100
                                     - DP.GetOnSL(t-1, s, a5)
                                     * DP.GetAIDSDeathsART(t-1, a5, s)
                                       / DP.GetOnART(t-1, s, a5));

        end;
      { Sum up for all ages }
      for s := DP_Male to DP_Female do
      begin
        sum := 0;
        sum2 := 0;
        sum3 := 0;
        sum6 := 0;
        sum7 := 0;
        sum8 := 0;
        sum9 := 0;
        for a5 := DP_A0_4 to DP_A80_up do
        begin
          sum := sum + DP.GetHIV(t, s, a5);
          sum2 := sum2 + DP.GetAIDSDeaths(t, a5, s);
          sum3 := sum3 + DP.GetNewInfections(t, a5, s);
          sum6 := sum6 + DP.GetOnART(t, s, a5);
          sum7 := sum7 + DP.GetOnSL(t, s, a5);
          sum8 := sum8 + DP.GetNewlyNeedART(t, s, a5);
          sum9 := sum9 + DP.GetUnmetNeedART(t, s, a5);
        end;
        DP.SetHIV(t, s, DP_AllAges, sum);
        DP.SetAIDSDeaths(t, DP_AllAges, s, sum2);
        DP.SetNewInfections(t, DP_AllAges, s, sum3);
        DP.SetOnART(t, s, DP_AllAges, sum6);
        DP.SetOnSL(t, s, DP_AllAges, sum7);
        DP.SetNewlyNeedART(t, s, DP_AllAges, sum8);
        DP.SetUnmetNeedART(t, s, DP_AllAges, sum9);
      end;

      { Sum up for both sexes }
      for a5 := DP_AllAges to DP_A80_up do
      begin
        DP.SetHIV(t, DP_BothSexes, a5, DP.GetHIV(t, DP_Male, a5) + DP.GetHIV
            (t, DP_Female, a5));
        DP.SetAIDSDeaths(t, a5, DP_BothSexes, DP.GetAIDSDeaths(t, a5, DP_Male)
            + DP.GetAIDSDeaths(t, a5, DP_Female));
        DP.SetAIDSDeathsNoART(t, a5, DP_BothSexes, DP.GetAIDSDeathsNoART
            (t, a5, DP_Male) + DP.GetAIDSDeathsNoART(t, a5, DP_Female));
        DP.SetAIDSDeathsART(t, a5, DP_BothSexes, DP.GetAIDSDeathsART
            (t, a5, DP_Male) + DP.GetAIDSDeathsART(t, a5, DP_Female));
        DP.SetNewInfections(t, a5, DP_BothSexes, DP.GetNewInfections
            (t, a5, DP_Male) + DP.GetNewInfections(t, a5, DP_Female));
        DP.SetOnART(t, DP_BothSexes, a5, DP.GetOnART(t, DP_Male, a5)
            + DP.GetOnART(t, DP_Female, a5));
        DP.SetOnSL(t, DP_BothSexes, a5, DP.GetOnSL(t, DP_Male, a5) + DP.GetOnSL
            (t, DP_Female, a5));
        DP.SetNewlyNeedART(t, DP_BothSexes, a5, DP.GetNewlyNeedART
            (t, DP_Male, a5) + DP.GetNewlyNeedART(t, DP_Female, a5));
        DP.SetUnmetNeedART(t, DP_BothSexes, a5, DP.GetUnmetNeedART
            (t, DP_Male, a5) + DP.GetUnmetNeedART(t, DP_Female, a5));
        DP.SetNewlyNeedART(t, DP_BothSexes, a5, DP.GetNewlyNeedART
            (t, DP_Male, a5) + DP.GetNewlyNeedART(t, DP_Female, a5));
      end; {a5 loop}
//
      // mark0407
      (* Save 'Number of HIV+ children not on ART by time of infection' for later display *)
      for s := DP_Male to DP_Female do
        for a := DP_A0 to DP_A14 do
          for d := DP_P_Perinatal to DP_P_BF12 do
          begin
            value := 0.0;
            for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
              value := value + pop1[t,s,a-1,c,d];
            DP.SetHIVChildNoARTByInfection(t,s,a,d,value);
          end;

      {for s := DP_BothSexes to DP_Female do
        for a := 0 to 14 do
          for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
            for d := DP_P_Perinatal to DP_P_BF12 do
              if s = DP_BothSexes then
                DP.SetHIVChildNoARTByInfection(t,s,a+1,c,d,pop1[a,t,DP_Male,c,d] +
                  pop1[a,t,DP_Female,c,d])
              else
                DP.SetHIVChildNoARTByInfection(t,s,a+1,c,d,pop1[a,t,s,c,d]);}

      { Save HIV pop 15-17 for later display }
      for s := DP_Male to DP_Female do
      begin
        DP.SetHIVPop1517(t, s, 0);
        for a := 15 to 17 do
          for d := DP_D_HIVPos to DP_D_Max do
            for c := DP_CD4_GT500 to DP_H_Max do
              DP.SetHIVPop1517(t, s, DP.GetHIVPop1517(t, s)
                  + pop1[t, s, a, c, d]);
      end; {s loop}
      DP.SetHIVPop1517(t, DP_BothSexes, DP.GetHIVPop1517(t, DP_Male)
          + DP.GetHIVPop1517(t, DP_Female));
      { Save HIV pop under 1 for later display }
      for s := DP_Male to DP_Female do
      begin
        DP.SetHIVPopUnder1(t, s, 0);
        for d := DP_P_Perinatal to DP_D_Max do
          for h := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
            DP.SetHIVPopUnder1(t, s, DP.GetHIVPopUnder1(t, s)
                + pop1[t, s, 0, h, d]);
      end; {s loop}
      DP.SetHIVPopUnder1(t, DP_BothSexes, DP.GetHIVPopUnder1(t, DP_Male)
          + DP.GetHIVPopUnder1(t, DP_Female));

      (* Used to determine Need for ART for age groups 0-4, 5-9 and 10-14 for the Need for ART by Age display *)
      agroup := 1;
      for ag := 0 to 2 do
      begin
        for s := DP_Male to DP_Female do
        begin
          value := 0.0;
          for a := ag + (ag * 4) to ag + (ag * 4) + 4 do // this for statement uses ag to create ranges 0-4, 5-9 and 10-14
            for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
            begin
              for d := DP_P_Perinatal to DP_P_BF12 do
                value := value + elig[a,s,c,d];
              for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
                value := value + pop1[t, s, a, c, d];
            end;
          case ag of
            0: agroup := DP_A0_4;
            1: agroup := DP_A5_9;
            2: agroup := DP_A10_14;
          end;
          DP.SetNeedART(t, s, agroup, value);
        end; {s loop}
        DP.SetNeedART(t, DP_BothSexes, agroup, DP.GetNeedART(t, DP_Male, agroup)
                                                 + DP.GetNeedART(t, DP_Female, agroup));
        DP.SetNeedART(t, DP_Male, DP_ALLAGES, DP.GetNeedART(t, DP_Male, DP_ALLAGES)
                      +DP.GetNeedART(t, DP_Male, agroup));
        DP.SetNeedART(t, DP_Female, DP_ALLAGES, DP.GetNeedART(t, DP_Female, DP_ALLAGES)
                      +DP.GetNeedART(t, DP_Female, agroup));
        DP.SetNeedART(t, DP_BothSexes, DP_ALLAGES, DP.GetNeedART(t, DP_BothSexes, DP_ALLAGES)
                      +DP.GetNeedART(t, DP_BothSexes, agroup));

      end; {ag loop}
      CalcNeedForARTEndYear;

      {                                                  }
      { Calculate median CD4 count of those starting ART }
      {                                                  }

      { Sum to get totals by CD4 category for all ages and both sexes }
      for c := DP_CD4_GT500 to DP_CD4_LT50 do
      begin
        NewARTCD4Distr[c] := 0;
        for a := 15 to DP_MaxSingleAges do
          for s := DP_Male to DP_Female do
            NewARTCD4Distr[c] := NewARTCD4Distr[c] + pop1[t,s,a,c,DP_D_ARTlt6m];
      end;

      { Find total new ART patients }
      sum := 0;
      for c := DP_CD4_GT500 to DP_CD4_LT50 do
        sum := sum + NewARTCD4Distr[c];

      { Calculate proportions }
      if sum >0 then
      begin
        for c := DP_CD4_GT500 to DP_CD4_LT50 do
          Prop3[c] := NewARTCD4Distr[c] / sum;
         { Find category with total just above half the total }
         sum1 := 0;
         c := DP_CD4_LT50 + 1;
         while sum1 < 0.5 do
         begin
           c := c - 1;
           sum1 := sum1 + Prop3[c];
         end;
         { Subtract last group, which exceeded the median }
         sum1 := sum1 - Prop3[c];
         MedCD4 := DP_CD4_lower_limit[c]
                   + (DP_CD4_lower_limit[c-1] - DP_CD4_lower_limit[c])
                   * (0.5 - sum1) / Prop3[c];
      end else MedCD4 := 0;
      DP.SetMedianCD4(t,MedCD4);

      //information needed for TIME modules
      if GB.ModuleActive(p, GB_TB) then
      begin

        (*EPW Start*)
        //Record HIV population, on ART, by cd4, for TB Impact model
        //CD4 distribution for children 5 years old. Use as average CD4 distribution for children.
        Age0to14CD4Distr[DP_CD4_GT500] := 0.58;
        Age0to14CD4Distr[DP_CD4_350_500] := 0.13;
        Age0to14CD4Distr[DP_CD4_250_349] := 0.08;
        Age0to14CD4Distr[DP_CD4_200_249] := 0.04;
        Age0to14CD4Distr[DP_CD4_100_199] := 0.09;
        Age0to14CD4Distr[DP_CD4_50_99] := 0.04;
        Age0to14CD4Distr[DP_CD4_LT50] := 0.04;

        t1:=10;//at the moment, record only the final t,t1 population structure

        for s := DP_Male to DP_Female do
        begin

          //Children
          for a := 0 to 14 do
          begin
            {Calculate index of 5 year age groups }
            a5 := trunc(a / 5) + 1;

            //HIV-
            DP_TI_SetHIVbyAge(p, t, t1, s, a5, DP_H_No, DP_HIVNeg,
              DP_TI_GetHIVbyAge(p,t,t1,s,a5,DP_H_No,DP_HIVNeg)+pop1[t, s, a, DP_H_No, DP_D_HIVNeg]);

            //0-4, no ART and on ART
            if(a5=DP_A0_4) then
            begin

              //HIV+, no ART
              NotOnART:=0;
              for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
              for d := DP_P_Perinatal to DP_P_BF12 do
                NotOnART:=NotOnART+pop1[t, s, a, c, d];

              for c := DP_CD4_GT500 to DP_CD4_LT50 do
                DP_TI_SetHIVbyAge(p, t, t1, s, a5, c, DP_HIVPos,
                  DP_TI_GetHIVbyAge(p,t,t1,s,a5,c,DP_HIVPos)+Age0to14CD4Distr[c]*NotOnART);

              //HIV+, on ART
              for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
              begin

                OnART:=0;
                for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
                  OnART:=OnART+pop1[t, s, a, c, d];

                d1 := d + DP_HIVARTlt6m - DP_D_ARTlt6m; // map d from AIM to TB model index values
                for c := DP_CD4_GT500 to DP_CD4_LT50 do
                  DP_TI_SetHIVbyAge(p, t, t1, s, a5, c, d1,
                    DP_TI_GetHIVbyAge(p, t,t1,s,a5,c,d1)+Age0to14CD4Distr[c]*OnART);
              end;
            end;{a5=DP_A0_4}

            //5-14, HIV+, no ART and on ART
            if(a5 in [DP_A5_9,DP_A10_14]) then
            for c := DP_CD4_Ped_GT1000 to DP_CD4_Ped_LT200  do
            begin
               case c of
                 DP_CD4_Ped_GT1000..DP_CD4_Ped_500_749:
                 begin
                    //no ART
                    DP_TI_SetHIVbyAge(p, t, t1, s, a5, DP_CD4_GT500, DP_HIVPos,
                      DP_TI_GetHIVbyAge(p, t,t1,s,a5,DP_CD4_GT500,DP_HIVPos)+
                      pop1[t, s, a, c, DP_P_Perinatal]+
                      pop1[t, s, a, c, DP_P_BF0]+
                      pop1[t, s, a, c, DP_P_BF7]+
                      pop1[t, s, a, c, DP_P_BF12]);

                    //ART
                    for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
                    begin
                      d1 := d + DP_HIVARTlt6m - DP_D_ARTlt6m; // map d from AIM to TB model index values
                      DP_TI_SetHIVbyAge(p, t, t1, s, a5, DP_CD4_GT500, d1,
                        DP_TI_GetHIVbyAge(p, t,t1,s,a5,DP_CD4_GT500,d1) + pop1[t, s, a, c, d]);
                    end;
                 end;

                 DP_CD4_Ped_350_499:
                 begin
                    //no ART
                    DP_TI_SetHIVbyAge(p, t, t1, s, a5, DP_CD4_350_500, DP_HIVPos,
                      DP_TI_GetHIVbyAge(p, t,t1,s,a5,DP_CD4_350_500,DP_HIVPos)+
                      pop1[t, s, a, c, DP_P_Perinatal]+
                      pop1[t, s, a, c, DP_P_BF0]+
                      pop1[t, s, a, c, DP_P_BF7]+
                      pop1[t, s, a, c, DP_P_BF12]);

                    //ART
                    for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
                    begin
                      d1 := d + DP_HIVARTlt6m - DP_D_ARTlt6m; // map d from AIM to TB model index values
                      DP_TI_SetHIVbyAge(p, t, t1, s, a5, DP_CD4_350_500, d1,
                        DP_TI_GetHIVbyAge(p, t,t1,s,a5,DP_CD4_350_500,d1) + pop1[t, s, a, c, d]);
                    end;
                 end;

                 DP_CD4_Ped_200_349:
                 begin
                      DP_TI_SetHIVbyAge(p, t, t1, s, a5, DP_CD4_250_349, DP_HIVPos,
                        DP_TI_GetHIVbyAge(p, t,t1,s,a5,DP_CD4_250_349,DP_HIVPos)+
                        2/3*(pop1[t, s, a, c, DP_P_Perinatal]+
                        pop1[t, s, a, c, DP_P_BF0]+
                        pop1[t, s, a, c, DP_P_BF7]+
                        pop1[t, s, a, c, DP_P_BF12]));

                      for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
                      begin
                        d1 := d + DP_HIVARTlt6m - DP_D_ARTlt6m; // map d from AIM to TB model index values
                        DP_TI_SetHIVbyAge(p, t, t1, s, a5, DP_CD4_250_349, d1,
                          DP_TI_GetHIVbyAge(p, t,t1,s,a5,DP_CD4_250_349,d1) + 2/3*pop1[t, s, a, c, d]);
                      end;

                      DP_TI_SetHIVbyAge(p, t, t1, s, a5, DP_CD4_200_249, DP_HIVPos,
                        DP_TI_GetHIVbyAge(p, t,t1,s,a5,DP_CD4_200_249,DP_HIVPos)+
                        1/3*(pop1[t, s, a, c, DP_P_Perinatal]+
                        pop1[t, s, a, c, DP_P_BF0]+
                        pop1[t, s, a, c, DP_P_BF7]+
                        pop1[t, s, a, c, DP_P_BF12]));

                      for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
                      begin
                        d1 := d + DP_HIVARTlt6m - DP_D_ARTlt6m;
                        DP_TI_SetHIVbyAge(p, t, t1, s, a5, DP_CD4_200_249, d1,
                          DP_TI_GetHIVbyAge(p, t,t1,s,a5,DP_CD4_200_249,d1) + 1/3*pop1[t, s, a, c, d]);
                      end;
                 end;

                 DP_CD4_Ped_LT200:
                 begin
                    //no ART
                    DP_TI_SetHIVbyAge(p, t, t1, s, a5, DP_CD4_100_199, DP_HIVPos,
                      DP_TI_GetHIVbyAge(p, t,t1,s,a5,DP_CD4_100_199,DP_HIVPos)+
                      0.35*(pop1[t, s, a, c, DP_P_Perinatal]+
                      pop1[t, s, a, c, DP_P_BF0]+
                      pop1[t, s, a, c, DP_P_BF7]+
                      pop1[t, s, a, c, DP_P_BF12]));
                    //ART
                    for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
                    begin
                      d1 := d + DP_HIVARTlt6m - DP_D_ARTlt6m;
                      DP_TI_SetHIVbyAge(p, t, t1, s, a5, DP_CD4_100_199, d1,
                        DP_TI_GetHIVbyAge(p, t,t1,s,a5,DP_CD4_100_199,d1) + 0.35*pop1[t, s, a, c, d]);
                    end;

                    //no ART
                    DP.SetHIVbyAge(t, t1, s, a5, DP_CD4_50_99, DP_HIVPos,
                      DP.GetHIVbyAge(t,t1,s,a5,DP_CD4_50_99,DP_HIVPos)+
                      0.21*(pop1[t, s, a, c, DP_P_Perinatal]+
                      pop1[t, s, a, c, DP_P_BF0]+
                      pop1[t, s, a, c, DP_P_BF7]+
                      pop1[t, s, a, c, DP_P_BF12]));
                    //ART
                    for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
                    begin
                      d1 := d + DP_HIVARTlt6m - DP_D_ARTlt6m;
                      DP_TI_SetHIVbyAge(p, t, t1, s, a5, DP_CD4_100_199, d1,
                       DP_TI_GetHIVbyAge(p, t,t1,s,a5,DP_CD4_100_199,d1) + 0.21*pop1[t, s, a, c, d]);
                    end;

                    //no ART
                    DP_TI_SetHIVbyAge(p, t, t1, s, a5, DP_CD4_50_99, DP_HIVPos,
                      DP_TI_GetHIVbyAge(p, t, t1,s,a5,DP_CD4_50_99,DP_HIVPos)+
                      0.44*(pop1[t, s, a, c, DP_P_Perinatal]+
                      pop1[t, s, a, c, DP_P_BF0]+
                      pop1[t, s, a, c, DP_P_BF7]+
                      pop1[t, s, a, c, DP_P_BF12]));
                    //ART
                    for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
                    begin
                      d1 := d + DP_HIVARTlt6m - DP_D_ARTlt6m;
                      DP_TI_SetHIVbyAge(p, t, t1, s, a5, DP_CD4_LT50, d1,
                         DP_TI_GetHIVbyAge(p, t,t1,s,a5,DP_CD4_LT50,d1) + 0.44*pop1[t, s, a, c, d]);
                    end;
                 end;
               end;{for c}
            end;{a5 in [DP_A5_9,DP_A10_14])}
          end;//children a loop

          //Adults
          for a := 15 to DP_MaxSingleAges do
          begin
            {Calculate index of 5 year age groups }
            a5 := trunc(a / 5) + 1;
            //HIV-
            DP_TI_SetHIVbyAge(p, t, t1, s, a5, DP_H_No, DP_HIVNeg,
              DP_TI_GetHIVbyAge(p, t,t1,s,a5,DP_H_No,DP_HIVNeg)+pop1[t, s, a, DP_H_No, DP_D_All]);

            //HIV+ not at ART
            d:=DP_D_HIVPos;
            for c := DP_CD4_GT500 to DP_CD4_LT50 do
              DP_TI_SetHIVbyAge(p, t, t1, s, a5, c, DP_HIVPos,
                DP_TI_GetHIVbyAge(p, t,t1,s,a5,c,DP_HIVPos)+pop1[t, s, a, c, d]);

            //HIV+ on ART
            for c := DP_CD4_GT500 to DP_CD4_LT50 do
            for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
            begin
              d1 := d + DP_HIVARTlt6m - DP_D_ARTlt6m; // map d from AIM to TB model index values
              DP_TI_SetHIVbyAge(p, t, t1, s, a5, c, d1,
                DP_TI_GetHIVbyAge(p, t,t1,s,a5,c,d1)+pop1[t, s, a, c, d]);
            end;
          end;//adult a loop
        end;//s loop

        //sum for both sexes
        for a := DP_A0_4 to DP_MAX_AGE do
        for c := DP_H_No to DP_CD4_LT50 do
        for d := DP_HIVNeg to DP_HIVARTgt12m do
          DP_TI_SetHIVbyAge(p, t, t1, DP_BothSexes, a, c, d,
            DP_TI_GetHIVbyAge(p, t, t1, DP_Male, a, c, d)
            +DP_TI_GetHIVbyAge(p, t, t1, DP_Female, a, c, d));

        //set aging rate
        for a := DP_A0_4 to (DP_MAX_AGE-1) do
        for s := DP_Male to DP_Female do
        begin
              //HIV-
              a5:=a*5-1;
              sum:=DP_TI_GetHIVbyAge(p, t, t1, s, a, DP_H_No, DP_HIVNeg);
              sum1:=pop1[t, s, a5, DP_H_No, DP_D_HIVNeg];

              if(sum>0) then
                value:=sum1/sum
              else
                value:=0;

              DP.SetHIVAgingRate(t, s, a, DP_HIVNeg,value);

              //HIV+
              sum:=0;
              for c := DP_CD4_GT500 to DP_CD4_LT50 do
                sum:=sum+DP_TI_GetHIVbyAge(p,t,t1,s,a,c,DP_HIVPos);

              sum1:=0;
              for c := DP_CD4_GT500 to DP_CD4_LT50 do
              if(a<=DP_A10_14) then
              begin
                for d := DP_P_Perinatal to DP_P_BF12 do
                sum1:=sum1+pop1[t, s, a5, c, d];
              end
              else
              begin
                d:=DP_D_HIVPos;
                sum1:=sum1+pop1[t, s, a5, c, d];
              end;

              if(sum>0) then
                value:=sum1/sum
               else
                value:=0;
              DP.SetHIVAgingRate(t, s, a, DP_HIVPos,value);

             //HIV+ on ART
             d:=DP_D_ARTgt12m;
             sum:=0;
             for c := DP_CD4_GT500 to DP_CD4_LT50 do
                sum:=sum+DP_TI_GetHIVbyAge(p,t,t1,s,a,c,DP_HIVARTgt12m);

             sum1:=0;
             for c := DP_CD4_GT500 to DP_CD4_LT50 do
                sum1:=sum1+pop1[t, s, a5, c, d];

             if(sum>0) then
               value:=sum1/sum
             else
               value:=0;

             DP.SetHIVAgingRate(t, s, a, DP_HIVARTlt6m,value);
             DP.SetHIVAgingRate(t, s, a, DP_HIVARTgt6m,value);
             DP.SetHIVAgingRate(t, s, a, DP_HIVARTgt12m,value);

        end;{a,s loop}


        (*Carel: ART info needed by TB module*)
        temp := 0;
        temp_all := 0;
        for c := DP_CD4_GT500 to DP_CD4_LT50 do
        for a := 15 to DP_MaxSingleAges do
        begin
          temp := temp + pop1[t,DP_Male,a,c,DP_D_ARTgt12M] + pop1[t,DP_Female,a,c,DP_D_ARTgt12M];
          for d := DP_D_ARTlt6M to DP_D_ARTgt12M do
            temp_all := temp_all + pop1[t,DP_Male,a,c,d] + pop1[t,DP_Female,a,c,d];
        end;{for c, a}

        if(temp_all > 0 ) then
          temp := temp / temp_all * 100
        else
          temp := 0;
        DP.SetOnARTgt12M(t, temp);

      end  //if GB.ModuleActive(p, GB_TB) then

    (*EPW End*)
    end; {t loop}
    for s := DP_BothSexes to DP_Female do
      finalize(EligibleAdults[s]);
  end;
{$ENDREGION}

{$REGION 'DPNewHIV'}
  procedure TDPProj.CalcNewAdultInfections(t : integer);
  var
    s, a, c, d, a5 : byte;
    sum1           : double;
  begin
    { Calculate new adult infections }

    (* Sum Array to get all durations and status *)
    SumPopulationByHIVStatus(t);

    (* If the projection has an SPT file, then set the last age to the AGERANGE
      from the file. Otherwise, set the 49. *)
    if DP.GetEppAgeRange > 0 then
      LastAge := DP.GetEppAgeRange
    else
      LastAge := 49;

    {Split incidence age pattern into single ages }
    for s := DP_Male to DP_Female do
    begin
      for a := 0 to DP_MaxSingleAges do
        sa[a] := 0;
      for a := DP_A0_4 to DP_MAX_AGE do
        p5[a] := DP.GetDistOfHIV(t,s,a);

      Calc_Single_Ages(GB,DP,p5, sa, s);
      for a := 0 to DP_MaxSingleAges do
      begin
        PrevRatio[a,s] := sa[a];
        if PrevRatio[a,s] < 0 then
        begin
          PrevRatio[a,s] := 0;
          if (a>14) and (a<60) then
            DP_AddLogEntry(GB, DP_DPAIM,'DPNewHIV','~2-10','PrevRatio[a,s]', DP_LTZ, PrevRatio[a,s],
                           'a = ' + IntToStr(a) + ',s = ' + IntToStr(s));
        end;
      end;
    end;

    DeathsAverted := 0;

    if (t = 2) and (not GB.GetInCalcStateMode(p)) then
      IncidenceAdjustmentFactor := 1;

    { If the EPP prevalence adjustment is turned on - do the adjustment }
    if DP.GetEPPPrevAdj then
    begin

      { If earlier than 2014 determine new infections to match prevalence }
      if DP.GetFirstYear + t - 1 <= DP.GetEPPPrevAdjustmentYear{2010} then
      begin
        { Calculate number of adults 15-49 }
        adultMales := 0;
        adultFemales := 0;
        for a := 15 to LastAge do
          for s := DP_Male to DP_Female do
            if s = DP_Male then adultMales := adultMales + pop1[t,s,a,DP_H_All,DP_D_All]
                           else adultFemales := adultFemales + pop1[t,s,a,DP_H_All,DP_D_All];
        adults := adultMales + adultFemales;

        { Calculate number of HIV+ adults }
        for s := DP_Male to DP_Female do
        begin
          adultHIV[s] := 0;
          for a := 15 to LastAge do
            for c := DP_CD4_GT500 to DP_H_Max do
              for d := DP_D_HIVPos to DP_D_Max do
                adultHIV[s] := adultHIV[s] + pop1[t,s,a,c,d];
        end;
        adultHIV[DP_BothSexes] := adultHIV[DP_male] + adultHIV[DP_female];

        { Determine total HIV infections required to attain input HIV prevalence }
        NeededHIV[DP_BothSexes] := adults * DP.GetAdultPrevalence(t) / 100
                                 - adultHIV[DP_BothSexes];
        if (DP.GetIncidenceByFit(DP.GetIncidenceOptions, t) <= 0) or (NeededHIV[DP_BothSexes] <= 0) then
        begin
          IncidenceAdjustmentFactor := 1;
          if (NeededHIV[DP_BothSexes] <= 0) then IncidenceAdjustmentFactor := 0;
          if DP.GetIncidenceByFit(DP.GetIncidenceOptions, t) < 0 then
            DP_AddLogEntry(GB, DP_DPAIM,'DPNewHIV','~2-70','DP.GetIncidenceByFit(DP.GetIncidenceOptions, t)',DP_LTEZ,
                            DP.GetIncidenceByFit(DP.GetIncidenceOptions, t),'t = ' + IntToStr(t));
          if NeededHIV[DP_BothSexes] < 0 then
            DP_AddLogEntry(GB, DP_DPAIM,'DPNewHIV','~2-80','NeededHIV[DP_BothSexes]',DP_LTEZ,
                            NeededHIV[DP_BothSexes],'t = ' + IntToStr(t));
        end
        else
          IncidenceAdjustmentFactor := (NeededHIV[DP_BothSexes]
                                      / (adults - adultHIV[DP_BothSexes]) * 100)
                                      / DP.GetIncidenceByFit(DP.GetIncidenceOptions, t);

        (* Save the IncidenceAdjustmentFactor for the year 2010 so it can be used
           for years greater than 2010 when running calculations from the calc year
           rather than 1970 *)
        if ((DP.GetFirstYear + t - 1) = DP.GetEPPPrevAdjustmentYear{2010}) and DP.GetEPPPrevAdj then
          DP.SetIncidenceAdjustmentFactor2012(IncidenceAdjustmentFactor)
        else DP.SetIncidenceAdjustmentFactor2012(1.0);
      end
      else
      begin
        (* Use the IncidenceAdjustmentFactor from 2010 for all years past 2010 *)
        IncidenceAdjustmentFactor := DP.GetIncidenceAdjustmentFactor2012;
      end;

      if IncidenceAdjustmentFactor > DP.GetEPPMaxAdjFactor then
        IncidenceAdjustmentFactor := DP.GetEPPMaxAdjFactor;

      if IncidenceAdjustmentFactor < 1 / DP.GetEPPMaxAdjFactor then
        IncidenceAdjustmentFactor := 1 / DP.GetEPPMaxAdjFactor;

    end
    else IncidenceAdjustmentFactor := 1;
    DP.SetIncidenceAdjustmentFactor(t,IncidenceAdjustmentFactor);

    { Determine total infections needed to match input incidence }

    { Calculate number of males and females in previous year }
    adultMales := 0;
    adultFemales := 0;
    for a := 15 to LastAge do
      for s := DP_Male to DP_Female do
        if s = DP_Male then adultMales := adultMales + pop1[t-1,s,a,DP_H_All,DP_D_All]
                       else adultFemales := adultFemales + pop1[t-1,s,a,DP_H_All,DP_D_All];

    { Calculate number of HIV+ males and females in previous year }
    for s := DP_Male to DP_Female do
    begin
      adultHIV[s] := 0;
      for a := 15 to LastAge do
      begin
        for c := DP_CD4_GT500 to DP_H_Max do
          for d := DP_D_HIVPos to DP_D_Max do
            adultHIV[s] := adultHIV[s] + pop1[t-1,s,a,c,d];
      end;
    end;

    { Calculate female incidence from total incidence and sex ratio }
    if (DP.GetHIVSexRatio({DP_Data, DP.GetIncEpidemicRGIdx, }t) = 0) then
    begin
      temp := 0;
      if DP.GetHIVSexRatio({DP_Data, DP.GetIncEpidemicRGIdx, }t) = 0 then
        DP_AddLogEntry(GB, DP_DPAIM,'DPNewHIV','~2-110','DP.GetHIVSexRatio(t)',
                        DP_EZ,DP.GetHIVSexRatio({DP_Data, DP.GetIncEpidemicRGIdx, }t),'t = ' + IntToStr(t));
    end
    else if (adultFemales - adultHIV[DP_female] + (adultMales - adultHIV[DP_male]) / DP.GetHIVSexRatio({DP_Data, DP.GetIncEpidemicRGIdx, }t)) = 0  then
    begin
      temp := 0;
      DP_AddLogEntry(GB, DP_DPAIM,'DPNewHIV','~2-100','(adultFemales - adultHIV[DP_female]) + (adultMales - adultHIV[DP_male])',
                      DP_EZ,(adultFemales - adultHIV[DP_female]) + (adultMales - adultHIV[DP_male]),' ');

    end
    else
    begin //getting new infections directly from Goals
      if (GB.ModuleActive(p, GB_HV)) then
      begin
        temp := 0*DP_GetHVTotalNewInfection(p,t) / (adultFemales + adultMales - adultHIV[DP_female] - adultHIV[DP_male]);
        (* Set sex ratio from Goals *)
        DP.SetHIVSexRatio({DP_Data, DP.GetIncEpidemicRGIdx, }t, DP_GetHVIncSexRatio(p,t));
        (* Turn off prevalence adjustment *)
        DP.SetEPPPrevAdj(False);
      end
      else
        temp:= DP.GetIncidenceByFit(DP.GetIncidenceOptions, t) / 100;

      temp := temp * IncidenceAdjustmentFactor
            * (adultFemales - adultHIV[DP_female] + adultMales - adultHIV[DP_male])
            / (adultFemales - adultHIV[DP_female] + (adultMales - adultHIV[DP_male]) / DP.GetHIVSexRatio({DP_Data, DP.GetIncEpidemicRGIdx, }t));
    end;
    NewHIV[DP_Female] := temp * (adultFemales - adultHIV[DP_Female]);
    NeededHIV[DP_Female] := temp
                            * (adultFemales - adultHIV[DP_Female])
                            + adultHIV[DP_Female];

    { Calculate male incidence from female incidence and prevalence ratio }
    if DP.GetHIVSexRatio({DP_Data, DP.GetIncEpidemicRGIdx, }t) = 0 then
    begin
    if (GB.ModuleActive(p, GB_HV)) then
      temp := 0*DP_GetHVTotalNewInfection(p,t) / (adultFemales + adultMales - adultHIV[DP_female] - adultHIV[DP_male])
    else
      temp:= DP.GetIncidenceByFit(DP.GetIncidenceOptions, t) / 100;

      temp := temp * IncidenceAdjustmentFactor
              * (adultFemales - adultHIV[DP_female] + adultMales - adultHIV[DP_male])
              / (adultMales - adultHIV[DP_male]);

      DP_AddLogEntry(GB, DP_DPAIM,'DPNewHIV','~2-120','DP.GetHIVSexRatio(t)',
                       DP_EZ,DP.GetHIVSexRatio({DP_Data, DP.GetIncEpidemicRGIdx, }t),'t = ' + IntToStr(t));
    end
    else
      temp := temp / DP.GetHIVSexRatio({DP_Data, DP.GetIncEpidemicRGIdx, }t);

    NewHIV[DP_Male] := temp * (adultMales - adultHIV[DP_Male]);
    NeededHIV[DP_Male] := temp
                            * (adultMales - adultHIV[DP_Male])
                            + adultHIV[DP_Male];
    NeededHIV[DP_BothSexes] := NeededHIV[DP_Female] + NeededHIV[DP_Male];

    NewHIV[DP_BothSexes] := NewHIV[DP_Male] + NewHIV[DP_Female];
    { Distribute new HIV infections by age and sex }

    for s := DP_Male to DP_Female do
    begin
      for a := DP_A0_4 to DP_A80_up do
          DP.SetNewInfections(t,a,s,0);

      for a := 0 to 14 do
        NewHIVbyAge[a] := 0;

       { Calculate new infections if incidence were 1 in 25-29 age group }
      sum1 := 0;
      for a := 15 to LastAge do
        sum1 := sum1 + pop1[t,s,a,DP_H_No,DP_D_HIVNeg] * PrevRatio[a,s];
      { Scale new infections to match total required }
      if sum1 >0 then AdjFactor := NewHIV[s] / sum1 else AdjFactor := 1;

      for a := 15 to DP_MaxSingleAges do
        NewHIVbyAge[a] := AdjFactor * PrevRatio[a,s] * pop1[t,s,a,DP_H_No,DP_D_HIVNeg];
      for a := 15 to DP_MaxSingleAges do
      begin
        case a of
          15 .. 24               : a10 := DP_CD4_15_24;
          25 .. 34               : a10 := DP_CD4_25_34;
          35 .. 44               : a10 := DP_CD4_35_44;
          45 .. DP_MaxSingleAges : a10 := DP_CD4_45_54;
        end;
  {          NewHIVbyAge[a] := NewHIV[s] * PrevRatio[a,s]; }
        a5 := trunc(a/5) + 1;
        (*When Goals is active NewHIVbyAge will be 0 here. NewHIVbyAge is set in InnerLoopCalcs for Goals*)
        if not (GB.ModuleActive(p, GB_HV)) then
        begin
          DP.SetNewInfections(t, a5, s, DP.GetNewInfections(t,a5,s) + NewHIVbyAge[a]);
          DP.SetNewInfectionsBySingleAge(t, s, a+1, NewHIVbyAge[a]); // DP_A<a> = a + 1
        end;
        for c := DP_CD4_GT500 to DP_CD4_LT50 do
        pop1[t,s,a,c,DP_D_HIVPos] := pop1[t,s,a,c,DP_D_HIVPos]
                           + NewHIVbyAge[a]
                           * DP.GetAdultDistNewInfectionsCD4(DP_Data,s,a10,c) / 100;
        pop1[t,s,a,DP_H_No,DP_D_HIVNeg] := pop1[t,s,a,DP_H_No,DP_D_HIVNeg] - NewHIVbyAge[a];
      end;
    end; {s loop}
  end;
{$ENDREGION}

{$REGION 'DPHIVbirths'}
  procedure TDPProj.CalcHIVBirths(t : integer);
  var
    a, d, c, a5,
    TreatOpt         : byte;
    birthsCurrAge,
//    MedianDurationBF,
    numPMTCT, sumARV, need,
    sum1, prev, df, IncRate   : double;
  begin
    { Calculate HIV+ births }
    { Calculate births to HIV+ mothers }
    births := 0;
    births1524 := 0;
    if DP.GetASFR(t, DP_AllAges) <> 0 then
    begin


      for a := 15 to 49 do
      begin
        a5 := trunc(a/5)+1;

        if pop1[t,DP_Female,a,DP_H_All, DP_D_All] > 0 then
          prev := (pop1[t,DP_Female,a,DP_H_All, DP_D_All]
                   - pop1[t,DP_Female,a,DP_H_No, DP_D_HIVNeg])
                  / pop1[t,DP_Female,a,DP_H_All, DP_D_All]
        else
          prev := 0;


        { df = discount factor }
        df := 0;
        { For HIV+ women not on ART }
        for c := DP_CD4_GT500 to DP_CD4_LT50 do
          df := df + DP.GetFRRByLocation(DP_Data)
                     * DP.GetFertCD4Discount(DP_Data, c)
                     * DP.GetHIVTFR(DP_Data,a5,t)
                     * (pop1[t,DP_Female,a,c,DP_D_HIVPos] + pop1[t-1,DP_Female,a,c,DP_D_HIVPos])/2 ;
        { For women on ART less than 6 months use same discount }
        for c := DP_CD4_GT500 to DP_CD4_LT50 do
          df := df + DP.GetFRRByLocation(DP_Data)
                   * DP.GetFertCD4Discount(DP_Data, c)
                   * DP.GetHIVTFR(DP_Data,a5,t)
                   * (pop1[t,DP_Female,a,c,DP_D_ARTlt6m] + pop1[t-1,DP_Female,a,c,DP_D_ARTlt6m])/2;
        { For HIV+ women on ART for more than 6 months there is no discount }
        for c := DP_CD4_GT500 to DP_CD4_LT50 do
          for d := DP_D_ART6to12M to DP_D_ARTgt12m do
              df := df + DP.GetFRRByLocation(DP_Data)
                       * DP.GetRatioWomenOnART(DP_Data, a5)
                       * (pop1[t,DP_Female,a,c,d] + pop1[t-1,DP_Female,a,c,d])/2;

        { Calculate weighted average }
        if (pop1[t,DP_Female,a,DP_H_All,DP_D_All]
            -pop1[t,DP_Female,a,DP_H_No,DP_D_HIVNeg]) > 0
            then df := df / ((pop1[t,DP_Female,a,DP_H_All,DP_D_All]
                              - pop1[t,DP_Female,a,DP_H_No,DP_D_HIVNeg]
                              + pop1[t-1,DP_Female,a,DP_H_All,DP_D_All]
                              - pop1[t-1,DP_Female,a,DP_H_No,DP_D_HIVNeg])/2)
            else df := 1;

        if DP.GetASFR(t, DP_AllAges) > 0 then
        begin
          birthsCurrAge := (pop1[t,DP_Female,a,DP_H_All,DP_D_All]
                            - pop1[t,DP_Female,a,DP_H_No,DP_D_HIVNeg]
                            + pop1[t-1,DP_Female,a,DP_H_All,DP_D_All]
                            - pop1[t-1,DP_Female,a,DP_H_No,DP_D_HIVNeg]) / 2
                         * DP.GetTFR(t)
                         * df / (df *prev + (1 - prev))
                         * DP.GetASFR(t, a5) / DP.GetASFR(t, DP_AllAges) / 5;
          births := births + birthsCurrAge;
          if a < 25 then
            births1524 := births1524 + birthsCurrAge;
        end;
      end; { a loop }
    end else births := 0;

    { Adjust births to HIV+ mothers for abortion }
    if DP.GetPregTermAbortionPerNum(t) = DP_Number then
      births := births - DP.GetPregTermAbortion(t)
    else births := births * (1 - DP.GetPregTermAbortion(t) / 100);
    if births < 0 then births := 0;
    DP.SetNeedPMTCT(t, births);
    DP.SetHIVPregWomen(t,births1524);

    { Mother to child transmission rate is weighted average of
      those receiving PMTCT and those not treated }
    { Calculate the percent of women by treatment option }
    TreatPercent[DP_NoProphylaxis] := 1.0;
    for TreatOpt := DP_SingleDoseNevir to DP_TotalTreat do
      TreatPercent[TreatOpt] := 0;
    { Entries might be number or percent - check number first }
    { If sumARV is zero then entries must be percent }
    sumARV := 0;
    for TreatOpt := DP_SingleDoseNevir to DP_TripleARTDurPreg_Late do
      sumARV := sumARV + DP.GetARVRegimen(DP_PrenatalProphylaxis, TreatOpt, DP_Number, t);
    if (sumARV <> 0) and (DP.GetNeedPMTCT(t) <> 0) then
    begin { PMTCT data entered as numbers }
      numPMTCT := max(0, sumARV + DP.GetPatientsReallocated(t));
      DP.SetOnPMTCT(t, numPMTCT);
      DP.SetPMTCTEffReg(t, numPMTCT * (1.0 - DP.GetARVRegimen(DP_PrenatalProphylaxis, DP_SingleDoseNevir, DP_Number, t) / sumARV));

      { Calculate the percentage of HIV+ pregnant women on each regimen }
      need := max(DP.GetNeedPMTCT(t), numPMTCT);
      for TreatOpt := DP_SingleDoseNevir to DP_TripleARTDurPreg_Late do
        TreatPercent[TreatOpt] := (DP.GetARVRegimen(DP_PrenatalProphylaxis, TreatOpt, DP_Number, t) / sumARV) * (numPMTCT / need);
    end
    else
    begin { PMTCT data entered as percentage }
      sumARV := 0;
      for TreatOpt := DP_SingleDoseNevir to DP_TripleARTDurPreg_Late do
        sumARV := sumARV + DP.GetARVRegimen(DP_PrenatalProphylaxis, TreatOpt, DP_Percent, t);
      numPMTCT := max(0, DP.GetNeedPMTCT(t) * sumARV / 100 + DP.GetPatientsReallocated(t));
      DP.SetOnPMTCT(t, numPMTCT);
      if sumARV > 0 then
        DP.SetPMTCTEffReg(t, numPMTCT * (1.0 - DP.GetARVRegimen(DP_PrenatalProphylaxis, DP_SingleDoseNevir, DP_Percent, t) / sumARV))
      else
        DP.SetPMTCTEffReg(t, 0);

      { Calculate the percentage of HIV+ pregnant women on each regimen }
      need := max(DP.GetNeedPMTCT(t), numPMTCT);
      if (need > 0) and (sumARV > 0) then
        for TreatOpt := DP_SingleDoseNevir to DP_TripleARTDurPreg_Late do
          TreatPercent[TreatOpt] := (DP.GetARVRegimen(DP_PrenatalProphylaxis, TreatOpt, DP_Percent, t) / sumARV) * (numPMTCT / need);
    end;

    sumARV := 0;
    for TreatOpt := DP_SingleDoseNevir to DP_TripleARTDurPreg_Late do
      sumARV := sumARV + TreatPercent[TreatOpt];
    TreatPercent[DP_NoProphylaxis] := max(0.0, 1.0 - sumARV); // take max to avoid negative numbers when sumARV almost 1

   { Add drop-out of those on ART }
   { For 2018 drop-out has been redefined as the percentage remaining
     on ART at the time of delivery }
    TreatPercent[DP_TripleARTBefPreg] := TreatPercent[DP_TripleARTBefPreg]
                                         * DP.GetPercentARTDelivery(DP_OnARTAtDelivery, t) / 100;
    TreatPercent[DP_TripleARTDurPreg] := TreatPercent[DP_TripleARTDurPreg]
                                         * DP.GetPercentARTDelivery(DP_StartingARTAtDelivery, t) / 100;

    { Re-calculate those on no prophylaxis }
    sum3 := 0;
    for TreatOpt := DP_SingleDoseNevir to DP_TripleARTDurPreg_Late do
      sum3 := sum3 + TreatPercent[TreatOpt];
    if sum3 < 1 then
      TreatPercent[DP_NoProphylaxis] := 1 - sum3
      else TreatPercent[DP_NoProphylaxis] := 0;

    sum3 := 0;
    for TreatOpt := DP_DualARV to DP_TripleARTDurPreg_Late do
      sum3 := sum3 + TreatPercent[TreatOpt];

    DP.SetHIVPosBFWomen3Months(t, births *
       ((1 - DP.GetInfantFeedingOptions(2, DP_NotInPMTCT,t) / 100) * (1 - sum3)
      + (1 - DP.GetInfantFeedingOptions(2, DP_InPMTCT,t) / 100) * sum3 ));

    DP.SetHIVPosBFWomen12Months(t, births *
       ((1 - DP.GetInfantFeedingOptions(7, DP_NotInPMTCT,t) / 100) * (1 - sum3)
      + (1 - DP.GetInfantFeedingOptions(7, DP_InPMTCT,t) / 100) * sum3 ));

    { First calculate the proportion of HIV+ pregnant women by CD4 count }
    sum3 := 0;
    sum1 := 0;
    sum2 := 0;
    for a := 15 to 49 do
    begin
      for c := DP_CD4_100_199 to DP_CD4_LT50 do
        sum3 := sum3 + pop1[t, DP_Female, a, c, DP_D_HIVPos];
      sum1 := sum1 + pop1[t, DP_Female, a, DP_CD4_200_249, DP_D_HIVPos]
                   + pop1[t, DP_Female, a, DP_CD4_250_349, DP_D_HIVPos];
      sum2 := sum2 + pop1[t, DP_Female, a, DP_CD4_350_500, DP_D_HIVPos]
                   + pop1[t, DP_Female, a, DP_CD4_GT500, DP_D_HIVPos];
    end;
    if (sum3+sum1+sum2) > 0 then
    begin
      PropLT200 := sum3 / (sum3 + sum1 + sum2);
      Prop200to350 := sum1 / (sum3+ sum1 + sum2);
      PropGE350 := sum2 / (sum3 + sum1 + sum2);
    end
    else
    begin
      PropLT200 := 0;
      Prop200to350 := 1;
      PropGE350 := 0;
    end;
    PropLT350 := PropLT200 + Prop200to350;

    (* Set the number of HIV pregnant women in need of treatment *)
    if DP.GetARVRegimen(DP_PrenatalProphylaxis,DP_TripleARTBefPreg,DP_Number,t) = 0 then
      TripleArtPrePreg := (DP.GetARVRegimen(DP_PrenatalProphylaxis,DP_TripleARTBefPreg,DP_Percent,t)/100) * DP.GetNeedPMTCT(t)
    else
      TripleARTPrePreg := DP.GetARVRegimen(DP_PrenatalProphylaxis,DP_TripleARTBefPreg,DP_Number,t);

    if DP.GetARVRegimen(DP_PrenatalProphylaxis,DP_TripleARTDurPreg,DP_Number,t) = 0 then
      TripleARTDurPreg := (DP.GetARVRegimen(DP_PrenatalProphylaxis,DP_TripleARTDurPreg,DP_Percent,t)/100) * DP.GetNeedPMTCT(t)
    else
      TripleARTDurPreg := DP.GetARVRegimen(DP_PrenatalProphylaxis,DP_TripleARTDurPreg,DP_Number,t);

    DP.SetHIVPregWomenCD4LT350(t,PropLT350 * (DP.GetNeedPMTCT(t) - TripleARTPrePreg - TripleARTDurPreg)
                                                 + TripleARTPrePreg + TripleARTDurPreg);

    { Adjust transmission rates for Opt A and Opt B if
      some have CD4 counts under 350
      Transmission is double if CD4 count is less than 350 }
    if TreatPercent[DP_OptA] + TreatPercent[DP_OptB] > PropGE350 then
    begin
      if PropGE350 <= 0 then ExcessRatio := 0
      else
        ExcessRatio := (TreatPercent[DP_OptA] + TreatPercent[DP_OptB])
                     / PropGE350 - 1;
      OptA_Transmission_Rate := DP.GetTransEffAssump(DP_Data, DP_OptionA, DP_Perinatal)
                                * (1 + ExcessRatio * 1);
      OptB_Transmission_Rate := DP.GetTransEffAssump(DP_Data, DP_OptionB, DP_Perinatal)
                                * (1 + ExcessRatio * 1);
    end
    else
    begin
      OptA_Transmission_Rate := DP.GetTransEffAssump(DP_Data, DP_OptionA, DP_Perinatal);
      OptB_Transmission_Rate := DP.GetTransEffAssump(DP_Data, DP_OptionB, DP_Perinatal);
    end;

    { Calculate transmission from women receiving some type of prophylaxis or treatment }

    { Include drop-out of those on ART
     The median months of pregnancy at first ante-natal visit is 4.2
     from DHS StatCompiler (12/21/2016) across all countries and all surveys
     So both those on ART before the current pregnancy and those starting
     more than 4 weeks before deliver are subject to 9.2 - 4.2 = 5 months
     of drop-out }
    PTR1 := TreatPercent[DP_SingleDoseNevir] * DP.GetTransEffAssump(DP_Data, DP_SingleDoseNev, DP_Perinatal)
            + TreatPercent[DP_DualARV] * DP.GetTransEffAssump(DP_Data, DP_WHO2006DualARV, DP_Perinatal)
            + TreatPercent[DP_OptA] * OptA_Transmission_Rate
            + TreatPercent[DP_OptB] * OptB_Transmission_Rate
            + TreatPercent[DP_TripleARTBefPreg]
              * DP.GetTransEffAssump(DP_Data, DP_ARTStartPrePreg, DP_Perinatal)
            + TreatPercent[DP_TripleARTDurPreg]
              * DP.GetTransEffAssump(DP_Data, DP_ARTStartDurPreg, DP_Perinatal)
            + TreatPercent[DP_TripleARTDurPreg_Late]
              * DP.GetTransEffAssump(DP_Data, DP_ARTStartDurPreg_Late, DP_Perinatal);

    PercentInProgram := TreatPercent[DP_SingleDoseNevir]
      + TreatPercent[DP_DualARV]
      + TreatPercent[DP_OptA]
      + TreatPercent[DP_OptB]
      + TreatPercent[DP_TripleARTBefPreg]
      + TreatPercent[DP_TripleARTDurPreg]
      + TreatPercent[DP_TripleARTDurPreg_Late];

    if PercentInProgram >0 then PTR2[t] := PTR1 / PercentInProgram;

    { Calculate transmission from women not receiving any prophylaxis }
    { First add transmission due to existing infections by CD4 count }
     PTR1 := PTR1 + TreatPercent[DP_NoProphylaxis]
                  * (PropLT200 * DP.GetTransEffAssump(DP_Data, DP_NoProphExistInfCD4LT200, DP_Perinatal)
                    + Prop200to350 * DP.GetTransEffAssump(DP_Data, DP_NoProphExistInfCD4200_350, DP_Perinatal)
                    + PropGE350 * DP.GetTransEffAssump(DP_Data, DP_NoProphExistInfCD4GT350, DP_Perinatal));
    { PTR3 is perintal transmission among HIV+ PW women. It will exclude incident infections. Used in BF calcs}
     PTR3 := PTR1;
    { Next add transmission due to incident infections }
    sum3 := 0; { new infections among women }
    sum2 := 0; { HIV-negative women }

    { sum2 is total population of HIV- women 15-49 weighted for ASFR }
    for a := 15 to 49 do
    begin
      a5 := trunc(a/5) + 1;
      sum2 := sum2 + pop1[t,DP_Female,a,DP_H_No,DP_D_HIVNeg]
                     * DP.GetASFR(t, a5) / DP.GetASFR(t, DP_AllAges);
    end;

    { sum3 is new infections weighted for ASFR }
    for a5 := DP_A15_19 to DP_A45_49 do
      sum3 := sum3 + DP.GetNewInfections(t, a5, DP_Female)
                   * DP.GetASFR(t, a5) / DP.GetASFR(t, DP_AllAges);
    { IncRate is incidence }
    if (sum2 - need) > 0 then IncRate := sum3 / sum2 else IncRate := 0;
    { v3 is new child infections due to incident infections in mothers}
    { Incident infections occur among women giving birth who are not HIV+ }
    v3 := IncRate * 9/12 * (DP.GetBirths(t, DP_AllAges, DP_BothSexes) - need)
             * DP.GetTransEffAssump(DP_Data, DP_NoProphIncidentInf, DP_Perinatal) / 100;
    if (need > 0) then PTR1 := PTR1 + v3 / need * 100;

    DP.SetPerinatalTransmission(t,PTR1);

    { Calculate HIV+ births }

    HIVbirths := births * DP.GetPerinatalTransmission(t) / 100;
    HIVbirthsToExistingInfections := births * PTR3 / 100;

      DP.SetMTCTRate6Wks(t,DP.GetPerinatalTransmission(t));

    { Calculate transmission during breastfeeding }
    { Calculate women with incident infection who are still BF }
    v2 := 0;
    for d := 1 to 18 do
    begin
      v2 := v2 + IncRate / 12 * 2 * (1 - DP.GetInfantFeedingOptions(d, DP_NotInPMTCT,t) / 100);
    end;
    {Calculate new infections due to incident infections }
    v4 := v2 * DP.GetTransEffAssump(DP_Data, DP_NoProphIncidentInf, DP_BreastfeedingGE350)/100;
    IncidentInfectionsBF := (DP.GetBirths(t, DP_AllAges, DP_BothSexes) - DP.GetNeedPMTCT(t)) * v4;

    { Incident infections = (HIV+ mothers - HIV+ births) * v4 / 100 }
    { We use HIV+ mothers rather than all births because v4 has already been adjusted for prevalence }

    { Calculate distribution of women by PMTCT option if numbers on each treatment option are entered }
    sum3 := DP.GetARVRegimen(DP_PrenatalProphylaxis, DP_SingleDoseNevir, DP_Number, t)
           + DP.GetARVRegimen(DP_PrenatalProphylaxis, DP_DualARV, DP_Number, t)
           + DP.GetARVRegimen(DP_PrenatalProphylaxis, DP_OptA, DP_Number, t)
           + DP.GetARVRegimen(DP_PrenatalProphylaxis, DP_OptB, DP_Number, t)
           {+ DP.GetARVRegimen(DP_PostnatalProphylaxis, DP_OptA, DP_Number, t)
           + DP.GetARVRegimen(DP_PostnatalProphylaxis, DP_OptB, DP_Number, t)}
           + DP.GetARVRegimen(DP_PrenatalProphylaxis, DP_TripleARTBefPreg, DP_Number, t)
           + DP.GetARVRegimen(DP_PrenatalProphylaxis, DP_TripleARTDurPreg, DP_Number, t)
           + DP.GetARVRegimen(DP_PrenatalProphylaxis, DP_TripleARTDurPreg_Late, DP_Number, t);
    if (sum3 <> 0) then
    begin
      { Adjust TreatPercent for OptA and OptB for postnatal prophylaxis }
      if(need<>0) then
      begin
        TreatPercent[DP_OptA_BF] := DP.GetARVRegimen(DP_PostnatalProphylaxis, DP_OptA, DP_Number, t) / need;
        TreatPercent[DP_OptB_BF] := DP.GetARVRegimen(DP_PostnatalProphylaxis, DP_OptB, DP_Number, t) / need;
      end;

      { Re-calculate no prophylaxis as the residual }
      sum3 := 0;
      for TreatOpt := DP_SingleDoseNevir to DP_TripleARTDurPreg_Late do
        sum3 := sum3 + TreatPercent[TreatOpt];

      { Adjust the sum for women receiving services outside the region (adds to sum) or
      receiving services in the region but living outside (subtracts from the sum) }
      sum3 := sum3 + DP.GetPatientsReallocated(t);

      if sum3 > 1 then sum3 := 1;
      TreatPercent[DP_NoProphylaxis] := 1 - sum3;
    end
    else
    { Calculate if distribution entries are percentages }
    begin
      for TreatOpt := DP_SingleDoseNevir to DP_TripleARTDurPreg_Late do
        TreatPercent[TreatOpt] := DP.GetARVRegimen(DP_PrenatalProphylaxis, TreatOpt, DP_Percent, t) / 100;
      TreatPercent[DP_OptA_BF] :=
         DP.GetARVRegimen(DP_PostnatalProphylaxis, DP_OptA, DP_Percent, t) / 100
         * (100 - DP.GetARVRegimen(DP_PrenatalProphylaxis, DP_TripleARTBefPreg, DP_Percent, t)
                - DP.GetARVRegimen(DP_PrenatalProphylaxis, DP_TripleARTDurPreg, DP_Percent, t))
                / 100;
     TreatPercent[DP_OptB_BF] :=
         DP.GetARVRegimen(DP_PostnatalProphylaxis, DP_OptB, DP_Percent, t) / 100
         * (100 - DP.GetARVRegimen(DP_PrenatalProphylaxis, DP_TripleARTBefPreg, DP_Percent, t)
                - DP.GetARVRegimen(DP_PrenatalProphylaxis, DP_TripleARTDurPreg, DP_Percent, t))
                / 100;
     TreatPercent[DP_ART_BF] :=
         DP.GetARVRegimen(DP_PostnatalProphylaxis, DP_ART0_12MthsBF, DP_Percent, t) / 100;


      sum3 := TreatPercent[DP_SingleDoseNevir]
              + TreatPercent[DP_DualARV]
              + TreatPercent[DP_OptA]
              + TreatPercent[DP_OptB]
              + TreatPercent[DP_TripleARTBefPreg]
              + TreatPercent[DP_TripleARTDurPreg];

    { Adjust the sum for women receiving services outside the region (adds to sum) or
      receiving services in the region but living outside (subtracts from the sum) }
      if HIVbirths > 0 then sum3 := sum3 + DP.GetPatientsReallocated(t) / HIVbirths;

      TreatPercent[DP_NoProphylaxis] := 1 - sum3;
      if TreatPercent[DP_NoProphylaxis] > 1 then TreatPercent[DP_NoProphylaxis] := 0;

    end;

    for TreatOpt := DP_SingleDoseNevir to DP_OptB_BF do
      DP.SetTreatPercent(t,TreatOpt,TreatPercent[TreatOpt]);

    {Re-calculate TreatPercent removing those who transmitted during perintatl period }
    for TreatOpt := DP_SingleDoseNevir to DP_TripleARTDurPreg_Late do
      TreatPercent[TreatOpt] := TreatPercent[TreatOpt] * (1- DP.GetTransEffAssump(DP_Data, TreatOpt+3, DP_Perinatal)/100);

    { Calculate transmission rate among breastfeeding women }
    CalcBFtransmission(p,t, 1,3);

    { Number of new infections is BFTR times number of HIV+ pregnant women, BFTR is already adjusted for those who have trasnmitted }
    NewInfFromBFLT6 := births * BFTR;
    CumNewInfFromBF := NewInfFromBFLT6;

    CalcBFtransmission(p,t, 4,6);
    NewInfFromBF6TO12 := births * BFTR;
    CumNewInfFromBF := CumNewInfFromBF + NewInfFromBF6TO12;

    if (CumNewInfFromBF + HIVbirths) <= 0 then PropNewInfFromBF := 0
    else
      PropNewInfFromBF := CumNewInfFromBF / (CumNewInfFromBF + HIVbirths);

    { Update PTR1 to include BF transmission in first year }
    if HIVbirths > 0 then
      DP.SetPerinatalTransmission(t,PTR1 * (1 + (CumNewInfFromBF + IncidentInfectionsBF)/ HIVbirths));

    HIVBirths := HIVBirths + DP.GetNosocomialInfections(t);

    DP.SetNewInfantInfections(t, HIVbirths + (CumNewInfFromBF + IncidentInfectionsBF));
  end;
{$ENDREGION}

{$REGION 'DPCTX'}
  procedure TDPProj.CalcNeedForCTX(t : integer);
  var
    s, c, d, a : byte;
  begin
     { Calculate need for CTX }
    { Need for CTX is all births to HIV+ mothers not tested with PCR
      until 18 months plus those tested and found to be HIV+ plus
      all HIV+ children aged 18 months to 14 years }
    NeedForCTX := 0;
    DP.SetChildCTXNeed1To4(t, DP_male, 0);
    DP.SetChildCTXNeed1To4(t, DP_female, 0);

    { births is all births to HIV+ mothers
      ( 1 - DP.GetPCR(t) / 100 is proportion of births not tested with PCR }
    { 1.5 inflates annual births to births in last 18 months }
    for s := DP_Male to DP_Female do
    begin
      TotalNeedForCTX := births * 1.5 * GetDP_SexBirthRatioPercent(DP,t,s);
      { Add HIV+ children aged 18 to 24 months }
      for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
        for d := DP_P_Perinatal to DP_D_Max do
        begin
          TotalNeedForCTX := TotalNeedForCTX
                           + 0.5 * pop1[t,s,1,c,d];
          DP.SetChildCTXNeed1To4(t, s, DP.GetChildCTXNeed1To4(t, s) + 0.5 * pop1[t,s,1,c,d]);
        end;
      { Add HIV+ children aged 2 to 4 years }
      for a := 2 to 4 do
        for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
          for d := DP_P_Perinatal to DP_D_Max do
          begin
            TotalNeedForCTX := TotalNeedForCTX + pop1[t,s,a,c,d];
            DP.SetChildCTXNeed1To4(t, s,
                  DP.GetChildCTXNeed1To4(t, s) + pop1[t,s,a,c,d]);
          end;
      { Add HIV+ children in need of treatment aged 5 to 14 years }
      for a := 5 to 14 do
        for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
          for d := DP_P_Perinatal to DP_P_BF12 do
            TotalNeedForCTX := TotalNeedForCTX + elig[a,s,c,d];
      DP.SetChildARTCalc(DP_NumNeedCot, s, t, TotalNeedForCTX);
    end; { s loop }
    DP.SetChildARTCalc(DP_NumNeedCot, DP_BothSexes, t,
                 DP.GetChildARTCalc(DP_NumNeedCot, DP_Male, t)
               + DP.GetChildARTCalc(DP_NumNeedCot, DP_Female, t));
    TotalNeedForCTX := DP.GetChildARTCalc(DP_NumNeedCot, DP_BothSexes, t);

    { Calculate children on CTX from input which may be % or # }
    if DP.GetChildARTByAgeGroupPerNum(DP_PerChildHIVPosCot, t) = DP_Percent
    then
      TotalOnCTX := TotalNeedForCTX * DP.GetChildTreatInputs(DP_PerChildHIVPosCot,t) / 100
    else
      TotalOnCTX := DP.GetChildTreatInputs(DP_PerChildHIVPosCot,t);

    if TotalNeedForCTX > 0
      then CTXcoverage := TotalOnCTX / TotalNeedForCTX
      else CTXcoverage := 0;
    if CTXcoverage > 1 then
    begin
      CTXcoverage := 1;
      DP_AddLogEntry(GB, DP_DPAIM,'DPCTX','~5-30','CTXcoverage', DP_GT1, CTXcoverage,' ');
    end;

    DP.SetChildARTCalc(DP_NumOnCot, DP_BothSexes, t, TotalOnCTX);
    if TotalNeedForCTX = 0 then
    begin
      DP.SetChildARTCalc(DP_NumOnCot, DP_Male, t, 0);
      DP.SetChildARTCalc(DP_NumOnCot, DP_Female, t, 0);
    end
    else
    begin
      DP.SetChildARTCalc(DP_NumOnCot, DP_Male, t,
                              TotalOnCTX
                              * DP.GetChildARTCalc(DP_NumNeedCot, DP_Male, t)
                              / TotalNeedForCTX);
      DP.SetChildARTCalc(DP_NumOnCot, DP_Female, t,
                              TotalOnCTX
                              * DP.GetChildARTCalc(DP_NumNeedCot, DP_Female, t)
                              / TotalNeedForCTX);
    end;
    { Calculate effect of cotrim }
    { We do not know how long anyone has been on cotrim }
    { We assume that only those for whom cotrim will be effective will be on it }
    CTXeff := 0;
    for d := 1 to 5 do
      CTXeff := CTXeff + DP.GetEffectTreatChild(DP_Data, DP_ChildEffNoArt, d);
    CTXeff := CTXeff / 5;
    { Calculate average ART effect for those on ART 1-4 years }
    CtxArtEff := 0;
    for d := 2 to 5 do
      CTXARTeff := CTXARTeff + DP.GetEffectTreatChild(DP_Data, DP_ChildEffWithART, d);
    CTXArtEff := CTXArtEff / 4;
  end;
{$ENDREGION}

{$REGION 'DPChild'}
  procedure TDPProj.CalcNumChildrenEligTreat(t : integer);
  var
    a, s, d, c, a1: byte;
  begin
      { Age groups for eligibility are 1=<11 months, 2=12-35 months, 3=35-59 months, 4=5+ years }

    { Calculations }
    { If t =2 then it is the first time through, read distribution files }
  //    if (not GB.GetInCalcStateMode(p)) and (t = 2) then
  //    begin
  //      AssignFile(InFile,GB_GetModDataPath(GB_DP) +'Cd4count.txt');
  //      Reset(InFile);
  //      Readln(Infile); { Read past title }
  //      Readln(Infile); { Read past column headings }
  //      for a := DP_A0 to DP_A14 do
  //      begin
  //        read(Infile, Percent_dist[a,1,t]); { Read past age }
  //        for c := 1 to 8 do
  //          read(Infile, Count_dist[a,c,t]);
  //        readln(Infile);
  //      end;
  //      Close(InFile);
  //      AssignFile(InFile,GB_GetModDataPath(GB_DP) +'Cd4percent.txt');
  //      Reset(InFile);
  //      Readln(Infile); { Read past title }
  //      Readln(Infile); { Read past column headings }
  //      for a := DP_A0 to DP_A14 do
  //      begin
  //        read(Infile, Percent_dist[a,1,t]); { Read past age }
  //        for c := 1 to 8 do
  //          read(Infile, Percent_dist[a,c,t]);
  //        readln(Infile);
  //      end;
  //      Close(InFile);
  //    end
  //    else
  //    begin
  //      { Distribution at age 0 is fixed }
  //      for c := 1 to Max_CD_categories do
  //        Count_dist[DP_A0, c, t] := Count_dist[DP_A0, c, t-1];
  //      for a := DP_A1 to DP_A14 do
  //      begin
  //        { Calculate distribution by CD4 count }
  //        { Highest category is just retention in that category from a-1, t-1 }
  //
  //        Count_dist[a, Max_CD_categories, t] := Retention_Rate_CD4_count[Max_CD_categories]
  //                                               * count_dist[a-1,Max_CD_categories,t-1];
  //        sum_start := Count_dist[a-1, Max_CD_categories, t-1];
  //        sum_end := Count_dist[a, Max_CD_categories, t];
  //        { Other categories are retention from a-1, t-1 plus the net number
  //          dropping down from higher categories }
  //        for c := Max_CD_categories-1 downto 1 do
  //        begin
  //          sum_start := sum_start + count_dist[a-1,c,t-1];
  //          Count_dist[a, c, t] := Retention_Rate_CD4_count[c]
  //                                 * (sum_start - sum_end);
  //          sum_end := sum_end + count_dist[a,c,t];
  //        end;
  //        { Calculate distribution by CD4 percent }
  //        { Distriubution at age 0 is fixed }
  //        for c := 1 to Max_CD_categories do
  //          Percent_dist[DP_A0, c, t] := Percent_dist[DP_A0, c, t-1];
  //        { Highest category is just retention in that category from a-1, t-1 }
  //        Percent_dist[a, Max_CD_categories, t] := Retention_Rate_CD4_percent[Max_CD_categories]
  //                                                 * percent_dist[a-1,Max_CD_categories,t-1];
  //        sum_start := Percent_dist[a-1, Max_CD_categories, t-1];
  //        sum_end := Percent_dist[a, Max_CD_categories, t];
  //        { Other categories are retention from a-1, t-1 plus the net number
  //          dropping down from higher categories }
  //        for c := Max_CD_categories-1 downto 1 do
  //        begin
  //          sum_start := sum_start + percent_dist[a-1,c,t-1];
  //
  //          Percent_dist[a, c, t] := Retention_Rate_CD4_percent[c]
  //                                 * (sum_start - sum_end);
  //          sum_end := sum_end + percent_dist[a,c,t];
  //        end;
  //      end; { a loop }
  //    end;

    { Calculate number eligible for treatment }
    for a := 0 to 14 do
      for s := DP_Male to DP_Female do
        for d := DP_D_All to DP_P_BF12 do
          for c:= DP_H_All to DP_CD4_Ped_LT200 do
            elig[a,s,c,d] :=0;

    for a := 0 to 14 do
    begin
      a1 := 0;
      case a of
        0     : a1 := 0;
        1..2  : a1 := 1;
        3..4  : a1 := 2;
        5..14 : a1 := 3;
      end;
      for s := DP_Male to DP_Female do
        for d := DP_P_Perinatal to DP_P_BF12 do
        begin
          if a < 5 then
          begin
            for c:= DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
            begin
              v1 := 0;
              { Include those below the age limit }
              if ((a+1)*12) <= DP.GetAgeHIVChildOnTreatment(t)
                then v1 := pop1[t, s, a, c, d]
                else
                  if (CD_percent_categories[c] <= DP.GetCD4ThreshHold(DP_Percent, a1, t))
                  then v1 := pop1[t, s, a, c, d];
              elig[a, s, c, d] := v1;
            end
          end else
          begin
            for c:= DP_CD4_Ped_GT1000 to DP_CD4_Ped_LT200 do
            begin
              v1 := 0;
              { Include those below the age limit }
              if ((a+1)*12) <= DP.GetAgeHIVChildOnTreatment(t)
                then v1 := pop1[t, s, a, c, d]
              else
                if CD_count_categories[c] <= DP.GetCD4ThreshHold(DP_Number, a1, t)
                then v1 := pop1[t, s, a, c, d];
              elig[a, s, c, d] := v1;
            end;
          end;
        end;
    end;

    { Sum up all asymptomatics }
    EligibleChildren := 0;
    for a := 0 to 14 do
      for s := DP_Male to DP_Female do
        for d := DP_P_Perinatal to DP_P_BF12 do
          if a < 5 then
          begin
            for c:= DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
            begin
              elig[a,s,DP_H_All,DP_D_All] := elig[a,s,DP_H_All,DP_D_All]
                                             + elig[a,s,c,d];
              EligibleChildren := EligibleChildren + elig[a,s,c,d];
            end
          end
          else
          begin
            for c:= DP_CD4_Ped_GT1000 to DP_CD4_Ped_LT200 do
            begin
              elig[a,s,DP_H_All,DP_D_All] := elig[a,s,DP_H_All,DP_D_All]
                                             + elig[a,s,c,d];
              EligibleChildren := EligibleChildren + elig[a,s,c,d];
            end;
          end;

    UnmetNeed := 0;
    UnmetNeedMale := 0;
    UnmetNeedFemale := 0;
    for a := 0 to 14 do
    begin
      UnmetNeedMale := UnmetNeedMale + elig[a,DP_Male,DP_H_All, DP_D_All];
      UnmetNeedFemale := UnmetNeedFemale + elig[a,DP_Female,DP_H_All, DP_D_All];
    end;
    UnmetNeed := UnmetNeedMale + UnmetNeedFemale;

    for s := DP_Male to DP_Female do
    begin
      temp := 0;
      for a := 1 to 4 do
      begin
        temp := temp + elig[a,s,DP_H_All,DP_D_All];
        for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
        begin
          for d := DP_P_Perinatal to DP_P_BF12 do
            temp := temp + elig[a,s,c,d];
          for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
            temp := temp + pop1[t,s,a,c,d];
        end;
      end;
      DP.SetChildARTNeed1To4(t, s, temp);
    end;
//    DP.SetChildARTNeed1To4(t, DP_BothSexes,
//                               DP.GetChildARTNeed1To4(t, DP_Male)
//                               + DP.GetChildARTNeed1To4(t, DP_Female));

    { Start new code for under 2}
    for s := DP_Male to DP_Female do
    begin
      temp := 0;
      for a := 0 to 2 do
      begin
        for c := DP_CD4_Per_GT30 to DP_CD4_Per_LT5 do
        begin
          for d := DP_P_Perinatal to DP_P_BF12 do
            temp := temp + pop1[t,s,a,c,d];
          for d := DP_D_ARTlt6m to DP_D_ARTgt12m do
            temp := temp + pop1[t,s,a,c,d];
        end;
      end;
      DP.SetARTNeedUnder2(t, s, temp);
    end;
    DP.SetARTNeedUnder2(t, DP_BothSexes,DP.GetARTNeedUnder2(t, DP_Male)
                               + DP.GetARTNeedUnder2(t, DP_Female));
  end;

{$ENDREGION}

function TDPProj.DP_Do_Calculations : boolean;
var
  ModRec    : GB_TModuleRec;
  TBOnStack : Boolean;
begin

//  WriteDPToJSON(DP);

  try
    GB.SetModuleCalcSet(p, GB_DP, true);
    TBOnStack := false;
    Try
      Result := ProjectPopulationBySingleYears;

      if not Result then
        GBStrMessageBox(RS(DP_stDataEr),RS(DP_stWarning), mb_OK)
      else
      begin
        DP.SetProjectionValid(True);

        if (GB.ModuleActive(p, GB_HV)) then
          DP_SetHVProjectionValid(p, True);

        if (GB.ModuleActive(p, GB_HA)) then
          DP_SetHAProjectionValid(p, True);
      end;
    except
      on E : Exception do
      begin
        Result := False;
        GBStrMessageBox(RS(DP_stDataEr)+'<LI>'+E.ClassName+' error raised, with message : '+E.Message,
                               RS(DP_stWarning), mb_OK);
      end;
    end;

    (* AIM Calculations *)
    GB.GetModRec (p, GB_AM, ModRec);
    if Result and (ModRec.Status = GB_ACTIVE) then
    begin
      Calc_SubPopulationTableOutputs(DP);

      if (GB.GetFirstYear(GB_FirstProj) <= GB_NativeYear)
        and (not DP.GetIsFitting) then       //added per John's instructions to speed up IRR and CSAVR fitting by ~10%
      begin
        CalcOrphans;
      end;


      (* Calculate impact if LiST is active or the indicator has been selected in AIM and it has not been
         calculated yet. *)
      if ((GB_CS_DeathsAverted in GB.GetModuleCalcSet(p)) or (GB_DP_RequestFromAIMResults in GB.GetModuleCalcSet(p)))
                  and (not (GB_DP_FlagImpact in GB.GetModuleCalcSet(p))) then
      begin
        try
          GB.SetModuleCalcSet(p,GB_DP_FlagImpact, true);

          (* taking DP off the stack so GB_Calculate will let DP calculate for the impact *)
          GB.SetModuleCalcSet(p,GB_DP, false);
          (* Do not calc TB when running impact calculations to save calculation time *)
          TBOnStack := GB_TB in GB.GetModuleCalcSet(p);
          GB.SetModuleCalcSet(p,GB_TB, True);
          DP_AIDSImpactCalculations(DP);
        finally
          GB.SetModuleCalcSet(p,GB_DP, true);
          GB.SetModuleCalcSet(p,GB_DP_FlagImpact, false);
          GB.SetModuleCalcSet(p,GB_CS_DeathsAverted, false);
          GB.SetModuleCalcSet(p,GB_DP_RequestFromAIMResults, false);
          GB.SetModuleCalcSet(p,GB_TB, TBONStack);
        end;
      end;
    end;

    (* Calculate result MVs after calculating other results *)
    if not DP.GetIsFitting then
      calcResultMVs(DP);

    (* if Generate AIM Exception Log is selected in the options then create the ExceptionLog file *)
    if GB.GetExceptLogBoolean then
      DP_WriteExceptionLog(p);
  finally
    GB.SetModuleCalcSet(p,GB_DP, false);
  end;

//  WriteDPToJSON(DP, '-calculated');
end;

procedure TDPProj.SumPopulationByHIVStatus(const t : integer);
var
  s, a, c, d : byte;
begin
  for s := DP_Female downto DP_BothSexes do
    for a := 0 to DP_MaxSingleAges do
    begin
      pop1[t,s,a,DP_H_All,DP_D_All] := 0;
      for c := DP_H_No to DP_H_Max do
      begin
        pop1[t,s,a,c,DP_D_All] := 0;
        for d := DP_D_HIVNeg to DP_D_Max do
          pop1[t,s,a,c,DP_D_All] := pop1[t,s,a,c,DP_D_All] + pop1[t,s,a,c,d];
        pop1[t,s,a,DP_H_All,DP_D_All] := pop1[t,s,a,DP_H_All,DP_D_All] + pop1[t,s,a,c,DP_D_All];
      end;  {c loop}
    end;  {a loop}
end;

procedure TDPProj.Calc_BothSexes;
var
  a, t : byte;
begin

  for t := 1 to DP.GetFinalIndex do
  begin

    DP.SetChildARTNeed1To4(t, DP_BothSexes, DP.GetChildARTNeed1To4(t, DP_Male) +
      DP.GetChildARTNeed1To4(t, DP_Female));

    for a := DP_A0 to DP_A80 do
    begin
      DP.SetOnARTBySingleAge(t, GB_BothSexes, a,
        DP.GetOnARTBySingleAge(t, GB_Male, a) +
        DP.GetOnARTBySingleAge(t, GB_Female, a));

      DP.SetNewInfectionsBySingleAge(t, DP_BothSexes, a,
        DP.GetNewInfectionsBySingleAge(t, DP_Male, a) +
        DP.GetNewInfectionsBySingleAge(t, DP_Female, a));

      DP.SetAIDSDeathsARTSingleAge(t, DP_BothSexes, a,
        DP.GetAIDSDeathsARTSingleAge(t, DP_Male, a) +
        DP.GetAIDSDeathsARTSingleAge(t, DP_Female, a));
    end;
  end;
end;

procedure TDPProj.Calc_ExtraIndicators;
var
  s, a, t : byte;
begin
  for t := 1 to DP.GetFinalIndex do
    for s := GB_BothSexes to GB_Female do
      for a := DP_A0 to DP_A80 do
        DP.SetAIDSDeathsNoARTSingleAge(t, s, a,
          Get_DP_AIDSDeathsBySingleAge(DP, t, s, a - 1, a - 1) -
          DP.GetAIDSDeathsARTSingleAge(t, s, a));

end;

end.




