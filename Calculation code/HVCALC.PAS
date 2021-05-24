unit HVcalc;

{ This unit does the calculations for the HIV Vaccine module }

interface

uses
  GBDefs, HVdefs, HVData;

{ Calculates results for a standard HV indicator }
procedure HV_CalcIndicator(IndicatorNum, FirstProj, LastProj, Sex, Interval : byte;
                            SelSet : GB_TByteSet;
                            DisplayValues: HV_PCalcResult;
                            FirstIndex, FinalIndex: integer);

procedure HV_CalcAllRisk(HV:THVData; IndId,Sex : byte;
                          SelSet : GB_TByteSet;
                          FirstIndex, FinalIndex: integer);

implementation

uses
  GBData, HVConst, HVLink;

{ Calculates results for a standard HV indicator }
procedure HV_CalcIndicator(IndicatorNum, FirstProj, LastProj, Sex,
  Interval : byte; SelSet : GB_TByteSet; DisplayValues : HV_PCalcResult;
  FirstIndex, FinalIndex : integer);
var
  p,sid,t,v,h,r,m : byte;
  numerator,
  denominator     : double;
  HV              : THVData;
begin
  (* Note: HV_CalcIndicator assumes that the all the parmeter Sets are *)
  (* valid. It is up to HVDsCnFm to set these sets correctly. *)

  { Fill DisplayValues that will be affected with 0's }
  for p := FirstProj to LastProj do
    for t := FirstIndex to FinalIndex do
      DisplayValues^[p, t] := 0.0;

   { Fill the display values array with data from the appropriate
    calculation or from the BCData object }
   for p := FirstProj to LastProj do
   begin
    HV := GetHVData(p);
    FinalIndex := GetGBFinalYearIdx(p);

    for t := FirstIndex to FinalIndex do
      Case IndicatorNum of
        HV_cmNewInfections :
        begin
          DisplayValues^[p,t] := HV.GetTotalNewInfection(t);
        end;
        HV_cmNumHIVPos:
        begin
          DisplayValues^[p,t] := HV.GetTotalAdultsHIV(HV_Male,t) + HV.GetTotalAdultsHIV(HV_Female,t);
        end;
        HV_cmAIDSDeaths:
        begin
          DisplayValues^[p,t] := HV.GetTotalAIDSDeaths(t);
        end;
        HV_cmAIDSDeathsART:
        begin
          for r := HV_AllRisk to HV_MSMIDU do
            DisplayValues^[p,t] := DisplayValues^[p,t] + HV.GetAIDSDeathsART(sex,r,t);
        end;
        HV_cmNumVaccs:
        begin
          for sid := 1 to HV_max_services do
            if sex = HV_BothSexes then
              DisplayValues^[p,t] := DisplayValues^[p,t] + HV.GetNewVaccinations(HV_Male,sid,t) + HV.GetNewVaccinations(HV_Female,sid,t)
            else
              DisplayValues^[p,t] := DisplayValues^[p,t] + HV.GetNewVaccinations(sex,sid,t);
        end;
        HV_cmNumVaccinated:
        begin
          if sex = HV_BothSexes then
            DisplayValues^[p,t] := HV.GetTotalVaccinated(HV_Male,t)+HV.GetTotalVaccinated(HV_Female,t)
          else
            DisplayValues^[p,t] := HV.GetTotalVaccinated(sex,t);
        end;
        HV_cmTotalPrEP:
        begin
          for r := HV_LRH to HV_MSMIDU do
            if ((HV_GetRNNumMSMRiskGroups(p) = 0) and (not (r in [HV_MSMLR..HV_MSMIDU]))) or
                                          ((HV_GetRNNumMSMRiskGroups(p) = 1) and (r <> HV_MSM)) then
              for m := HV_RN_PrEPOral to HV_RN_PrEPRing do
                if sex = HV_BothSexes then
                  DisplayValues^[p,t] := DisplayValues^[p,t] +
                                             (HV.GetAdults(HV_Male,r,HV_Negative,HV_RN_AllVacc,t)
                                              *HV_GetRNPrEPCoverage(p,HV_Male,r,t)/100
                                              *HV_GetRNMethodMix(p,HV_Male,r,m,t)/100)

                                          + (HV.GetAdults(HV_Female,r,HV_Negative,HV_RN_AllVacc,t)
                                             *HV_GetRNPrEPCoverage(p,HV_Female,r,t)/100
                                             *HV_GetRNMethodMix(p,HV_Female,r,m,t)/100)
                else
                  DisplayValues^[p,t] := DisplayValues^[p,t] +
                                            (HV.GetAdults(sex,r,HV_Negative,HV_RN_AllVacc,t)
                                             *HV_GetRNPrEPCoverage(p,sex,r,t)/100
                                             *HV_GetRNMethodMix(p,sex,r,m,t)/100);
        end;
        HV_cmPrEPSingleMethod:
        begin
          for r := HV_LRH to HV_MSMIDU do
            if ((HV_GetRNNumMSMRiskGroups(p) = 0) and (not (r in [HV_MSMLR..HV_MSMIDU]))) or
                                          ((HV_GetRNNumMSMRiskGroups(p) = 1) and (r <> HV_MSM)) then
              for m := HV_RN_PrEPOral to HV_RN_PrEPRing do
              begin
                if m in SelSet then
                  if sex = HV_BothSexes then
                    DisplayValues^[p,t] := DisplayValues^[p,t] +
                                               (HV.GetAdults(HV_Male,r,HV_Negative,HV_RN_AllVacc,t)
                                                *HV_GetRNPrEPCoverage(p,HV_Male,r,t)/100
                                                *HV_GetRNMethodMix(p,HV_Male,r,m,t)/100)

                                            + (HV.GetAdults(HV_Female,r,HV_Negative,HV_RN_AllVacc,t)
                                               *HV_GetRNPrEPCoverage(p,HV_Female,r,t)/100
                                               *HV_GetRNMethodMix(p,HV_Female,r,m,t)/100)
                  else
                    DisplayValues^[p,t] := DisplayValues^[p,t] +
                                              (HV.GetAdults(sex,r,HV_Negative,HV_RN_AllVacc,t)
                                               *HV_GetRNPrEPCoverage(p,sex,r,t)/100
                                               *HV_GetRNMethodMix(p,sex,r,m,t)/100);
              end;
        end;
        HV_cmIncidence:
        begin
           DisplayValues^[p,t] := HV.GetIncidence(t)*100;
        end;
        HV_cmPrevalence:
        begin
          if Sex = HV_BothSexes then
            DisplayValues^[p,t] := ((HV.GetTotalAdultsHIV(HV_Male,t)+HV.GetTotalAdultsHIV(HV_Female,t))/HV.GetPopulations(HV_BothSexes,t))*100
          else
            DisplayValues^[p,t] := ( (HV.GetTotalAdultsHIV(Sex,t))/HV.GetPopulations(Sex,t))*100;
        end;
        HV_cmPrevRG:
        begin
          for sid := 1 to HV_max_services do
            if sid in SelSet then
              if sex = HV_BothSexes then
                DisplayValues^[p,t] := DisplayValues^[p,t] + (HV.GetCalcPrevalence(HV_Male,sid,t)*HV.GetPercentPop(HV_Male,t) + HV.GetCalcPrevalence(HV_Female,sid,t)*HV.GetPercentPop(HV_Female,t))*100
              else
                DisplayValues^[p,t] := DisplayValues^[p,t] + HV.GetCalcPrevalence(sex,sid,t)*100;
        end;
        HV_cmNewInfectionsRG:
        begin
          for v := HV_RN_Unv to HV_RN_NoProt do
            for sid := 1 to HV_max_services do
              if Sid in SelSet then
                if sex = HV_BothSexes then
                  DisplayValues^[p,t] := DisplayValues^[p,t] + HV.GetNewInfections(HV_Male,sid,v,t) + HV.GetNewInfections(HV_Female,sid,v,t)
                else
                DisplayValues^[p,t] := DisplayValues^[p,t] + HV.GetNewInfections(sex,sid,v,t);
        end;
        HV_cmAIDSDeathsRG:
        begin
          for v := HV_RN_Unv to HV_RN_NoProt do
            for sid := 1 to HV_max_services do
              if Sid in SelSet then
                if Sid = HV_None then
                  DisplayValues^[p,t] := 0
                else
                  if sex = HV_BothSexes then
                    DisplayValues^[p,t] := DisplayValues^[p,t] + HV.GetAIDSDeaths(HV_Male,sid,v,t) + HV.GetAIDSDeaths(HV_Female,sid,v,t)
                  else
                    DisplayValues^[p,t] := DisplayValues^[p,t] + HV.GetAIDSDeaths(sex,sid,v,t);
        end;
        HV_cmPopbyRG:
        begin
          for h := HV_Negative to HV_CD4_LT50_ART do
            if h <> HV_AllHIV then
              for v := HV_RN_Unv to HV_RN_NoProt do
                for sid := 1 to HV_max_services do
                  if Sid in SelSet then
                    if sex = HV_BothSexes then
                      DisplayValues^[p,t] := DisplayValues^[p,t] + HV.GetAdults(HV_Male,sid,h,v,t) + HV.GetAdults(HV_Female,sid,h,v,t)
                    else
                     DisplayValues^[p,t] := DisplayValues^[p,t] + HV.GetAdults(sex,sid,h,v,t);
        end;
        HV_cmARTCoveragebyRG:
        begin
          for sid := 1 to HV_max_services do
            if Sid in SelSet then
              DisplayValues^[p,t] := DisplayValues^[p,t] + HV.GetARTCoverageByRG(sex,sid,t);
        end;
        HV_cmTotalARTCoverageAR:
        begin
          DisplayValues^[p,t] := HV.GetTotalARTCoverage(sex,t);
        end;
        HV_cmNumART:
        begin
          for v := HV_RN_Unv to HV_RN_NoProt do
            for sid := 1 to HV_max_services do
              for h := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
                DisplayValues^[p,t] := DisplayValues^[p,t] + HV.GetAdults(HV_Male,sid,h,v,t) +
                                               HV.GetAdults(HV_Female,sid,h,v,t);
        end;
        HV_cmHIVPosRG:
        begin
          for sid := 1 to HV_max_services do
          begin
            for h := HV_Primary to HV_CD4_LT50 do
              if Sid in SelSet then
                 if sex = HV_BothSexes then
                   DisplayValues^[p,t] := DisplayValues^[p,t] + HV.GetAdults(HV_Male,sid,h,HV_RN_AllVacc,t) + HV.GetAdults(HV_Female,sid,h,HV_RN_AllVacc,t)
                 else
                   DisplayValues^[p,t] := DisplayValues^[p,t] + HV.GetAdults(sex,sid,h,HV_RN_AllVacc,t);

            for h := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
              if Sid in SelSet then
                 if sex = HV_BothSexes then
                   DisplayValues^[p,t] := DisplayValues^[p,t] + HV.GetAdults(HV_Male,sid,h,HV_RN_AllVacc,t) + HV.GetAdults(HV_Female,sid,h,HV_RN_AllVacc,t)
                 else
                   DisplayValues^[p,t] := DisplayValues^[p,t] + HV.GetAdults(sex,sid,h,HV_RN_AllVacc,t);
          end;{sid loop}
        end;
        HV_cmNumVaccsRG:
        begin
          for sid := 1 to HV_max_services do
            if sid in SelSet then
              if sex = HV_BothSexes then
                DisplayValues^[p,t] := DisplayValues^[p,t] + HV.GetNewVaccinations(HV_Male,sid,t) + HV.GetNewVaccinations(HV_Female,sid,t)
              else
                DisplayValues^[p,t] := DisplayValues^[p,t] + HV.GetNewVaccinations(sex,sid,t);
        end;
        HV_cmPrEPRG:
        begin
          for sid := 1 to HV_max_services do
            if sid in SelSet then
              if sid = 1 then
                DisplayValues^[p,t] := 0
              else
              begin
                for v := HV_RN_PrEPOral to HV_RN_PrEPRing do
                    if sex = HV_BothSexes then
                      DisplayValues^[p,t] := DisplayValues^[p,t] +
                                                 (HV.GetAdults(HV_Male,sid,HV_Negative,HV_RN_AllVacc,t)
                                                  *HV_GetRNPrEPCoverage(p,HV_Male,sid,t)/100
                                                  *HV_GetRNMethodMix(p,HV_Male,sid,v,t)/100)

                                              + (HV.GetAdults(HV_Female,sid,HV_Negative,HV_RN_AllVacc,t)
                                                 *HV_GetRNPrEPCoverage(p,HV_Female,sid,t)/100
                                                 *HV_GetRNMethodMix(p,HV_Female,sid,v,t)/100)
                    else
                      DisplayValues^[p,t] := DisplayValues^[p,t] +
                                                (HV.GetAdults(sex,sid,HV_Negative,HV_RN_AllVacc,t)
                                                 *HV_GetRNPrEPCoverage(p,sex,sid,t)/100
                                                 *HV_GetRNMethodMix(p,sex,sid,v,t)/100);
              end;
        end;
        HV_cmIA:
        begin
          DisplayValues^[p,t] := HV.GetInfectionsAverted(t);
        end;
        HV_cmCumIA:
        begin
          DisplayValues^[p,t] := HV.GetCumInfectionsAverted(t);
        end;
        HV_cmDA:
        begin
          DisplayValues^[p,t] := HV.GetDeathsAverted(t);
        end;
        HV_cmCumDA:
        begin
          DisplayValues^[p,t] := HV.GetCumDeathsAverted(t);
        end;
        HV_cmCumCostIA:
        begin
         DisplayValues^[p,t] := HV_GetRNCumCostIA(p,t);
        end;
        HV_cmCumCostDA:
        begin
          DisplayValues^[p,t] := HV_GetRNCumCostDA(p,t);
        end;
        HV_cmTotalCosts2:
        begin
          DisplayValues^[p,t] := HV_GetRNTotalCosts(p,t);
        end;
        HV_cmNewlyOnARTTotal:
        begin
          for h := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
            if sex = HV_BothSexes then
              DisplayValues^[p,t] := DisplayValues^[p,t] + HV.GetNewlyOnART(HV_Male,h,t)
                                              + HV.GetNewlyOnART(HV_Female,h,t)
            else
              DisplayValues^[p,t] := DisplayValues^[p,t] + HV.GetNewlyOnART(sex,h,t);
        end;
        HV_cmNewlyEligibleARTTotal:
        begin
          for h := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
            if sex = HV_BothSexes then
              DisplayValues^[p,t] := DisplayValues^[p,t] + HV.GetNewlyEligibleART(HV_Male,h,t)
                                              + HV.GetNewlyEligibleART(HV_Female,h,t)
            else
              DisplayValues^[p,t] := DisplayValues^[p,t] + HV.GetNewlyEligibleART(sex,h,t);
        end;
        HV_cmMultiplier:
        begin
          DisplayValues^[p,t] := HV.GetrMultAll(HV_BothSexes,HV_AllRisk,t);
        end;
        HV_cmMultiplierbyRG:
        begin
          for sid := 1 to HV_max_services do
            if sid in SelSet then
              DisplayValues^[p,t] := DisplayValues^[p,t] + HV.GetrMultAll(sex,sid,t);
        end;
        HV_cmIncidenceRG:
        begin
          for sid := 1 to HV_max_services do
          begin
            if Sid in SelSet then
            begin
              numerator := 0;
              denominator := 0;
              for v := HV_RN_Unv to HV_RN_NoProt do
              begin
                if sex = HV_BothSexes then
                begin
                  numerator := numerator + HV.GetNewInfections(HV_Male,sid,v,t) +
                                           HV.GetNewInfections(HV_Female,sid,v,t);
                  denominator := denominator + HV.GetAdults(HV_Male,sid,HV_Negative,v,t) +
                                               HV.GetAdults(HV_Female,sid,HV_Negative,v,t);
                end
                else
                begin
                  numerator := numerator + HV.GetNewInfections(sex,sid,v,t);
                  denominator := denominator + HV.GetAdults(sex,sid,HV_Negative,v,t);
                end;
              end;{for v}
              if denominator > 0 then
                DisplayValues^[p,t] := (numerator/denominator) * 100
              else
                DisplayValues^[p,t] := 0;
            end;
          end;{for p,t,sid}
        end;
      end;
   end;
end;


procedure HV_CalcAllRisk(HV:THVData; IndId,Sex : byte;SelSet : GB_TByteSet;FirstIndex, FinalIndex: integer);
var
  sid,t,v,h,m,r,
  SubSet,
  StartSid,
  EndSid          : byte;
  value,
  numerator,
  denominator     : double;
begin
  HV.SetMaxAllRisk(0);
  for t := FirstIndex to FinalIndex do
  begin
    SubSet := 0;
    if IndId = HV_cmPrEPAllMethods then
    begin
      StartSid := HV_RN_PrEPOral;
      EndSid := HV_MaxPrEPMethods;
    end
    else
    begin
      StartSid := 1;
      EndSid := HV_max_services;
    end;
    for sid := StartSid to EndSid do
    begin
      if sid in SelSet then
      begin
        Case IndId of
          HV_cmPopAR:
              begin
               value := 0;
               for h := HV_Negative to HV_CD4_LT50_ART do
                  if h <> HV_AllHIV then
                    begin
                      for v := HV_RN_Unv to HV_RN_NoProt do
                      begin
                        if sex = HV_BothSexes then
                          value := value + (HV.GetAdults(HV_Male,sid,h,v,t) + HV.GetAdults(HV_Female,sid,h,v,t))
                        else
                          value := value + (HV.GetAdults(sex,sid,h,v,t));
                      end;
                    end;
               HV.SetAllRiskOutput(SubSet,t,value);
               if value > HV.GetMaxAllRisk then
                HV.SetMaxAllRisk(value);
               inc(SubSet);
              end;
          HV_cmNewInfectionsAR:
              begin
               value := 0;
               for v := HV_RN_Unv to HV_RN_NoProt do
               begin
                 if sex = HV_BothSexes then
                    value := value + HV.GetNewInfections(HV_Male,sid,v,t) + HV.GetNewInfections(HV_Female,sid,v,t)
                 else
                    value := value + HV.GetNewInfections(sex,sid,v,t);
               end;
               HV.SetAllRiskOutput(SubSet,t,value);
               if value > HV.GetMaxAllRisk then
                HV.SetMaxAllRisk(value);
               inc(SubSet);
              end;
          HV_cmPrevAR:
              begin
               if sex = HV_BothSexes then
                  value := (HV.GetCalcPrevalence(HV_Male,sid,t)*HV.GetPercentPop(HV_Male,t) + HV.GetCalcPrevalence(HV_Female,sid,t)*HV.GetPercentPop(HV_Female,t))*100
              else
                   value := HV.GetCalcPrevalence(sex,sid,t)*100;

               HV.SetAllRiskOutput(SubSet,t,value);
               if value > HV.GetMaxAllRisk then
                HV.SetMaxAllRisk(value);
               inc(SubSet);
              end;
          HV_cmHIVPosAR:
              begin
               value := 0;
               for h := HV_Primary to HV_CD4_LT50 do
                 for v := HV_RN_UnV to HV_RN_NoProt do
                 begin
                   if sex = HV_BothSexes then
                      value := value + HV.GetAdults(HV_Male,sid,h,v,t) + HV.GetAdults(HV_Female,sid,h,v,t)
                   else
                      value := value + HV.GetAdults(sex,sid,h,v,t);
                 end;

               for h := HV_CD4_GT500_ART to HV_CD4_LT50_ART do
                 for v := HV_RN_UnV to HV_RN_NoProt do
                 begin
                   if sex = HV_BothSexes then
                      value := value + HV.GetAdults(HV_Male,sid,h,v,t) + HV.GetAdults(HV_Female,sid,h,v,t)
                   else
                      value := value + HV.GetAdults(sex,sid,h,v,t);
                 end;
               HV.SetAllRiskOutput(SubSet,t,value);
               if value > HV.GetMaxAllRisk then
                HV.SetMaxAllRisk(value);
               inc(SubSet);
              end;
          HV_cmAIDSDeathsAR:
              begin
               value := 0;
               for v := HV_RN_Unv to HV_RN_NoProt do
               begin
                 if Sid = HV_None then
                    Value := 0
                 else
                 begin
                   if sex = HV_BothSexes then
                      value := value + HV.GetAIDSDeaths(HV_Male,sid,v,t) + HV.GetAIDSDeaths(HV_Female,sid,v,t)
                   else
                      value := value + HV.GetAIDSDeaths(sex,sid,v,t);
                  end
               end;
               HV.SetAllRiskOutput(SubSet,t,value);
               if value > HV.GetMaxAllRisk then
                HV.SetMaxAllRisk(value);
               inc(SubSet);
              end;
          HV_cmNumVaccsAR:
              begin
                if sex = HV_BothSexes then
                   value := HV.GetNewVaccinations(HV_Male,sid,t) + HV.GetNewVaccinations(HV_Female,sid,t)
                else
                   value := HV.GetNewVaccinations(sex,sid,t);
                HV.SetAllRiskOutput(SubSet,t,value);
                if value > HV.GetMaxAllRisk then
                  HV.SetMaxAllRisk(value);
                inc(SubSet);
              end;
          HV_cmPrEPAR:
              begin
                if sid = 1 then
                  value := 0
                else
                begin
                  value := 0;
                  for m := HV_RN_PrEPOral to HV_RN_PrEPRing do
                    if sex = HV_BothSexes then
                      value := value + (HV.GetAdults(HV_Male,sid,HV_Negative,HV_RN_AllVacc,t)
                                                  *HV_GetRNPrEPCoverage(HV.Proj,HV_Male,sid,t)/100
                                                  *HV_GetRNMethodMix(HV.Proj,HV_Male,sid,m,t)/100)

                                     + (HV.GetAdults(HV_Female,sid,HV_Negative,HV_RN_AllVacc,t)
                                                 *HV_GetRNPrEPCoverage(HV.Proj,HV_Female,sid,t)/100
                                                 *HV_GetRNMethodMix(HV.Proj,HV_Female,sid,m,t)/100)
                    else
                      value := value + (HV.GetAdults(sex,sid,HV_Negative,HV_RN_AllVacc,t)
                                                 *HV_GetRNPrEPCoverage(HV.Proj,sex,sid,t)/100
                                                 *HV_GetRNMethodMix(HV.Proj,sex,sid,m,t)/100);
                end;
                HV.SetAllRiskOutput(SubSet,t,value);
                if value > HV.GetMaxAllRisk then
                  HV.SetMaxAllRisk(value);
                inc(SubSet);
              end;
          HV_cmPrEPAllMethods:
              begin
                value := 0;
                for r := HV_LRH to HV_MSMIDU do
                                      if ((HV_GetRNNumMSMRiskGroups(HV.Proj) = 0) and (not (r in [HV_MSMLR..HV_MSMIDU]))) or
                                          ((HV_GetRNNumMSMRiskGroups(HV.Proj) = 1) and (r <> HV_MSM)) then
                  if sex = HV_BothSexes then
                    value := value + (HV.GetAdults(HV_Male,r,HV_Negative,HV_RN_AllVacc,t)
                                                *HV_GetRNPrEPCoverage(HV.Proj,HV_Male,r,t)/100
                                                *HV_GetRNMethodMix(HV.Proj,HV_Male,r,sid,t)/100)

                                   + (HV.GetAdults(HV_Female,r,HV_Negative,HV_RN_AllVacc,t)
                                               *HV_GetRNPrEPCoverage(HV.Proj,HV_Female,r,t)/100
                                               *HV_GetRNMethodMix(HV.Proj,HV_Female,r,sid,t)/100)
                  else
                    value := value + (HV.GetAdults(sex,r,HV_Negative,HV_RN_AllVacc,t)
                                               *HV_GetRNPrEPCoverage(HV.Proj,sex,r,t)/100
                                               *HV_GetRNMethodMix(HV.Proj,sex,r,sid,t)/100);
                HV.SetAllRiskOutput(SubSet,t,value);
                if value > HV.GetMaxAllRisk then
                  HV.SetMaxAllRisk(value);
                inc(SubSet);
              end;
          HV_cmMultiplierAR:
              begin
                value := HV.GetrMultAll(HV_Male,sid,t) + HV.GetrMultAll(HV_Female,sid,t);
                HV.SetAllRiskOutput(SubSet,t,value);
                if value > HV.GetMaxAllRisk then
                  HV.SetMaxAllRisk(value);
                inc(SubSet);
              end;
          HV_cmIncidenceAR:
              begin
                numerator := 0;
                denominator := 0;
                for v := HV_RN_Unv to HV_RN_NoProt do
                begin
                  if sex = HV_BothSexes then
                  begin
                    numerator := numerator + HV.GetNewInfections(HV_Male,sid,v,t) +
                                             HV.GetNewInfections(HV_Female,sid,v,t);
                    denominator := denominator + HV.GetAdults(HV_Male,sid,HV_Negative,v,t) +
                                                 HV.GetAdults(HV_Female,sid,HV_Negative,v,t);
                  end
                  else
                  begin
                    numerator := numerator + HV.GetNewInfections(sex,sid,v,t);
                    denominator := denominator + HV.GetAdults(sex,sid,HV_Negative,v,t);
                  end;
                end;{for v}
                if denominator > 0 then
                  value := (numerator/denominator) * 100
                else
                  value := 0;
                HV.SetAllRiskOutput(SubSet,t,value);
                if value > HV.GetMaxAllRisk then
                  HV.SetMaxAllRisk(value);
                inc(SubSet);
              end;
        end; {case end}
      end; {sid in Select}
    end; {sid loop}
  end;{t loop}
end;

end.
