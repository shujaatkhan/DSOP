function scriptkappa = ScriptKappa(mt,mLowerBoundLife,DeltaGothicHLife,KappaMin,PeriodsUntilT,chiIntData)

if PeriodsUntilT == 0;
    scriptkappa = ones(size(mt));
else
    scriptkappa = KappafromChi(mt,mLowerBoundLife,DeltaGothicHLife,KappaMin,PeriodsUntilT,chiIntData);
%     for i =1:length(mt)
%         scriptc = ScriptC(mt(i),mLowerBoundLife,DeltaGothicHLife,KappaMin,PeriodsUntilT,chiIntData,bl);
%         if scriptc == mt(i)+bl||scriptc==10^-6
%             scriptkappa(i) = 1;
%         end
%     end
            
end

end