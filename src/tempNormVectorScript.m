allSubFluxDiffMat = [];
for i=1:10
    [subLabels,subFluxSums,subFluxDiffSums,subFluxStdSums,subFluxDiffScaledSums] = segmentFluxBySubsystem(bigModelTableFlux,normFluxesNormalArr{i}(:,3),1,normFluxesNormalArr{i+1}(:,3),1);
    allSubFluxDiffMat(:,end+1) = subFluxDiffSums;
end
corrWithDiffArr = [];
[~,~,plainDiffArr,~,~] = segmentFluxBySubsystem(bigModelTableFlux,normFluxesNormalArr{1}(:,3),1,normFluxesNormalArr{11}(:,3),1);
for i=1:size(allSubFluxDiffMat,1)
    [rho pval] = corr(allSubFluxDiffMat(i,:)',(1:10)','type','Spearman');
    corrWithDiffArr(end+1) = rho^2;
end
plainDiffArrTemp = plainDiffArr(~isnan(plainDiffArr) & ~isinf(plainDiffArr) & ~isnan(corrWithDiffArr) & ~isinf(corrWithDiffArr));
corrWithDiffArrTemp = corrWithDiffArr(~isnan(plainDiffArr) & ~isinf(plainDiffArr) & ~isnan(corrWithDiffArr) & ~isinf(corrWithDiffArr));
corrWithDiffArr = corrWithDiffArrTemp;
plainDiffArr = plainDiffArrTemp;
corr(corrWithDiffArr',plainDiffArr','type','Spearman')
writeData({plainDiffArr,corrWithDiffArr},[transferDir filesep 'plainVsCorrDiff.txt'],'\t',{'plainDiff','corrDiff'});
[corrWithDiffArr,sortIdxs] = sort(corrWithDiffArr,'descend');
subLabels = addIdxStrings(subLabels(sortIdxs));
writeData({subLabels,corrWithDiffArr},[transferDir filesep 'corrWithDiff.txt'],'\t',{'sub','corr'});
