totalOxArr1 = [];
totalOxArr2 = [];
totalFermentArr1 = [];
totalFermentArr2 = [];

for z=1:length(modelNames)
    model1 = modelNamesToModels(modelNames{z});
    if 0
        pseudoFlux1 = pseudoFluxDistArr1{z};
	pseudoFlux2 = pseudoFluxDistArr2{z};
	if length(pseudoFlux1)>0
	    [totalOx1 totalFerm1] = measureOxFermFunc(model1,pseudoFlux1);
            [totalOx2 totalFerm2] = measureOxFermFunc(model1,pseudoFlux2);
	    totalOxArr1(z) = totalOx1;
	    totalOxArr2(z) = totalOx2;
	    totalFermentArr1(z) = totalFerment1;
	    totalFermentArr2(z) = totalFerment2;
	end
    end
end
	
totalOxArr1Temp = totalOxArr1;
nonzeroIdxs = find(totalOxArr1Temp~=0);
for i=1:length(totalOxArr1Temp)
    if totalOxArr1Temp(i)==0
	replaceVal = totalOxArr1Temp(nonzeroIdxs(randi(length(nonzeroIdxs),1)));
        totalOxArr1Temp(i) = replaceVal;
    end
end
totalFermentArr1Temp = totalFermentArr1;
nonzeroIdxs = find(totalFermentArr1Temp~=0);
for i=1:length(totalFermentArr1Temp)
    if totalFermentArr1Temp(i)==0
	replaceVal = totalFermentArr1Temp(nonzeroIdxs(randi(length(nonzeroIdxs),1)));
        totalFermentArr1Temp(i) = replaceVal;
    end
end
totalOxArr2Temp = totalOxArr2;
nonzeroIdxs = find(totalOxArr2Temp~=0);
for i=1:length(totalOxArr2Temp)
    if totalOxArr2Temp(i)==0
	replaceVal = totalOxArr2Temp(nonzeroIdxs(randi(length(nonzeroIdxs),1)));
        totalOxArr2Temp(i) = replaceVal;
    end
end
totalFermentArr2Temp = totalFermentArr2;
nonzeroIdxs = find(totalFermentArr2Temp~=0);
for i=1:length(totalFermentArr2Temp)
    if totalFermentArr2Temp(i)==0
	replaceVal = totalFermentArr2Temp(nonzeroIdxs(randi(length(nonzeroIdxs),1)));
        totalFermentArr2Temp(i) = replaceVal;
    end
end
totalRatioArr1 = totalOxArr1Temp./totalFermentArr1Temp;
totalRatioArr2 = totalOxArr2Temp./totalFermentArr2Temp;
totalRatioArr1 = totalRatioArr1(~isinf(totalRatioArr1) & ~isnan(totalRatioArr1));
totalRatioArr2 = totalRatioArr2(~isinf(totalRatioArr2) & ~isnan(totalRatioArr2));
