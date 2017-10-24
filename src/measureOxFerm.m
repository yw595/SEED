if 0
%nadhIdxs = cellfun(@(x) ~isempty(regexp(x,'Nicotinamide|NADH')), mergedModel.metNames);
%fadhIdxs = cellfun(@(x) ~isempty(regexp(x,'Nicotinamide|FADH')), mergedModel.metNames);

fermentKEGGs = {'C00186','C00469'};%mergedModel.metKEGGs(nadhIdxs);
oxKEGGs = {'C00007'};

totalOxArr1 = [];
totalOxArr2 = [];
totalFermentArr1 = [];
totalFermentArr2 = [];
for z=1:length(modelNames)
model1 = modelNamesToModels(modelNames{z});
fermentRxnIdxs = [];
oxRxnIdxs = [];
for i=1:length(model1.rxns)
    consumeMetKEGGs = model1.metKEGGs(model1.S(:,i)<0);
    produceMetKEGGs = model1.metKEGGs(model1.S(:,i)>0);
    if length(intersect(produceMetKEGGs,fermentKEGGs))>0 && length(intersect(consumeMetKEGGs,fermentKEGGs))==0
        fermentRxnIdxs(end+1) = i;
    end
    if length(intersect(produceMetKEGGs,oxKEGGs))==0 && length(intersect(consumeMetKEGGs,oxKEGGs))>0
        oxRxnIdxs(end+1) = i;
    end
end
pseudoFlux1 = pseudoFluxDistArr1{z};
pseudoFlux2 = pseudoFluxDistArr2{z};
if length(pseudoFlux1)>0
totalOxArr1(z) = sum(abs(pseudoFlux1(oxRxnIdxs)));
totalOxArr2(z) = sum(abs(pseudoFlux2(oxRxnIdxs)));
totalFermentArr1(z) = sum(abs(pseudoFlux1(fermentRxnIdxs)));
totalFermentArr2(z) = sum(abs(pseudoFlux2(fermentRxnIdxs)));
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
