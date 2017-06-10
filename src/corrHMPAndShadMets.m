%configSEED;

if 0
uniqAllShadMets = unique(allShadMets);
numShadSpecies = zeros(length(uniqAllShadMets),1);
xArr = []; yArr = []; speciesArr = {}; metArr = {}; presabsArr = {};
for i=1:length(modelNames)
    for j=1:length(uniqAllShadMets)
        xArr(end+1) = j;
        yArr(end+1) = i;
        speciesArr{end+1} = modelNamesToModels(modelNames{i}).modelName;
        metArr{end+1} = uniqAllShadMets{j};
        presabsArr{end+1} = 'absent';
        if strcmp(uniqAllShadMets{j},allShadMets{i})
            presabsArr{end} = 'present';
            numShadSpecies(j) = numShadSpecies(j)+1;
        end
    end
end

numSharedShad = [];
for i=1:length(allShadMets)
    numSharedShad(i) = numShadSpecies(strcmp(uniqAllShadMets,allShadMets{i}));
end

writeData({xArr,yArr,speciesArr,metArr,presabsArr},'/home/fs01/yw595/shadMetHeatMap.txt','\t',{'x','y','species','met','presabs'});
end

if 0
sampleFamilies = allAbundsData{1};
nonDestroyedArr = [];
totalBiomassArr = [];
for i=2:length(allAbundsData)
    i
    sampleData = cellfun(@(x) str2num(x), allAbundsData{i});
    totalBiomass = sum(sampleData);
    presentFamilies = sampleFamilies(sampleData~=0);
    uniqModelNames = {};
    for j=1:length(presentFamilies)
        jthPresentSpecies = closestFamiliesData{1}(strcmp(presentFamilies{j},closestFamiliesData{2}));
        [~,intersectIdxs,~] = intersect(modelNamesShort,jthPresentSpecies);
        jthPresentModels = modelNames(intersectIdxs);
        uniqModelNames = union(uniqModelNames,jthPresentModels);
    end
    [~,intersectIdxs,~] = intersect(modelNames,uniqModelNames);
    uniqBiomassRates = allBiomassRates(intersectIdxs);
    totalNondestroyed = sum(uniqBiomassRates);    
    nonDestroyedArr(end+1) = totalNondestroyed;
    totalBiomassArr(end+1) = totalBiomass;
end
end

if 1
HMPFamilies = allAbundsData{1};
relBiomassMat = [];
realMetabCoopArr = zeros(size(allAbundsDataMatrix,2),1);
numInterArr = zeros(size(allAbundsDataMatrix,2),1);

closestRatesArr = [];
closestTwoRatesMat = [];
closestDiffMat = [];
newDiffMatrix = abs(newFBAMatrix(:,:,1)-newFBAMatrix(:,:,2));
for i=1:length(allAbundsDataNames)
    matchIdxs = find(strcmp(closestFamiliesRevData{1}, allAbundsDataNames{i}));
    if length(matchIdxs) == 1
        closestModel = closestFamiliesRevData{2}{matchIdxs};
        matchIdxs2 = find(strcmp(modelNamesShort,closestModel)); matchIdxs2 = matchIdxs2(1);
        closestRatesArr(i) = mean(allBiomassRates(matchIdxs2));
    else
        closestRatesArr(i) = 0;
    end
    for j=1:length(allAbundsDataNames)
        matchIdxsJ = find(strcmp(closestFamiliesRevData{1}, allAbundsDataNames{j}));
        if length(matchIdxsJ) == 1
            closestModelJ = closestFamiliesRevData{2}{matchIdxsJ};
            matchIdxs2J = find(strcmp(modelNamesShort,closestModelJ));
            
            closestTwoRatesMat(i,j) = mean(mean(allTwoBiomasses(matchIdxs2,matchIdxs2J)));
            closestDiffMat(i,j) = mean(mean(newDiffMatrix(matchIdxs2,matchIdxs2J)));
        else
            closestTwoRatesMat(i,j) = 0;
            closestDiffMat(i,j) = 0;
        end
    end
end

if 1
realAbundsCorrMatrix = zeros(size(allAbundsDataMatrix,1),size(allAbundsDataMatrix,1));
for i=1:size(allAbundsDataMatrix,1)
    i
    for j=1:size(allAbundsDataMatrix,1)
        if i~=j
            realAbundsCorrMatrix(i,j) = corr(allAbundsDataMatrix(i,:)',allAbundsDataMatrix(j,:)','Type','Spearman');
        end
    end
end
end
[rhoTwo pvalTwo] = corr(realAbundsCorrMatrix(:),closestTwoRatesMat(:),'Type','Spearman');
writeData({realAbundsCorrMatrix(:),closestTwoRatesMat(:)},'/home/fs01/yw595/abundCorrVsBiomGain.txt','\t',{'abundCorr','biomGain'});
[rhoDiff pvalDiff] = corr(realAbundsCorrMatrix(:),closestDiffMat(:),'Type','Spearman');
writeData({realAbundsCorrMatrix(:),closestDiffMat(:)},'/home/fs01/yw595/abundCorrVsDiff.txt','\t',{'abundCorr','diff'});

for i=1:size(allAbundsDataMatrix,2)
    relBiomassMat(:,i) = allAbundsDataMatrix(:,i)./closestRatesArr';
    numNonzero = sum(allAbundsDataMatrix(:,i)~=0);
    numInterArr(i) = numNonzero*(numNonzero-1)/2;
    for j=1:size(allAbundsDataMatrix,1)
        for k=1:size(allAbundsDataMatrix,1)
            jthRate = 0;
            if allAbundsDataMatrix(j,i)~=0
                jthRate = closestRatesArr(j);
            end
            kthRate = 0;
            if allAbundsDataMatrix(k,i)~=0
                kthRate = closestRatesArr(k);
            end
            twoRate = 0;
            if allAbundsDataMatrix(k,i)~=0 && ...
                    allAbundsDataMatrix(j,i)~=0
                twoRate = closestTwoRatesMat(j,k);
            end
            realMetabCoopArr(i) = realMetabCoopArr(i) + jthRate + ...
                kthRate + twoRate;
        end
    end
end
realMetabCoopArr = realMetabCoopArr./numInterArr;
speciesAbundsArr = sum(allAbundsDataMatrix)';
speciesAbundsArr = speciesAbundsArr(~isinf(realMetabCoopArr) & ~isnan(realMetabCoopArr));
realMetabCoopArr = realMetabCoopArr(~isinf(realMetabCoopArr) & ~isnan(realMetabCoopArr));
[rhoCoop pvalCoop] = corr(realMetabCoopArr, speciesAbundsArr, 'Type','Spearman');
writeData({realMetabCoopArr,speciesAbundsArr},'/home/fs01/yw595/metabCoopVsSpeciesAbunds.txt','\t',{'metabCoop','speciesAbunds'});
end

if 0
relBiomassMatScaled = [];
ithIdx = 0;
for i=1:size(relBiomassMat,1)
    if ~isnan(relBiomassMat(i,1)) && ~isinf(relBiomassMat(i,1))
        disp(i)
        ithIdx = ithIdx+1;
        maxVal = max(relBiomassMat(i,:)');
        %disp(maxVal);
        for j=1:size(relBiomassMat,2)
            %if ~isnan(relBiomassMat(i,j)) && ~isinf(relBiomassMat(i,j))
            if maxVal ~= 0 && relBiomassMat(i,j)~=0               
                relBiomassMatScaled(ithIdx,j) = relBiomassMat(i,j)/maxVal;
            end
        end
    else
        ithIdx = ithIdx+1;
    end
end

corrMatrix = [];
for i1=1:size(relBiomassMatScaled,1)
    disp(i1)
    for i2=1:size(relBiomassMatScaled,1)
        if ~all(relBiomassMatScaled(i1,:)==0) && ~all(relBiomassMatScaled(i2,:)==0) && i1~=i2
            corrMatrix(i1,i2) = corr(relBiomassMatScaled(i1,:)',relBiomassMatScaled(i2,:)','type','Spearman');
        end
    end
end
end

biomassAcrossSampsArr = {};
rhoArr = [];
pvalArr = [];
biomassAcrossSampsRelArr = {};
rhoRelArr = [];
pvalRelArr = [];
for i=1:length(allAbundsData{1})
    i
    aTemp = cellfun(@(x) str2num(x{i}), allAbundsData,'UniformOutput',0);
    biomassAcrossSampsArr{end+1} = cellfun(@(x) x(1),aTemp(2:end));
    biomassAcrossSampsRelArr{end+1} = biomassAcrossSampsArr{end}./totalBiomassArr;
    [rho pval] = corr(biomassAcrossSampsArr{end}',nonDestroyedArr','type','Spearman');
    rhoArr(end+1) = rho; pvalArr(end+1) = pval;
    [rho pval] = corr(biomassAcrossSampsRelArr{end}',nonDestroyedArr','type','Spearman');
    rhoRelArr(end+1) = rho; pvalRelArr(end+1) = pval;
end

writeData({pvalArr,rhoArr},'/home/fs01/yw595/indSpeciesRho.txt','\t',{'pval','rho'});
writeData({pvalRelArr,rhoRelArr},'/home/fs01/yw595/indSpeciesRhoRel.txt','\t',{'pval','rho'});

writeData({nonDestroyedArr,totalBiomassArr},'/home/fs01/yw595/nondestroyedBiomassCorr.txt','\t',{'x','y'});

writeData({numSharedShad,closestAbunds},'/home/fs01/yw595/ShadAbundCorr.txt','\t',{'x','y'});
writeData({numSharedShad,closestOccurs},'/home/fs01/yw595/ShadOccurCorr.txt','\t',{'x','y'});

[uniqShared uniqSharedIdxs] = unique(numSharedShad);

writeData({numSharedShad(uniqSharedIdxs),closestAbunds(uniqSharedIdxs)},'/home/fs01/yw595/ShadAbundCorrUnique.txt','\t',{'x','y'});
writeData({numSharedShad(uniqSharedIdxs),closestOccurs(uniqSharedIdxs)},'/home/fs01/yw595/ShadOccurCorrUnique.txt','\t',{'x','y'});








