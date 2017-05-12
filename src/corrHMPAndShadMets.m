if 0
closestFamiliesData = textscan(fopen('/home/fs01/yw595/closestFamilies.txt'),'%s%s','Delimiter','|','HeaderLines',0);
abundsData = textscan(fopen('/home/fs01/yw595/MATLAB/SEED/src/HMPFamiliesAbundsAndOccurs.txt'),'%s%s%s','Delimiter','\t','HeaderLines',0);
allAbundsData = textscan(fopen('/home/fs01/yw595/MATLAB/SEED/src/HMPAllAbunds.txt'),repmat('%s',1,3840),'Delimiter','\t','HeaderLines',0);

modelNames = keys(modelNamesToModels);
a = cellfun(@(x) strsplit(modelNamesToModels(x).modelName,'_'), modelNames, 'UniformOutput', 0);
b = {};
for i=1:length(a)
    temp = a{i};
    b{end+1} = [temp{1} ' ' temp{2}];
end

closestAbunds = [];
closestOccurs = [];

for i=1:length(a)
    speciesName = b{i};
    if sum(strcmp(closestFamiliesData{1},speciesName))~=0
        closestFamily = closestFamiliesData{2}{strcmp(closestFamiliesData{1},speciesName)};
        disp(closestFamily);
        closestAbund = abundsData{2}{strcmp(abundsData{1},closestFamily)};
        closestOccur = abundsData{3}{strcmp(abundsData{1},closestFamily)};
        closestAbunds(i) = str2num(closestAbund);
        closestOccurs(i) = str2num(closestOccur);
        disp(closestAbund);
        disp(closestOccur);
    end
end

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
        [~,intersectIdxs,~] = intersect(b,jthPresentSpecies);
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

if 0
closestFamiliesRevData = textscan(fopen('/home/fs01/yw595/closestFamiliesRev.txt'),'%s%s','Delimiter','|','HeaderLines',0);
HMPFamilies = allAbundsData{1};
relBiomassMat = [];
for i=2:length(allAbundsData)
    i
    for j=1:length(allAbundsData{1})
        matchIdxs = find(strcmp(closestFamiliesRevData{1},allAbundsData{1}{j}));
        if length(matchIdxs) == 1
            closestModel = closestFamiliesRevData{2}{matchIdxs};
            matchIdxs2 = find(strcmp(b,closestModel));
            closestRate = mean(allBiomassRates(matchIdxs2));
            relBiomassMat(j,i-1) = str2num(allAbundsData{i}{j}) / closestRate;
        end
    end
end
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

useForslund = 0;
if useForslund
    modelNames = keys(modelNamesToModels);
    a = cellfun(@(x) strsplit(modelNamesToModels(x).modelName,'_'), modelNames, 'UniformOutput', 0);
    b = {};
    for i=1:length(a)
        temp = a{i};
        b{end+1} = [temp{1} ' ' temp{2}];
    end
    closestTaxaRevData = textscan(fopen('/home/fs01/yw595/closestTaxaForslundRev.txt'),'%s%s','Delimiter','|','HeaderLines',0);
    ForslundData = textscan(fopen('/home/fs01/yw595/MATLAB/SEED/src/ForslundSI3.txt'),'%s%s','Delimiter','|','HeaderLines',0);
    totalBiomass1 = 0; totalBiomass2 = 0;
    totalBiomass1Arr = []; totalBiomass2Arr = [];
    totalFragRank1 = 0; totalFragRank2 = 0;
    totalFragRank1Arr = []; totalFragRank2Arr = [];
    ForslundTaxa1 = ForslundData{1}(strcmp(ForslundData{2},'up'));
    ForslundTaxa2 = ForslundData{1}(strcmp(ForslundData{2},'down'));
    for i=1:length(ForslundTaxa1)
        d = closestTaxaRevData{2}{strcmp(ForslundTaxa1{i},closestTaxaRevData{1})};
        tempIdx = find(strcmp(d,b));
        if sum(tempIdx) > 1
            tempIdx =tempIdx(1);
        end
        totalBiomass1Arr(end+1) = allBiomassRates(tempIdx);
        totalBiomass1 = totalBiomass1 + allBiomassRates(tempIdx);
        for j=1:length(ForslundTaxa1)
            d2 = closestTaxaRevData{2}{strcmp(ForslundTaxa1{j},closestTaxaRevData{1})};
            tempIdx2 = find(strcmp(d2,b));
            if sum(tempIdx2) > 1
                tempIdx2 =tempIdx2(1);
            end
            totalFragRank1 = totalFragRank1 + find(sortIdxsExtra==(tempIdx*length(b)+tempIdx2));
            totalFragRank1Arr(end+1) = find(sortIdxsExtra==(tempIdx*length(b)+tempIdx2));
        end
    end
    for i=1:length(ForslundTaxa2)
        d = closestTaxaRevData{2}{strcmp(ForslundTaxa2{i},closestTaxaRevData{1})};
        tempIdx = find(strcmp(d,b));
        if sum(tempIdx) > 1
            tempIdx = tempIdx(1);
        end
        totalBiomass2Arr(end+1) = allBiomassRates(tempIdx);
        totalBiomass2 = totalBiomass2 + allBiomassRates(tempIdx);
        for j=1:length(ForslundTaxa2)
            d2 = closestTaxaRevData{2}{strcmp(ForslundTaxa2{j},closestTaxaRevData{1})};
            tempIdx2 = find(strcmp(d2,b));
            if sum(tempIdx2) > 1
                tempIdx2 =tempIdx2(1);
            end
            totalFragRank2 = totalFragRank2 + find(sortIdxsExtra==(tempIdx*length(b)+tempIdx2));
            totalFragRank2Arr(end+1) = find(sortIdxsExtra==(tempIdx*length(b)+tempIdx2));
        end        
    end
end

useDavid = 1;
if useDavid
    modelNames = keys(modelNamesToModels);
    a = cellfun(@(x) strsplit(modelNamesToModels(x).modelName,'_'), modelNames, 'UniformOutput', 0);
    b = {};
    for i=1:length(a)
        temp = a{i};
        b{end+1} = [temp{1} ' ' temp{2}];
    end
    closestTaxaData = textscan(fopen('/home/fs01/yw595/closestTaxa.txt'),'%s%s','Delimiter','|','HeaderLines',0);
    closestTaxaRevData = textscan(fopen('/home/fs01/yw595/closestTaxaRev.txt'),'%s%s','Delimiter','|','HeaderLines',0);
    DavidTaxa1Data = textscan(fopen('/home/fs01/yw595/MATLAB/SEED/src/DavidTaxa1.txt'),'%s','Delimiter','\n','HeaderLines',0);
    DavidTaxa2Data = textscan(fopen('/home/fs01/yw595/MATLAB/SEED/src/DavidTaxa2.txt'),'%s','Delimiter','\n','HeaderLines',0);
    totalBiomass1 = 0; totalBiomass2 = 0;
    totalBiomass1Arr = []; totalBiomass2Arr = [];
    totalFragRank1 = 0; totalFragRank2 = 0;
    totalFragRank1Arr = []; totalFragRank2Arr = [];
    uniqTaxa = {};
    useAlt = 1;
    if useAlt
        DavidTaxa1 = DavidTaxa1Data{1};
        DavidTaxa2 = DavidTaxa2Data{1};
        for i=1:length(DavidTaxa1)
            d = closestTaxaRevData{2}{strcmp(DavidTaxa1{i},closestTaxaRevData{1})};
            totalBiomass1Arr(end+1) = allBiomassRates(strcmp(d,b));
            totalBiomass1 = totalBiomass1 + allBiomassRates(strcmp(d,b));
            for j=1:length(DavidTaxa1)
                d2 = closestTaxaRevData{2}{strcmp(DavidTaxa1{j},closestTaxaRevData{1})};
                totalFragRank1 = totalFragRank1 + find(sortIdxsExtra==(find(strcmp(d,b))*length(b)+find(strcmp(d2,b))));
                totalFragRank1Arr(end+1) = find(sortIdxsExtra==(find(strcmp(d,b))*length(b)+find(strcmp(d2,b))));
            end
        end
        for i=1:length(DavidTaxa2)
            d = closestTaxaRevData{2}{strcmp(DavidTaxa2{i},closestTaxaRevData{1})};
            totalBiomass2Arr(end+1) = allBiomassRates(strcmp(d,b));
            totalBiomass2 = totalBiomass2 + allBiomassRates(strcmp(d,b));
            for j=1:length(DavidTaxa2)
                d2 = closestTaxaRevData{2}{strcmp(DavidTaxa2{j},closestTaxaRevData{1})};
                totalFragRank2 = totalFragRank2 + find(sortIdxsExtra==(find(strcmp(d,b))*length(b)+find(strcmp(d2,b))));
                totalFragRank2Arr(end+1) = find(sortIdxsExtra==(find(strcmp(d,b))*length(b)+find(strcmp(d2,b))));
            end
        end
    else
        for i=1:length(a)
            speciesName = b{i};
            if sum(strcmp(closestTaxaData{1},speciesName))~=0
                closestTaxon = closestTaxaData{2}{strcmp(closestTaxaData{1},speciesName)};
                disp(i)
                uniqTaxa{end+1} = speciesName;
                if ~isempty(strcmp(closestTaxon,DavidTaxa1Data{1}))
                    totalBiomass1 = totalBiomass1 + allBiomassRates(i);
                    disp('1')
                    disp(speciesName)
                end
                if ~isempty(strcmp(closestTaxon,DavidTaxa2Data{1}))
                    totalBiomass2 = totalBiomass2 + allBiomassRates(i);
                    disp('2')
                    disp(speciesName)
                end
                disp(closestTaxon);
            end
        end
    end
end

if useDavid || useForslund
    totalBiomassRand1Arr = bootstrap(totalBiomass1Arr,50,3);
    totalBiomassRand2Arr = bootstrap(totalBiomass2Arr,50,3);
    totalFragRankRand1Arr = bootstrap(totalFragRank1Arr,50,floor(length(totalFragRank1Arr)/10));
    totalFragRankRand2Arr = bootstrap(totalFragRank2Arr,50,floor(length(totalFragRank2Arr)/10));
    % for i=1:50
    %     baseperm = totalBiomass2Arr;
    %     ithperm = baseperm(randperm(length(baseperm)));
    %     totalBiomassRand2Arr(end+1) = sum(ithperm(1:length(totalBiomass2Arr)-3));
    %     %totalBiomassRand2Arr(end+1) = sum(ithperm(length(totalBiomass1Arr)+1:end));
    % end

    groupArr = {};
    groupFragArr = {};
    for i=1:length(totalBiomassRand1Arr)
        if useDavid
            groupArr{end+1} = 'Plant';
        end
        if useForslund
            groupArr{end+1} = 'Normal';
        end
    end
    for i=1:length(totalBiomassRand2Arr)
        if useDavid
            groupArr{end+1} = 'Animal';
        end
        if useForslund
            groupArr{end+1} = 'Diabetic';
        end
    end
    for i=1:length(totalFragRankRand1Arr)
        if useDavid
            groupFragArr{end+1} = 'Plant';
        end
        if useForslund
            groupFragArr{end+1} = 'Normal';
        end
    end
    for i=1:length(totalFragRankRand2Arr)
        if useDavid
            groupFragArr{end+1} = 'Animal';
        end
        if useForslund
            groupFragArr{end+1} = 'Diabetic';
        end
    end
    if useDavid
        writeData({[totalBiomassRand1Arr,totalBiomassRand2Arr],groupArr},'/home/fs01/yw595/plantVsAnimal.txt','\t',{'biomass','group'});
        writeData({[totalFragRankRand1Arr,totalFragRankRand2Arr],groupFragArr},'/home/fs01/yw595/plantVsAnimalFrag.txt','\t',{'fragility','group'});
    end
    if useForslund
        writeData({[totalBiomassRand1Arr,totalBiomassRand2Arr],groupArr},'/home/fs01/yw595/diabetesVsNorm.txt','\t',{'biomass','group'});
        writeData({[totalFragRankRand1Arr,totalFragRankRand2Arr],groupFragArr},'/home/fs01/yw595/diabetesVsNormFrag.txt','\t',{'fragility','group'});
    end
end
nonsense = nonsense+1;

biomassAcrossSampsArr = {};
rhoArr = [];
pvalArr = [];
biomassAcrossSampsRelArr = {};
rhoRelArr = [];
pvalRelArr = [];
for i=1:length(allAbundsData{1})
    i
    a = cellfun(@(x) str2num(x{i}), allAbundsData,'UniformOutput',0);
    biomassAcrossSampsArr{end+1} = cellfun(@(x) x(1),a(2:end));
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








