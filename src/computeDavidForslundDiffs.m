useForslund = 0;
if useForslund
    closestTaxaRevData = textscan(fopen('/home/fs01/yw595/closestTaxaForslundRev.txt'),'%s%s','Delimiter','|','HeaderLines',0);
    ForslundData = textscan(fopen('/home/fs01/yw595/MATLAB/SEED/src/ForslundSI3.txt'),'%s%s','Delimiter','|','HeaderLines',0);
    totalBiomass1 = 0; totalBiomass2 = 0;
    totalBiomass1Arr = []; totalBiomass2Arr = [];
    totalFragRank1 = 0; totalFragRank2 = 0;
    totalFragRank1Arr = []; totalFragRank2Arr = [];
    ForslundTaxa1 = ForslundData{1}(strcmp(ForslundData{2},'up'));
    ForslundTaxa2 = ForslundData{1}(strcmp(ForslundData{2},'down'));
    ForslundTaxa1Intersect = {};
    ForslundTaxa2Intersect = {};
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
            ForslundTaxa1Intersect{end+1} = [ForslundTaxa1{i} ' ' ForslundTaxa1{j}];
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
            ForslundTaxa2Intersect{end+1} = [ForslundTaxa2{i} ' ' ForslundTaxa2{j}];
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
    DavidTaxa1Intersect = {};
    DavidTaxa2Intersect = {};
    if useAlt
        DavidTaxa1 = DavidTaxa1Data{1};
        DavidTaxa2 = DavidTaxa2Data{1};
        for i=1:length(DavidTaxa1)
            d = closestTaxaRevData{2}{strcmp(DavidTaxa1{i},closestTaxaRevData{1})};
            totalBiomass1Arr(end+1) = allBiomassRates(strcmp(d,b));
            totalBiomass1 = totalBiomass1 + allBiomassRates(strcmp(d,b));
            for j=1:length(DavidTaxa1)
                d2 = closestTaxaRevData{2}{strcmp(DavidTaxa1{j},closestTaxaRevData{1})};
                DavidTaxa1Intersect{end+1} = [DavidTaxa1{i} ' ' DavidTaxa1{j}];
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
                DavidTaxa2Intersect{end+1} = [DavidTaxa2{i} ' ' DavidTaxa2{j}];
                totalFragRank2 = totalFragRank2 + find(sortIdxsExtra==(find(strcmp(d,b))*length(b)+find(strcmp(d2,b))));
                totalFragRank2Arr(end+1) = find(sortIdxsExtra==(find(strcmp(d,b))*length(b)+find(strcmp(d2,b))));
            end
        end
    else
        for i=1:length(modelNamesShort)
            speciesName = modelNamesShort{i};
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
    if useDavid
        Taxa1 = DavidTaxa1;
        Taxa2 = DavidTaxa2;
        Taxa1Intersect = DavidTaxa1Intersect;
        Taxa2Intersect = DavidTaxa2Intersect;
    elseif useForslund
        Taxa1 = ForslundTaxa1;
        Taxa2 = ForslundTaxa2;
        Taxa1Intersect = ForslundTaxa1Intersect;
        Taxa2Intersect = ForslundTaxa2Intersect;
    end
    [totalBiomassRand1Arr biomassTaxaRand1Arr] = bootstrap(totalBiomass1Arr,50,3,Taxa1);
    [totalBiomassRand2Arr biomassTaxaRand2Arr] = bootstrap(totalBiomass2Arr,50,3,Taxa2);
    [totalFragRankRand1Arr fragrankTaxaRand1Arr] = bootstrap(totalFragRank1Arr,50,floor(length(totalFragRank1Arr)/10),Taxa1Intersect);
    [totalFragRankRand2Arr fragrankTaxaRand2Arr] = bootstrap(totalFragRank2Arr,50,floor(length(totalFragRank2Arr)/10),Taxa2Intersect);
    jaccard1Arr = cellfun(@(x) length(intersect(x,Taxa2))/length(union(x,Taxa2)),biomassTaxaRand1Arr);
    jaccard2Arr = cellfun(@(x) length(intersect(x,Taxa1))/length(union(x,Taxa1)),biomassTaxaRand2Arr);
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
        writeData({jaccard1Arr,totalBiomassRand1Arr},'/home/fs01/yw595/plantVsAnimalJaccardCorr.txt','\t',{'jaccard','biomass'});
    end
    if useForslund
        writeData({[totalBiomassRand1Arr,totalBiomassRand2Arr], groupArr},[transferDir filesep 'diabetesVsNorm.txt'],'\t',{'biomass','group'});
        writeData({[totalFragRankRand1Arr,totalFragRankRand2Arr],groupFragArr},'/home/fs01/yw595/diabetesVsNormFrag.txt','\t',{'fragility','group'});
    end
end