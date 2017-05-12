allKeys = keys(modelNamesToModels);
if 1
FI = fopen([inputDir filesep 'compoundsCleaned.csv']);
dataFields = textscan(FI,repmat('%s',1,10),'Delimiter',',');
cpdData = [dataFields{:}]; fclose(FI);
cpdIDs = cpdData(:,1); cpdNames = cpdData(:,2); cpdAbbrvs = cpdData(:,3); cpdKEGGs = cpdData(:,5);

for i=1:length(allKeys)
    i
    testmodel = modelNamesToModels(allKeys{i});
    for j=1:length(testmodel.mets)
        reconcmet = testmodel.mets{j};
        if ~isempty(regexp(reconcmet,'\['))
            reconcmet = reconcmet(1:end-3);
        end
        matchIdx = find(strcmp(reconcmet,cpdIDs));
        if ~isempty(matchIdx)
            matchIdx = matchIdx(1);
            testmodel.metKEGGs{j} = cpdKEGGs{matchIdx};
        else
            testmodel.metKEGGs{j} = '';
        end
    end
    modelNamesToModels(allKeys{i}) = testmodel;
end
end
allBiomNames = {};
for i=1:length(allKeys)
    testmodel = modelNamesToModels(allKeys{i});
    biomidx = find(cellfun(@(x) length(regexpi(x,'BIOM')),testmodel.rxnNames));
    %allKeys{i}
    %testmodel.rxnNames(biomidx)
    for j=1:length(biomidx)
        allBiomNames{end+1} = testmodel.rxnNames{biomidx(j)};
    end
end

%from unique results of allBiomNames
sortedBiomNames = {'Biomass','Biomass2','biomass objective function','EX Biomass c','EX Biomass e','Biomass production, carbon limited','Biomass production, nitrogen limited','Model-specific reaction, used to group lipid formation for biomass production (carbon limited)','Model-specific reaction, used to group lipid formation for biomass production (nitrogen limited)','biomass SC5 notrace'};

allShadMets = {};
allBiomassRates = [];
for i=1:length(allKeys)
    testmodel = modelNamesToModels(allKeys{i});
    for j=1:length(sortedBiomNames)
        biomidx = find(strcmp(testmodel.rxnNames,sortedBiomNames{j}));
        if ~isempty(biomidx)
            %biomidx = find(strcmp(testmodel.rxnNames,sortedBiomNames{j}));
            if length(biomidx)==1
                testmodel = changeObjective(testmodel,testmodel.rxns{biomidx});
            else
                testmodel = changeObjective(testmodel,testmodel.rxns(biomidx),ones(length(biomidx)));
            end
            sol = optimizeCbModel(testmodel);
            [maxShad, maxIdx] = max(abs(sol.y));
            disp(maxShad)
            disp(testmodel.modelName);
            disp(testmodel.metNames{maxIdx});
            disp(testmodel.metKEGGs{maxIdx});
            allShadMets{i} = testmodel.metKEGGs{maxIdx};
            allBiomassRates(i) = sol.f;
            if strcmp(allShadMets{i},'')
                allShadMets{i} = testmodel.metNames{maxIdx};
            end
            break;
        end
    end
end

numShadSpecies = zeros(length(uniqAllShadMets),1);
xArr = []; yArr = []; speciesArr = {}; metArr = {}; presabsArr = {};
uniqAllShadMets = unique(allShadMets)
for i=1:length(allKeys)
    for j=1:length(uniqAllShadMets)
        xArr(end+1) = j;
        yArr(end+1) = i;
        speciesArr{end+1} = modelNamesToModels(allKeys{i}).modelName;
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