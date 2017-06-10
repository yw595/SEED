%configSEED;
outputDir1 = [outputDir filesep 'examineBiomass'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1);
end

allBiomNames = {};
for i=1:length(modelNames)
    testmodel = modelNamesToModels(modelNames{i});
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
allBiomassDists = {};
for i=1:length(modelNames)
    testmodel = modelNamesToModels(modelNames{i});
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
            %disp(testmodel.modelName);
            disp(testmodel.metNames{maxIdx});
            disp(testmodel.metKEGGs{maxIdx});
            allShadMets{i} = testmodel.metKEGGs{maxIdx};
            allBiomassRates(i) = sol.f;
            allBiomassDists{i} = sol.x;
            if strcmp(allShadMets{i},'')
                allShadMets{i} = testmodel.metNames{maxIdx};
            end
            break;
        end
    end
end

save([outputDir1 filesep 'examineBiomass.mat'],'allBiomassRates','allShadMets','allBiomassDists');

