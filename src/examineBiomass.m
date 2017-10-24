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

allShadMets = {};
allBiomassRates = [];
allBiomassDists = {};
for i=1:length(modelNames)
    testmodel = modelNamesToModels(modelNames{i});
    testmodel = assignSortedBiom(testmodel);
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

save([outputDir1 filesep 'examineBiomass.mat'],'allBiomassRates','allShadMets','allBiomassDists');

