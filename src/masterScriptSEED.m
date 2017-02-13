configSEED;
load([outputDir filesep 'makeBigModelAccum' filesep 'makeBigModelAccum.mat']);
%modelNames = keys(modelNamesToModels);
%modelSizes = cellfun(@(x) length(modelNamesToModels(x).rxns), modelNames);
%[~, sortIdxs] = sort(modelSizes);
%modelNames = modelNames(sortIdxs);
modelnamesArr = {}; normArr = []; transArr = []; exchArr = []; unclassArr = [];
for i=1:length(bigModelsAccum)
    ithNum = num2str(i);
    while length(ithNum) < 3
        ithNum = ['0' ithNum];
    end
    modelnamesArr{i} = ithNum;
    transArr(i) = sum(strcmp(bigModelsAccum{i}.subSystems, 'Transport'));
    exchArr(i) = sum(strcmp(bigModelsAccum{i}.subSystems,'Exchange'));
    unclassArr(i) = sum(strcmp(bigModelsAccum{i}.subSystems,''));
    normArr(i) = sum(~strcmp(bigModelsAccum{i}.subSystems,'Transport') & ~strcmp(bigModelsAccum{i}, 'Exchange') & ~strcmp(bigModelsAccum{i}.subSystems,''));
end
writeData({modelnamesArr,normArr,transArr,exchArr,unclassArr},'/home/fs01/yw595/subsystemsAdded.txt','\t',{'modelname','norm','trans','exch','unclass'});
    



