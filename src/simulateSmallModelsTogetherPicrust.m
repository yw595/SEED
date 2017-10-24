configSEED;
outputDir1 = [outputDir filesep 'simulateSmallModelsTogetherPicrust'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1);
end

if 0
togetherModel = modelNamesToModels(modelNames{1});
togetherModel.description = 'together';
togetherModel.rxnECNums = {};
for k=1:length(togetherModel.rxns)
	togetherModel.rxnECNums{end+1} = {};
end
togetherModel.rxns = cellfun(@(x) [x '_1'],togetherModel.rxns,'UniformOutput',0);
for i=2:length(modelNames)
	i
	ithModel = modelNamesToModels(modelNames{i});
ithModel.description = ithModel.modelName;
ithModel.rxnECNums = {};
for k=1:length(ithModel.rxns)
	ithModel.rxnECNums{end+1} = {};
end
ithModel.rxns = cellfun(@(x) [x '_' num2str(i)],ithModel.rxns,'UniformOutput',0);
togetherModel = mergeModels(togetherModel,ithModel,'',0);
end
end

if 0
togetherModel2 = fluxModelFunc(togetherModel);
togetherModel2 = assignSortedBiom(togetherModel2);
cIdxs = find(togetherModel2.c);
if length(cIdxs)==1
    pseudoFluxDistTogether1 = pseudoMapPicrust(togetherModel2,1,1,togetherModel2.rxnNames{cIdxs});
    pseudoFluxDistTogether2 = pseudoMapPicrust(togetherModel2,0,1,togetherModel2.rxnNames{cIdxs});
else
    pseudoFluxDistTogether1 = pseudoMapPicrust(togetherModel2,1,1,togetherModel2.rxnNames{cIdxs(1)});
    pseudoFluxDistTogether2 = pseudoMapPicrust(togetherModel2,0,1,togetherModel2.rxnNames{cIdxs(1)});
end
end

if 0
pseudoFluxDistArr1 = {};
pseudoFluxDistArr2 = {};
subLabelsArr = {};
for i=1:length(modelNames)
    try
    model1 = modelNamesToModels(modelNames{i});
    model1.description = model1.modelName;
    model1.rxnECNums = {};
    for k=1:length(model1.rxns)
	model1.rxnECNums{end+1} = {};
    end
    model1.rxns = cellfun(@(x) [x '_1'],model1.rxns,'UniformOutput',0);
    model1 = fluxModelFunc(model1);
    model1 = assignSortedBiom(model1);
    cIdxs = find(model1.c);
    if length(cIdxs)==1
  	pseudoFluxDist1 = pseudoMapPicrust(model1,1,1,model1.rxnNames{cIdxs});
        pseudoFluxDist2 = pseudoMapPicrust(model1,0,1,model1.rxnNames{cIdxs});
    else
	pseudoFluxDist1 = pseudoMapPicrust(model1,1,1,model1.rxnNames{cIdxs(1)});
        pseudoFluxDist2 = pseudoMapPicrust(model1,0,1,model1.rxnNames{cIdxs(1)});
    end
    pseudoFluxDist1
    pseudoFluxDistArr1{i} = pseudoFluxDist1;
    pseudoFluxDistArr2{i} = pseudoFluxDist2;
    [subLabels,~,~] = segmentFluxBySubsystem(model1,pseudoFluxDist1,1,pseudoFluxDist2,1);
    subLabelsArr{i} = subLabels;
    catch
        disp('CARCH')
    end
end
end
      
[totalSubsNames totalSubsRanks] = makeTotalSubsRanks(subLabelsArr);
[totalSubsRanks,sortIdxs] = sort(totalSubsRanks,'descend');
totalSubsNames = totalSubsNames(sortIdxs);
totalSubsNames = addIdxStrings(totalSubsNames);
writeData({totalSubsNames,totalSubsRanks},'/mnt/vdb/home/ubuntu2/modelTogetherTotalSubs.txt','\t',{'subname','subrank'});
