function [mergedModel, minRxnsAdded, addedModel] = mergeSmallest(origModel,modelsTemp,sortNum)

if ~exist('sortNum','var')
    sortNum=1;
end
mergedModel = origModel;
for i=1:length(modelsTemp)
    modelTemp = modelsTemp{i};
    diffnums(i) = length(setdiff(modelTemp.rxns,origModel.rxns));
end
diffnums(diffnums==0) = max(diffnums)+1;
[sortdiffnums sortIdxs] = sort(diffnums,'ascend');
minRxnsAdded = sortdiffnums(sortNum);
minIdx = sortIdxs(sortNum);
addedModel = modelsTemp{minIdx};
mergedModel = mergeModels(origModel,modelsTemp{minIdx});

end