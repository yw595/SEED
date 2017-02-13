function returnModel = makeUnique(oldModel,rxnsToECs)

returnModel = makeEmptyModel();
seenECs = {};
for i=1:length(oldModel.subSystems)
    currentECs = oldModel.rxnECNums{i};
    notInSeenECs = 0;
    for j=1:length(currentECs)
        if ~any(strcmp(currentECs{j},seenECs))
            notInSeenECs = 1;
            seenECs{end+1} = currentECs{j};
        end
    end
    if notInSeenECs
        i
        returnModel = mergeModels(returnModel,oldModel,oldModel.rxns{i});
        returnModel.rxnECNums{strcmp(returnModel.rxns,oldModel.rxns{i})} = rxnsToECs(oldModel.rxns{i});
        returnModel = checkModelDims(returnModel);
    end
end

end