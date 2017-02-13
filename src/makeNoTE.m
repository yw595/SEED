function returnModel = makeNoTE(oldModel,rxnsToECs)

returnModel = makeEmptyModel();
for i=1:length(oldModel.subSystems)
    if ~strcmp(oldModel.subSystems{i},'Transport') && ~strcmp(oldModel.subSystems{i},'Exchange')
        if mod(i,100)==0
            disp(i)
        end
        if isKey(rxnsToECs,oldModel.rxns{i})
            returnModel = mergeModels(returnModel,oldModel,oldModel.rxns{i});
            returnModel.rxnECNums{strcmp(returnModel.rxns,oldModel.rxns{i})} = rxnsToECs(oldModel.rxns{i});
            returnModel = checkModelDims(returnModel);
        else
            disp('Key not found')
            disp(i)
        end
    end
end

end