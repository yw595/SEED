function returnModel = mergeModels(origModel,modelTemp,specificRxn)

returnModel = origModel;
for j=1:length(modelTemp.rxns)
    if exist('specificRxn','var')
        matchRxn = strcmp(modelTemp.rxns{j},specificRxn);
    else
        matchRxn = 1;
    end
    if ~any(strcmp(returnModel.rxns,modelTemp.rxns{j})) && matchRxn
        returnModel.rxns{end+1} = modelTemp.rxns{j};
        returnModel.rxnNames{end+1} = modelTemp.rxnNames{j};
        returnModel.subSystems{end+1} = modelTemp.subSystems{j};
        mets = modelTemp.mets(modelTemp.S(:,j)~=0);
        metNames = modelTemp.metNames(modelTemp.S(:,j)~=0);
        coeffs = modelTemp.S(modelTemp.S(:,j)~=0,j);
        for k=1:length(mets)
            if ~any(strcmp(returnModel.mets,mets{k}))
                returnModel.mets{end+1} = mets{k};
                returnModel.metNames{end+1} = metNames{k};
            end
            returnModel.S(strcmp(returnModel.mets,mets{k}),strcmp(returnModel.rxns,modelTemp.rxns{j})) = coeffs(k);
            returnModel.lb(strcmp(returnModel.rxns,modelTemp.rxns{j})) = modelTemp.lb(j);
            returnModel.ub(strcmp(returnModel.rxns,modelTemp.rxns{j})) = modelTemp.ub(j);
        end
    end
end

end