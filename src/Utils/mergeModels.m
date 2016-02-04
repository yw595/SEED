function returnModel = mergeModels(origModel,modelTemp,specificRxn)

returnModel = origModel;
count = 0;
for j=1:length(modelTemp.rxns)
    if exist('specificRxn','var')
        matchRxn = strcmp(modelTemp.rxns{j},specificRxn);
    else
        matchRxn = 1;
    end
    if ~any(strcmp(returnModel.rxns,modelTemp.rxns{j})) && matchRxn
        count = count+1;
        %disp(count);
        returnModel.rxns{end+1} = modelTemp.rxns{j};
        returnModel.rxnNames{end+1} = modelTemp.rxnNames{j};
        returnModel.subSystems{end+1} = modelTemp.subSystems{j};
        returnModel.lb(strcmp(returnModel.rxns,modelTemp.rxns{j})) = modelTemp.lb(j);
        returnModel.ub(strcmp(returnModel.rxns,modelTemp.rxns{j})) = modelTemp.ub(j);
        returnModel.rev(strcmp(returnModel.rxns,modelTemp.rxns{j})) = modelTemp.rev(j);
        if ~isempty(regexp(modelTemp.rxnNames{j},'EX'))
            returnModel.subSystems{end} = 'Exchange';
        elseif ~isempty(regexpi(modelTemp.rxnNames{j},'transport'))
            returnModel.subSystems{end} = 'Transport';
        end
        mets = modelTemp.mets(modelTemp.S(:,j)~=0);
        metNames = modelTemp.metNames(modelTemp.S(:,j)~=0);
        coeffs = modelTemp.S(modelTemp.S(:,j)~=0,j);
        for k=1:length(mets)
            if ~any(strcmp(returnModel.mets,mets{k}))
                returnModel.mets{end+1} = mets{k};
                returnModel.metNames{end+1} = metNames{k};
            end
            returnModel.S(strcmp(returnModel.mets,mets{k}),strcmp(returnModel.rxns,modelTemp.rxns{j})) = coeffs(k);
        end
    end
end

end