function newModel = addMustEx(oldModel)

newModel = oldModel;
for i=1:length(newModel.mets)
    rxnIdxs = find(newModel.S(i,:)~=0);
    if (all(newModel.S(i,rxnIdxs)' < 0 & newModel.lb(rxnIdxs) >= 0)) || (all(newModel.S(i,rxnIdxs)' > 0 & newModel.ub(rxnIdxs) <= 0)) 
        cleanedMetName = strrep(strrep(newModel.mets{i},']',')'),'[','(');
        newModel.rxns{end+1} = ['MUST_EX_' cleanedMetName];
        newModel.rxnNames{end+1} = ['MUST_EX_' cleanedMetName];
        newModel.subSystems{end+1} = 'MUST_EX';
        newModel.S(i,end+1) = 1;
        newModel.lb(end+1) = -1000;
        newModel.ub(end+1) = 1000;
        newModel.rev(end+1) = 1;
        newModel.rxnECNums{end+1} = {};
        newModel.c(end+1) = 0;
    end
end

end