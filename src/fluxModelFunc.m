function fluxAssignModel = fluxModelFunc(origModel)

fluxAssignModel = origModel;
for i=1:length(fluxAssignModel.rxns)
    fluxAssignModel.rxns{i} = strrep(strrep(fluxAssignModel.rxns{i},')',''),'(','');
end
fluxAssignModel.genes = fluxAssignModel.rxns;
fluxAssignModel.grRules = fluxAssignModel.rxns;
fluxAssignModel.rules = fluxAssignModel.rxns;
fluxAssignModel.rxnGeneMat = speye(length(fluxAssignModel.rxns));
