function restrictModel = restrictPrevTE(origModel, restrictFluxDist)

restrictModel = origModel;
for i=1:length(restrictModel.rxns)
    if strcmp(restrictModel.subSystems{i},'Transport') || strcmp(restrictModel.subSystems{i},'Exchange')
        restrictModel.lb(i) = restrictFluxDist(i);
        restrictModel.ub(i) = restrictFluxDist(i);
    end
end