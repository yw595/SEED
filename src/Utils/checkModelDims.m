function newModel = checkModelDims(oldModel)

newModel = oldModel;
if size(newModel.rxns,1) < size(newModel.rxns,2)
    newModel.rxns = newModel.rxns';
end
if size(newModel.rxnNames,1) < size(newModel.rxnNames,2)
    newModel.rxnNames = newModel.rxnNames';
end
if size(newModel.mets,1) < size(newModel.mets,2)
    newModel.mets = newModel.mets';
end
if size(newModel.metNames,1) < size(newModel.metNames,2)
    newModel.metNames = newModel.metNames';
end
if size(newModel.subSystems,1) < size(newModel.subSystems,2)
    newModel.subSystems = newModel.subSystems';
end
if size(newModel.lb,1) < size(newModel.lb,2)
    newModel.lb = newModel.lb';
end
if size(newModel.ub,1) < size(newModel.ub,2)
    newModel.ub = newModel.ub';
end

end