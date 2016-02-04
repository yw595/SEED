function newModel = makeEmptyModel()

newModel = struct();
newModel.S = [];
newModel.rxns = {};
newModel.mets = {};
newModel.lb = [];
newModel.ub = [];
newModel.rxnNames = {};
newModel.metNames = {};
newModel.subSystems = {};

end