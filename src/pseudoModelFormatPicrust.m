function formattedModel = pseudoModelFormatPicrust(origModel)

minModel = makeMinimalFBAModel(origModel);
minModelTrue = subselectModel(minModel,'minModelTrue',minModel.rxns(minModel.lb~=0 | minModel.ub~=0));
[reducedModelArrBig subsystemsAddedArrBig] = makeReducedModelArr(minModelTrue,origModel,1);
formattedModel = reducedModelArrBig{end};
