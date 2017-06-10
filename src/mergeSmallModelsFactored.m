function newAllTwo = mergeSmallModelsFactored(ithIdx,jthIdx,modelNamesToModels,cpdKEGGs,cpdIDs)
newTwoModel = makeEmptyModel();
ithModel = modelNamesToModels(modelNames{ithIdx});
jthModel = modelNamesToModels(modelNames{jthIdx});
newTwoModel = mergeModels(ithModel,jthModel);

newTwoModel = checkModelDims(newTwoModel);
[newIthBiom,~,biomassToUseI] = examineBiomassFactored(ithModel);
[newJthBiom,~,biomassToUseJ] = examineBiomassFactored(jthModel);
[minBiom, minIdx] = max([newIthBiom,newJthBiom]);
biomassToUseI
biomassToUseJ
if minIdx==1
    newAllTwo = examineBiomassFactored(newTwoModel,biomassToUseI);
else
    newAllTwo = examineBiomassFactored(newTwoModel,biomassToUseJ);
end
newAllTwo - minBiom
end