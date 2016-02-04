configSEED;

bigModelTableFlux = addMustEx(bigModelTable);
[bigModelTableFlux, ~, addedModel] = mergeSmallest(bigModelTableFlux,modelNamesToModelsValues,1);
[biomassMets biomassCoeffs] = makeMergedBiomass({addedModel},bigModelTableFlux.mets);
bigModelTableFlux = addReaction(bigModelTableFlux,'BIOMASS', {bigModelTableFlux.mets{end}},[-1],1);
bigModelTableFlux.c = zeros(length(bigModelTableFlux.rxns),1);
bigModelTableFlux = changeObjective(bigModelTableFlux,'BIOMASS');
expressionIDs = rxnsToExpressNormKeys;
for i=1:length(expressionIDs)
    expressionData(i) = rxnsToExpressNorm(expressionIDs{i});
end
expressionSDs = ones(size(expressionData));
for i=1:length(bigModelTableFlux.rxns)
    bigModelTableFlux.rxns{i} = strrep(strrep(bigModelTableFlux.rxns{i},')',''),'(','');
end
bigModelTableFlux.genes = bigModelTableFlux.rxns;
bigModelTableFlux.grRules = bigModelTableFlux.rxns;
bigModelTableFlux.rules = bigModelTableFlux.rxns;
bigModelTableFlux.rxnGeneMat = eye(length(bigModelTableFlux.rxns));
bigModelTableFlux = changeObjective(bigModelTableFlux,'BIOMASS');
bigModelTableFlux.description = 'bigModel';
falconfluxes = runFluxMethod(expressionData,expressionIDs,'test',bigModelTableFlux,'FALCON',expressionSDs);