function picrustFluxes = pseudoMapPicrust(model,useRandNorm,geneExpVal,biomassName)

%useRandNorm = 0;
if useRandNorm
    expressionData = geneExpVal*abs(normrnd(0,1,length(model.rxns),1));
else
    expressionData = geneExpVal*ones(length(model.rxns),1);
end
expressionSDs = ones(size(model.rxns));
expressionIDs = model.rxns;
if exist('biomassName','var')
    picrustFluxes = runFluxMethod(expressionData,expressionIDs,'testfalcon',model,'FALCON',expressionSDs,biomassName);
else
    picrustFluxes = runFluxMethod(expressionData,expressionIDs,'testfalcon',model,'FALCON',expressionSDs);
end
