function [speciesexpr] = simulateFuncFlux(useFBA,modelExprFile,outputDir1,i,j,ucrFolder,AGORAModel,speciesNumberPrefix,doSimulation)
if useFBA==0
    if exist(modelExprFile,'file')
        speciesexpr = importdata(modelExprFile);
        expressionData = speciesexpr;
    else
        return
    end
else
    expressionData = ones(size(AGORAModel.rxns));
end			    
expressionSDs = ones(size(AGORAModel.rxns)).*min(expressionData,1);
expressionIDs = AGORAModel.rxns;
AGORAModel.rxnECNums = {};
for k=1:length(AGORAModel.rxns)
    AGORAModel.rxnECNums{end+1} = {};
end
AGORAModel = fluxModelFunc(AGORAModel);
AGORAModel = assignSortedBiom(AGORAModel);
if doSimulation
    try
	if useFBA
	    picrustFluxes = runFluxMethod(expressionData,expressionIDs,'testfalcon',AGORAModel,'FBA',expressionSDs);
	    writeData({picrustFluxes},[outputDir1 filesep num2str(i) '_' num2str(j) speciesNumberPrefix ucrFolder '_FBA.flux'],'\t');
	else
	    picrustFluxes = runFluxMethod(expressionData,expressionIDs,'testfalcon',AGORAModel,'FALCON',expressionSDs);
	    writeData({picrustFluxes},[outputDir1 filesep num2str(i) '_' num2str(j) speciesNumberPrefix ucrFolder '.flux'],'\t');
	end
    catch
	disp('PICRUST FAIL')
    end
end
