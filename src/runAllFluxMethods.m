function normFluxes = runAllFluxMethods(model,expressionData,expressionSDs,expressionIDs)

    fbasol = optimizeCbModel(model);
    fbafluxes = fbasol.x;
    efluxes = runFluxMethod(expressionData,expressionIDs,'testeflux',model,'EFlux');
    falconfluxes = runFluxMethod(expressionData,expressionIDs,'testfalcon',model,'FALCON',expressionSDs);
    gxfbafluxes = runFluxMethod(expressionData,expressionIDs,'testgxfba',model,'GXFBA');
    normFluxes(:,1) = fbafluxes;
    normFluxes(:,2) = efluxes;
    normFluxes(:,3) = falconfluxes;
    normFluxes(:,4) = gxfbafluxes;
    for i=1:4
        normFluxes(:,i) = normFluxes(:,i)/max(abs(normFluxes(:,i)));
    end

    modelRELATCH = model;
    expressionIDsRELATCH = expressionIDs;
    for i=1:length(expressionIDs)
        %expressionIDsRELATCH{i} = num2str(i);
    end
    %uniqueRules = unique(modelRELATCH.rules);
    modelRELATCH.genes = unique(modelRELATCH.genes);
    for i=1:length(modelRELATCH.genes)
        %modelRELATCH.genes{i} = num2str(find(strcmp(modelRELATCH.genes{i},expressionIDs)));
        modelRELATCH.rules{i} = [ '(x(' num2str(find(strcmp(modelRELATCH.rules{i},modelRELATCH.genes))) '))'];  
        if strcmp(modelRELATCH.rules{i},'(x())')
            modelRELATCH.rules{i} = '';
        end
        %modelRELATCH.grRules{i} = num2str(find(strcmp(modelRELATCH.grRules{i},expressionIDs)));
    end
    %relatchfluxes = runFluxMethod(expressionData,expressionIDs,'testrelatch',model,'RELATCH');

    %madefluxes = runFluxMethod(expressionData,expressionIDs,'testmade',model,'MADE');
end