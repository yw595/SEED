function assignedModel = assignSortedBiom(origModel)

assignedModel = origModel;
%from unique results of allBiomNames
sortedBiomNames = {'Biomass','Biomass2','biomass objective function','EX Biomass c','EX Biomass e','Biomass production, carbon limited','Biomass production, nitrogen limited','Model-specific reaction, used to group lipid formation for biomass production (carbon limited)','Model-specific reaction, used to group lipid formation for biomass production (nitrogen limited)','biomass SC5 notrace'};
for k=1:length(sortedBiomNames)
    biomidx = find(strcmp(assignedModel.rxnNames,sortedBiomNames{k}));
    if ~isempty(biomidx)
	if length(biomidx)==1
	    assignedModel = changeObjective(assignedModel,assignedModel.rxns{biomidx});
	else
	    assignedModel = changeObjective(assignedModel,assignedModel.rxns(biomidx),ones(length(biomidx))); 
	end
	break;
    end
end
