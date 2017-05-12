function [biomassRate, shadMet] = examineBiomassFactored(testmodel)

sortedBiomNames = {'Biomass','Biomass2','biomass objective function','EX Biomass c','EX Biomass e','Biomass production, carbon limited','Biomass production, nitrogen limited','Model-specific reaction, used to group lipid formation for biomass production (carbon limited)','Model-specific reaction, used to group lipid formation for biomass production (nitrogen limited)','biomass SC5 notrace'};

for i=1:length(sortedBiomNames)
    biomidx = find(strcmp(testmodel.rxnNames,sortedBiomNames{i}));
    if ~isempty(biomidx)
        %biomidx = find(strcmp(testmodel.rxnNames,sortedBiomNames{j}));
        if length(biomidx)==1
            testmodel = changeObjective(testmodel,testmodel.rxns{biomidx});
        else
            testmodel = changeObjective(testmodel,testmodel.rxns(biomidx),ones(length(biomidx)));
        end
        sol = optimizeCbModel(testmodel);
        [maxShad, maxIdx] = max(abs(sol.y));
        %disp(maxShad)
        %disp(testmodel.modelName);
        %disp(testmodel.metNames{maxIdx});
        %disp(testmodel.metKEGGs{maxIdx});
        shadMet = testmodel.metKEGGs{maxIdx};
        biomassRate = sol.f;
        if strcmp(shadMet,'')
            shadMet = testmodel.metNames{maxIdx};
        end
        break;
    end
end




