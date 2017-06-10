function [biomassRate, shadMet, biomassToUse] = examineBiomassFactored(testmodel,biomassToUse)

sortedBiomNames = {'Biomass','Biomass2','biomass objective function','EX Biomass c','EX Biomass e','Biomass production, carbon limited','Biomass production, nitrogen limited','Model-specific reaction, used to group lipid formation for biomass production (carbon limited)','Model-specific reaction, used to group lipid formation for biomass production (nitrogen limited)','biomass SC5 notrace'};

if ~exist('biomassToUse','var')
    for i=1:length(sortedBiomNames)
        biomidx = find(strcmp(testmodel.rxnNames,sortedBiomNames{i}));
        if ~isempty(biomidx)
            biomassToUse = i;
            break;
        end
    end
else
    biomidx = find(strcmp(testmodel.rxnNames,sortedBiomNames{biomassToUse}));
end

if length(biomidx)==1
    testmodel = changeObjective(testmodel,testmodel.rxns{biomidx});
else
    testmodel = changeObjective(testmodel,testmodel.rxns(biomidx),ones(length(biomidx)));
end
sol = optimizeCbModel(testmodel);
if isfield(sol,'y')
    [maxShad, maxIdx] = max(abs(sol.y));
    shadMet = testmodel.metKEGGs{maxIdx};
    if strcmp(shadMet,'')
        shadMet = testmodel.metNames{maxIdx};
    end
end
%disp(maxShad)
%disp(testmodel.modelName);
%disp(testmodel.metNames{maxIdx});
%disp(testmodel.metKEGGs{maxIdx});
biomassRate = sol.f;




