function [biomassRate, shadMet, biomassToUse] = examineBiomassFactored(testmodel,biomassToUse)

testmodel = assignSortedBiom(testmodel);

sol = optimizeCbModel(testmodel);
if isfield(sol,'y')
    [maxShad, maxIdx] = max(abs(sol.y));
    shadMet = testmodel.metKEGGs{maxIdx};
    if strcmp(shadMet,'')
        shadMet = testmodel.metNames{maxIdx};
    end
end
biomassRate = sol.f;




