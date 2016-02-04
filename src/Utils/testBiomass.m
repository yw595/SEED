function [sols, bigModelBio] = testBiomass(bigModel,biomassMets,biomassCoeffs,cumulative)

sols = {};
for i=1:length(biomassMets)
    disp(i)
    if i==1 | ~cumulative
        bigModelBio = bigModel;
        bigModelBio.rxns{end+1} = 'BIOMASS';
        bigModelBio.rxnNames{end+1} = 'BIOMASS';
        bigModelBio.subSystems{end+1} = 'BIOMASS';
        biomassIdx = length(bigModelBio.rxns);
    end
    bigModelBio.S(ismember(bigModelBio.mets,biomassMets{i}),biomassIdx) = biomassCoeffs(i);
    bigModelBio.lb(biomassIdx)=-1000;bigModelBio.ub(biomassIdx)=1000;
    bigModelBio.c = zeros(length(bigModelBio.rxns),1);
    bigModelBio = changeObjective(bigModelBio,'BIOMASS');
    sols{end+1} = optimizeCbModel(bigModelBio);
end

end