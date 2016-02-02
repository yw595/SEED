function [biomassMets biomassCoeffs] = makeMergedBiomass(allModels,restrictMets)

biomassMets = {}; biomassCoeffs = [];
for i=1:length(allModels);
    modelTemp = allModels{i};
    biomassIdxs = find(cellfun(@(x) ~isempty(regexpi(x,'(biomass|grow)')),modelTemp.rxnNames));
    if ~isempty(biomassIdxs)
        for j=1:length(biomassIdxs)
            biomassMetIdxs = find(modelTemp.S(:,biomassIdxs(j))~=0);
            for k=1:length(biomassMetIdxs)
                if any(strcmp(modelTemp.mets{biomassMetIdxs(k)},restrictMets))
                    biomassMets{end+1} = modelTemp.mets{biomassMetIdxs(k)};
                    biomassCoeffs(end+1) = modelTemp.S(biomassMetIdxs(k),biomassIdxs(j));
                end
            end
        end
    end
end
[biomassMets uniqIdxs ~] = unique(biomassMets);
biomassCoeffs = biomassCoeffs(uniqIdxs);

end