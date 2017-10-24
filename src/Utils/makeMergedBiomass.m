function [biomassMets biomassCoeffs] = makeMergedBiomass(allModels,restrictMets)

biomassMets = {}; biomassCoeffs = [];
for i=1:length(allModels);
    modelTemp = allModels{i};
    biomassIdxs = find(cellfun(@(x) ~isempty(regexpi(x,'(biomass|grow)')),modelTemp.rxnNames));
    if ~isempty(biomassIdxs)
        for j=1:length(biomassIdxs)
            biomassMetIdxs = find(modelTemp.S(:,biomassIdxs(j))~=0);
            %disp(modelTemp.mets(biomassMetIdxs))
            %disp(restrictMets)
            for k=1:length(biomassMetIdxs)
		candMet = modelTemp.mets{biomassMetIdxs(k)};
                candMet2 = '';
                if any(strcmp(candMet,restrictMets))
		    candMet2 = candMet;
                end
                if any(strcmp(candMet(1:end-3),restrictMets))
		    candMet2 = candMet(1:end-3);
                end
		if ~strcmp(candMet2,'')
                    biomassMets{end+1} = candMet2;%modelTemp.mets{biomassMetIdxs(k)};
                    biomassCoeffs(end+1) = modelTemp.S(biomassMetIdxs(k),biomassIdxs(j));
                end
	    end
        end
    end
end
[biomassMets uniqIdxs ~] = unique(biomassMets);
biomassCoeffs = biomassCoeffs(uniqIdxs);

end
