allTwoModels = {};
allTwoBiomasses = [];
allModelKeys = keys(modelNamesToModels);
parfor z=1:length(allModelKeys)^2
    i = floor(z/length(allModelKeys))+1;
    j = mod(z,length(allModelKeys))+1;
    %for i=1:length(allModelKeys)
    %parfor j=1:length(allModelKeys)
        initCobraToolbox;
        if i~=j
            i
            j
            newTwoModel = makeEmptyModel();
            ithModel = modelNamesToModels(allModelKeys{i});
            ithModel.rxnECNums = {};
            jthModel = modelNamesToModels(allModelKeys{j});
            jthModel.rxnECNums = {};
            newTwoModel = mergeModels(ithModel,jthModel);
            newTwoModel.metKEGGs = {};
            for k=1:length(newTwoModel.mets)
                reconcmet = newTwoModel.mets{k};
                if ~isempty(regexp(reconcmet,'\['))
                    reconcmet = reconcmet(1:end-3);
                end
                matchIdx = find(strcmp(reconcmet,cpdIDs));
                if ~isempty(matchIdx)
                    matchIdx = matchIdx(1);
                    newTwoModel.metKEGGs{k} = cpdKEGGs{matchIdx};
                else
                    newTwoModel.metKEGGs{k} = '';
                end
            end
            newTwoModel = checkModelDims(newTwoModel);
            allTwoBiomasses(z) = examineBiomassFactored(newTwoModel);
            %end
    end
end


for i=1:length(allTwoModels)
    for j=1:length(allTwoModels)
        if i~=j
        end
    end
end