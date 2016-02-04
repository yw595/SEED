function writeOutModel(model,outputFile)

KEGGmets = model.mets;
KEGGmets = cellfun(@(x) strrep(x,'cpd','C'),KEGGmets,'UniformOutput',0);
for i=1:length(KEGGmets)
    braceIdx = regexp(KEGGmets{i},'[');
    currKEGGmet = KEGGmets{i};
    if ~isempty(braceIdx)
        KEGGmets{i} = currKEGGmet(1:braceIdx-1);
    end
end
metField = {}; stoichField = []; rxnField = {}; subSystemField = {}; KEGGField = {}; longMetField = {}; longRxnField = {};
for i=1:length(model.rxns)
    metIdxs = find(model.S(:,i));
    for j=1:length(metIdxs)
        rxnField{end+1} = model.rxns{i};
        stoichField(end+1) = full(model.S(metIdxs(j),i));
        metField{end+1} = model.mets{metIdxs(j)};
        subSystemField{end+1} = model.subSystems{i};
        KEGGField{end+1} = KEGGmets{metIdxs(j)};
    end
end
writeData({rxnField,stoichField,metField,subSystemField,KEGGField},outputFile,'\t');
end