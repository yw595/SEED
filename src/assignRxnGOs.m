function assignModel = assignRxnGOs(origModel,keggIDs)

assignModel = origModel;
for i=1:length(assignModel.rxns)
    rxnKEGG = assignModel.rxnKEGGs{i};
    matchIdx = strcmp(keggIDs,['R' rxnKEGG(4:end)]);
    if sum(matchIdx)==1
        assignModel.rxnGOs{i} = goTerms{matchIdx};
    end
end
