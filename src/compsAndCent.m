function [ecnums comps cent] = compsAndCent(model)

connMatrixTable = makeConnMatrix(model);
rxnMatrix = makeRxnConnMatrix(model,1);
centsRxnTable = betweenness_centrality(sparse(rxnMatrix));
centsRxnTable = centsRxnTable*2/((length(centsRxnTable)-1)*(length(centsRxnTable)-2));
ECsToCents = containers.Map;
for i=1:length(model.rxnNames)
    rxnECs = model.rxnECNums{i};
    for j=1:length(rxnECs)
        if isKey(ECsToCents,rxnECs{j})
            temp = ECsToCents(rxnECs{j});
            temp(end+1) = centsRxnTable(i);
            ECsToCents(rxnECs{j}) = temp;
        else
            ECsToCents(rxnECs{j}) = [centsRxnTable(i)];
        end
    end
end
ECsToCentsKeys = keys(ECsToCents);
disp(length(ECsToCentsKeys))
correspondingCents = [];
for i=1:length(ECsToCentsKeys)
    correspondingCents(i) = mean(ECsToCents(ECsToCentsKeys{i}));
end

[a b] = components(sparse(rxnMatrix));
comps = max(b);
cent = correspondingCents;
ecnums = ECsToCentsKeys;
%writeData({ECsToCentsKeys,correspondingCents},[baseDir filesep 'ECsToCents.txt'],'\t');

end