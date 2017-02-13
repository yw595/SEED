function dupRxns = testMultiplicity(model)

dupRxns = {};
uniqueRxns = unique(model.rxnNames);
for i=1:length(uniqueRxns)
    idxs = find(strcmp(model.rxnNames,uniqueRxns{i}));
    multiplicity = length(idxs);
    if multiplicity > 1
        dupRxns{end+1} = model.rxnNames{idxs(1)};
        disp(model.rxns(idxs))
        disp(model.rxnNames(idxs))
        disp(multiplicity)
        for j=1:length(idxs)
            tempECNums = model.rxnECNums(idxs);
            for k=1:length(tempECNums)
                disp(tempECNums{k})
            end
            disp(printRxnEq(model,model.rxns{idxs(j)}));
        end
    end
end

end