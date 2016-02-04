function [rxnsToECs ECsToRxns] = updateTwoMaps(rxnsToECs, ECsToRxns, rxnName, ECNums)

if ~iscell(ECNums)
    ECNums = {ECNums};
end
if ~isKey(rxnsToECs,rxnName)
    rxnsToECs(rxnName) = ECNums;
else
    temp = rxnsToECs(rxnName);
    for j=1:length(ECNums)
        temp{end+1} = ECNums{j};
    end
    temp = unique(temp);
    rxnsToECs(rxnName) = temp;
end
for j=1:length(ECNums)
    if ~isKey(ECsToRxns,ECNums{j})
        ECsToRxns(ECNums{j}) = {rxnName};
    else
        temp = ECsToRxns(ECNums{j});
        temp{end+1} = rxnName;
        temp = unique(temp);
        ECsToRxns(ECNums{j}) = temp;
    end
end

end