function rxnsToExpress = mapExpToRxns(ECsToRxns,filename)

FI = fopen(filename);
dataFields = textscan(FI,'%s%s','Delimiter','\t');
fclose(FI);
dataFields = [dataFields{:}];
ECs = keys(ECsToRxns);
rxnsToExpress = containers.Map;
for i=1:length(ECs)
    matchIdx = strcmp(dataFields(:,1),ECs{i});
    matchIdx = find(matchIdx);
    if matchIdx > 0
        rxns = ECsToRxns(ECs{i});
        if ~iscell(rxns)
            rxns = {rxns};
        end
        for j=1:length(rxns)
            if isKey(rxnsToExpress,rxns{j})
                rxnsToExpress(rxns{j}) = rxnsToExpress(rxns{j}) + str2num(dataFields{matchIdx,2});
            else 
                rxnsToExpress(rxns{j}) = str2num(dataFields{matchIdx,2});
            end    
        end
    end
end

end