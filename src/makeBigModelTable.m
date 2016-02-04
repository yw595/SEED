configSEED;
outputDir1 = [outputDir filesep 'makeBigModelTable'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1);
end
FI = fopen([inputDir filesep 'reactionsCleaned.csv']);
dataFields = textscan(FI,repmat('%s',1,9),'Delimiter',',');
rxnData = [dataFields{:}]; fclose(FI);
FI = fopen([inputDir filesep 'compoundsCleaned.csv']);
dataFields = textscan(FI,repmat('%s',1,10),'Delimiter',',');
cpdData = [dataFields{:}]; fclose(FI);
rxnsToECsTable = containers.Map; ECsToRxnsTable = containers.Map;
for i=1:length(rxnData)
    rxnName = rxnData{i,1};
    ECCell = strrep(rxnData{i,3},'|',':');
    [regex1 regex2] = regexp(ECCell,'(\d|-)+\.(\d|-)+\.(\d|-)+\.(\d|-)+');
    if ~isempty(regex1)
        ECNums = arrayfun(@(x,y) ECCell(x:y), regex1,regex2, 'UniformOutput',0);
        [rxnsToECsTable ECsToRxnsTable] = updateTwoMaps(rxnsToECsTable, ECsToRxnsTable, rxnName, ECNums);
    end
end

warning('off')
bigModelTable = bigModelAccum;
cpdIDs = cpdData(:,1); cpdNames = cpdData(:,2); cpdAbbrvs = cpdData(:,3);
rxnIDs = rxnData(:,1); rxnNames = rxnData(:,2); equations = rxnData(:,7);
equations = cellfun(@(x) strrep(strrep(x,'(',''),')',''), equations, 'UniformOutput',0);
badIdxs = find(strcmp(equations,'<=> '));
for i=1:length(rxnIDs)
    if mod(i,100)==0
        disp(i)
        disp(length(bigModelTable.rxns));
    end
    if ~any(strcmp(rxnIDs{i},bigModel.rxns)) && isKey(rxnsToECsTable,rxnIDs{i}) && any(ismember(GreenblumEC,rxnsToECsTable(rxnIDs{i}))) && ~any(i==badIdxs)
        [metList coeffList revFlag] = parseRxnFormula(equations{i});
        bigModelTable = addReaction(bigModelTable,rxnIDs{i}, metList,coeffList,revFlag);
    end
end
save([outputDir1 filesep 'makeBigModelTable.mat'],'bigModelTable','rxnsToECsTable','ECsToRxnsTable');