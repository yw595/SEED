configSEED;
IGNOREDASH = 1;
outputDir1 = [outputDir filesep 'makeBigModelTable'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1);
end

rxnsToECsTable = containers.Map; ECsToRxnsTable = containers.Map;
for i=1:length(rxnData)
    rxnName = rxnData{i,1};
    ECCell = strrep(rxnData{i,3},'|',':');
    if ~IGNOREDASH
        [regex1 regex2] = regexp(ECCell,'(\d|-)+\.(\d|-)+\.(\d|-)+\.(\d|-)+');
    else
        [regex1 regex2] = regexp(ECCell,'(\d)+\.(\d)+\.(\d)+\.(\d)+');
    end
    if ~isempty(regex1)
        ECNums = arrayfun(@(x,y) ECCell(x:y), regex1,regex2, 'UniformOutput',0);
        [rxnsToECsTable ECsToRxnsTable] = updateTwoMaps(rxnsToECsTable, ECsToRxnsTable, rxnName, ECNums);
    end
end

warning('off')
bigModelTable = bigModelAccum;
badIdxs = find(strcmp(equations,'<=> '));
for i=1:length(rxnIDs)
    if mod(i,100)==0
        disp(i)
        disp(length(bigModelTable.rxns));
    end
    if ~any(strcmp(rxnIDs{i},bigModelTable.rxns)) && (~isempty(regexp(rxnNames{i},'EX')) || ~isempty(regexpi(rxnNames{i},'transport')))
        [metList coeffList revFlag] = parseRxnFormula(equations{i});
        % for j=1:length(metList)
        %     bracketIdx = regexp(metList{j},'[');
        %     if ~isempty(bracketIdx)
        %         metList{j} = metList{j}(1:bracketIdx-1);
        %     end
        % end
        bigModelTable = addReaction(bigModelTable,rxnIDs{i}, metList,coeffList,revFlag);
        if ~isempty(find(strcmp(bigModelTable.rxns,rxnIDs{i})))
            bigModelTable.rxnECNums{strcmp(bigModelTable.rxns,rxnIDs{i})} = {};
        end
        if ~isempty(regexp(rxnNames{i},'EX'))
            bigModelTable.subSystems{end} = 'Exchange';
        else
            bigModelTable.subSystems{end} = 'Transport';
        end
    end
    if ~any(strcmp(rxnIDs{i},bigModelTable.rxns)) && isKey(rxnsToECsTable,rxnIDs{i}) && any(ismember(GreenblumEC,rxnsToECsTable(rxnIDs{i}))) && ~any(i==badIdxs)
        [metList coeffList revFlag] = parseRxnFormula(equations{i});
        % for j=1:length(metList)
        %     bracketIdx = regexp(metList{j},'[');
        %     if ~isempty(bracketIdx)
        %         metList{j} = metList{j}(1:bracketIdx-1);
        %     end
        % end
        bigModelTable = addReaction(bigModelTable,rxnIDs{i}, metList,coeffList,revFlag);
        bigModelTable.rxnECNums{strcmp(bigModelTable.rxns,rxnIDs{i})} = rxnsToECsTable(rxnIDs{i});
    end
end

bigModel = addMustEx(bigModelTable);
bigModel.metKEGGs = {};
for i=1:length(bigModel.mets)
    reconcmet = bigModel.mets{i};
    if ~isempty(regexp(reconcmet,'\['))
        reconcmet = reconcmet(1:end-3);
    end
    matchIdx = find(strcmp(reconcmet,cpdIDs));
    if ~isempty(matchIdx)
        matchIdx = matchIdx(1);
        bigModel.metKEGGs{i} = cpdKEGGs{matchIdx};
    else
        bigModel.metKEGGs{i} = '';
    end
end
bigModel.metKEGGs = bigModel.metKEGGs';
modelNamesToModelsValues = values(modelNamesToModels);
testModel = modelNamesToModelsValues{1};
bigModelAdded = mergeModels(bigModel,testModel);
save([outputDir1 filesep 'makeBigModelTable.mat'],'bigModelTable','rxnsToECsTable','ECsToRxnsTable','bigModel','bigModelAdded','testModel');