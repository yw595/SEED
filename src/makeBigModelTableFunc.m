function expandedModel = makeBigModelTableFunc(origModel,rxnData,equations,rxnIDs,rxnNames,GreenblumEC,cpdIDs,cpdKEGGs)

  IGNOREDASH = 1;
  
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
expandedModel = origModel;
badIdxs = find(strcmp(equations,'<=> '));
for i=1:length(rxnIDs)
    if mod(i,100)==0
        disp(i)
        disp(length(expandedModel.rxns));
    end
    if ~any(strcmp(rxnIDs{i},expandedModel.rxns)) && (~isempty(regexp(rxnNames{i},'EX')) || ~isempty(regexpi(rxnNames{i},'transport')))
        [metList coeffList revFlag] = parseRxnFormula(equations{i});
        % for j=1:length(metList)
        %     bracketIdx = regexp(metList{j},'[');
        %     if ~isempty(bracketIdx)
        %         metList{j} = metList{j}(1:bracketIdx-1);
        %     end
        % end
        expandedModel = addReaction(expandedModel,rxnIDs{i}, metList,coeffList,revFlag);
        if ~isempty(find(strcmp(expandedModel.rxns,rxnIDs{i})))
            expandedModel.rxnECNums{strcmp(expandedModel.rxns,rxnIDs{i})} = {};
        end
        if ~isempty(regexp(rxnNames{i},'EX'))
            expandedModel.subSystems{end} = 'Exchange';
        else
            expandedModel.subSystems{end} = 'Transport';
        end
    end
    if ~any(strcmp(rxnIDs{i},expandedModel.rxns)) && isKey(rxnsToECsTable,rxnIDs{i}) && any(ismember(GreenblumEC,rxnsToECsTable(rxnIDs{i}))) && ~any(i==badIdxs)
        [metList coeffList revFlag] = parseRxnFormula(equations{i});
        % for j=1:length(metList)
        %     bracketIdx = regexp(metList{j},'[');
        %     if ~isempty(bracketIdx)
        %         metList{j} = metList{j}(1:bracketIdx-1);
        %     end
        % end
        expandedModel = addReaction(expandedModel,rxnIDs{i}, metList,coeffList,revFlag);
        expandedModel.rxnECNums{strcmp(expandedModel.rxns,rxnIDs{i})} = rxnsToECsTable(rxnIDs{i});
    end
end

expandedModel = addMustEx(expandedModel);
expandedModel.metKEGGs = {};
for i=1:length(expandedModel.mets)
    reconcmet = expandedModel.mets{i};
    if ~isempty(regexp(reconcmet,'\['))
        reconcmet = reconcmet(1:end-3);
    end
    matchIdx = find(strcmp(reconcmet,cpdIDs));
    if ~isempty(matchIdx)
        matchIdx = matchIdx(1);
        expandedModel.metKEGGs{i} = cpdKEGGs{matchIdx};
    else
        expandedModel.metKEGGs{i} = '';
    end
end
expandedModel.metKEGGs = expandedModel.metKEGGs';
