%bigModelReconc = bigModelAdded;
%reconcSubArr = unique(bigModelReconc.subSystems);
%reconcSubMatchArr = zeros(length(reconcSubArr),1);

if 0

bigModelReconc = bigModelAdded;
bigModelReconc.prevmets = bigModelReconc.mets;
bigModelReconc.prevmetNames = bigModelReconc.metNames;
bigModelReconc.prevrxns = bigModelReconc.rxns;
bigModelReconc.prevrxnNames = bigModelReconc.rxnNames;

origRecon2MetKEGGs = {};
origRecon2MetKEGGsHash = {};
for i=1:length(origRecon2.rxns)
    origRecon2MetKEGGs{i} = origRecon2.metKeggID(origRecon2.S(:,i)~=0);
    hashArr = [];
    for j=1:length(origRecon2MetKEGGs{i})
        if length(origRecon2MetKEGGs{i}{j}) > 0
            hashArr(j) = str2num(origRecon2MetKEGGs{i}{j}(2:end));
        else
            hashArr(j) = 0;
        end
    end
    origRecon2MetKEGGsHash{i} = hashArr;
end
matchRxnNum = 0;
for i=1:length(bigModelReconc.rxns)
    i
    searchRxnNames = cellfun(@(x) x(1:regexp(x,'_')-1), bigModelReconc.rxnNames,'UniformOutput',0);
    ithMetKEGGs = bigModelReconc.metKEGGs(bigModelReconc.S(:,i)~=0);
    ithHashArr = [];
    for j=1:length(ithMetKEGGs)
        if length(ithMetKEGGs{j}) == 6
            ithHashArr(j) = str2num(ithMetKEGGs{j}(2:end));
        else
            ithHashArr(j) = 0;
        end
    end
    ithHashArr
    commonKEGGsArrHash = cellfun(@(x) intersect(ithHashArr,x), origRecon2MetKEGGsHash,'UniformOutput',0);
    commonKEGGsArrHashMatchLen = cellfun(@(x) length(x(x~=0)),commonKEGGsArrHash);
    [maxCommon1 maxCommonIdx1] = max(commonKEGGsArrHashMatchLen);
    toSearchIdxs = find(commonKEGGsArrHashMatchLen==maxCommon1);
    commonKEGGsArr = cellfun(@(x) intersect(ithMetKEGGs,x), origRecon2MetKEGGs(toSearchIdxs),'UniformOutput',0);
    %commonKEGGsArr = cellfun(@(x) intersect(ithMetKEGGs,x), origRecon2MetKEGGs,'UniformOutput',0);
    [maxCommon maxCommonIdx] = max(cellfun(@(x) length(x(~strcmp(x,''))),commonKEGGsArr));
    maxCommonIdx = toSearchIdxs(maxCommonIdx);
    if maxCommon > 0
        disp(bigModelReconc.rxnNames{i})
        disp(origRecon2.rxnNames{maxCommonIdx})
        matchRxnNum = matchRxnNum+1;
        disp(matchRxnNum)
        sameMatchNum = sum(strcmp(searchRxnNames,origRecon2.rxnNames{maxCommonIdx}));
        bigModelReconc.rxnNames{i} = [origRecon2.rxnNames{maxCommonIdx} '_' num2str(sameMatchNum+1)];
        bigModelReconc.rxns{i} = [origRecon2.rxns{maxCommonIdx} '_' num2str(sameMatchNum+1)];
    end
    % for j=1:length(origRecon2MetKEGGs)
    %     commonKEGGs = intersect(ithMetKEGGs,origRecon2MetKEGGs{j});
    %     commonKEGGs = commonKEGGs(~strcmp('',commonKEGGs));
    %     if length(commonKEGGs) > maxCommon
    %         maxCommon = maxCommon+1;
    %         bigModelReconc.rxns{i} = origRecon2.rxns{j};
    %     end
    % end
    % if maxCommon ~= 0
    %     matchRxnNum = matchRxnNum+1;
    % end
end
end
    
if 0
bigModelReconc = bigModelAdded;
reconcSubArr = unique(bigModelReconc.subSystems);
reconcSubMatchArr = zeros(length(reconcSubArr),1);
matchRxnNum = 0;
for i=1:length(bigModelAdded.rxnECNums)
    isReconc = 0;
    for j=1:length(bigModelAdded.rxnECNums{i})
        ECNum = bigModelAdded.rxnECNums{i}{j};
        matchIdxs = find(strcmp(origRecon2.rxnECNumbers,ECNum));
        if length(matchIdxs) > 0
            matchIdx = 1;
            for k=1:length(matchIdxs)
                if isempty(regexp(origRecon2.rxns{matchIdxs(matchIdx)},'rxn'))
                    matchIdx = k;
                end
            end
            bigModelReconc.rxnNames{i} = origRecon2.rxnNames{matchIdxs(1)};
            bigModelReconc.rxns{i} = origRecon2.rxns{matchIdxs(1)};
            isReconc = 1;
            disp('reconcile')
            break;
        end
        disp([num2str(i) ' ' num2str(j) ' ' num2str(length(matchIdxs))])
    end
    if isReconc
        subIdx = find(strcmp(reconcSubArr,bigModelAdded.subSystems{i}));
        reconcSubMatchArr(subIdx) = reconcSubMatchArr(subIdx)+1;
        matchRxnNum = matchRxnNum+1;
    end
end

for i=1:length(reconcSubArr)
    reconcSubMatchArr(i) = reconcSubMatchArr(i)/sum(strcmp(bigModelAdded.subSystems,reconcSubArr{i}));
end

writeData({reconcSubArr,reconcSubMatchArr},'/home/fs01/yw595/reconcBigModelRxnStat.txt','\t',{'subsystem','fractionreconc'});

end

if 0

fullMatrix = full(bigModelReconc.S);
connectivities = cellfun(@(x) sum(fullMatrix(strcmp(bigModelReconc.mets,x),:)~=0),bigModelReconc.mets);
connectivitiesArr = {'1','2','3','4','5','6','7','8','9','10','>10'};
connectivitiesCountArr = zeros(length(connectivitiesArr),1);
matchMetNum = 0;
for i=1:length(bigModelReconc.metNames)
    keggID = bigModelReconc.mets{i};
    keggID = keggID(4:8); keggID = ['C' keggID];
    matchIdxs = find(strcmp(origRecon2.metKeggID,keggID));
    if length(matchIdxs) > 0
        i
        matchMetNum = matchMetNum+1;
        bigModelReconc.metNames{i} = origRecon2.metNames{matchIdxs(1)};
        bigModelReconc.mets{i} = origRecon2.mets{matchIdxs(1)};
        connectivity = connectivities(i);
        if connectivity > 10
            connectivity = 11;
        end
        connectivitiesCountArr(connectivity) = connectivitiesCountArr(connectivity)+1;
    end
end

for i=1:length(connectivitiesCountArr)
    if i==11
        connectivitiesCountArr(i) = connectivitiesCountArr(i)/sum(connectivities>10);
    else
        connectivitiesCountArr(i) = connectivitiesCountArr(i)/sum(connectivities==i);
    end
end

writeData({connectivitiesArr,connectivitiesCountArr},'/home/fs01/yw595/reconcBigModelMetStat.txt','\t',{'conn','fractionconn'});

end

%nonsense = nonsense+1;

if 1
mocatFI = fopen('/home/fs01/yw595/MATLAB/SEED/input/MGMData/mocatCounts.txt');
mocatData = textscan(mocatFI,'%s%s','HeaderLines',0);
mocatECs = mocatData{1}; mocatCounts = mocatData{2};
bigModelReconc.express = zeros(length(bigModelReconc.rxns),1);
bigModelReconc.metKeggID = cell(length(bigModelReconc.rxns),1);
for i=1:length(bigModelReconc.rxnECNums)
    express = 0;
    for j=1:length(bigModelReconc.rxnECNums{i})
        ECNum = bigModelReconc.rxnECNums{i}{j};
        matchIdxs = find(strcmp(mocatECs,ECNum));
        if ~isempty(matchIdxs)
            %mocatCounts{matchIdxs(1)}
            %str2num(mocatCounts{matchIdxs(1)})
            express = express + str2num(mocatCounts{matchIdxs(1)});
        end
        
        %express = express + sum(cellfun(@(x) str2num(x), mocatCounts));
    end
    bigModelReconc.express(i) = express;
    bigModelReconc.metKeggID{i} = '';
end
end

if 0
cpdMatchNum = 0;
for i=1:length(bigModelReconc.metKEGGs)
    i
    searchMetNames = cellfun(@(x) x(1:regexp(x,'_')), bigModelReconc.metNames,'UniformOutput',0);
    if length(bigModelReconc.metKEGGs{i})==6
        matchIdx = find(strcmp(origRecon2.metKeggID,bigModelReconc.metKEGGs{i}));
        if ~isempty(matchIdx)
            cpdMatchNum = cpdMatchNum+1;
            if length(matchIdx) > 1
                matchIdx = matchIdx(1);
            end
            sameMatchNum = sum(strcmp(searchMetNames,origRecon2.metNames{matchIdx}));
            bigModelReconc.mets{i} = [origRecon2.mets{matchIdx} '_' num2str(sameMatchNum+1)];
            bigModelReconc.metNames{i} = [origRecon2.metNames{matchIdx} '_' num2str(sameMatchNum+1)];
        end
    end
end
end

if 1
for i=1:length(bigModelReconc.rxns)
    if strcmp(bigModelReconc.rxns{i}(end-1:end),'_1')
        bigModelReconc.rxns{i} = bigModelReconc.rxns{i}(1:end-2);
        bigModelReconc.rxnNames{i} = bigModelReconc.rxnNames{i}(1:end-2);
    end
end
for i=1:length(bigModelReconc.mets)
    if strcmp(bigModelReconc.mets{i}(end-1:end),'_1')
        bigModelReconc.mets{i} = bigModelReconc.mets{i}(1:end-2);
        bigModelReconc.metNames{i} = bigModelReconc.metNames{i}(1:end-2);
    end
end

end

printModel(bigModelReconc,bigModelReconc.express,'/home/fs01/yw595/reducedFBABigModelReconc.txt','/home/fs01/yw595/frequentMetsBigModelReconc.txt');
save([outputDir filesep 'reconcBigModelRecon2' filesep 'reconcBigModelRecon2.mat'],'bigModelReconc');




