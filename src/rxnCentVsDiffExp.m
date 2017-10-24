configSEED;

outputDir1 = [outputDir filesep 'rxnCentVsDiffExp'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1);
end

%make binary rxnMatrix, reactions linked by mets, excluding union of topMetIdxs, most connected .001, and predefined excludeMets
bigModel = bigModelAdded;
connMatrixTable = makeConnMatrix(bigModel);
rxnMatrix = zeros(length(bigModel.rxns),length(bigModel.rxns));
[~, topMetIdxsOrig] = sort(sum(abs(bigModel.S)~=0,2),'descend');
topMetIdxs = topMetIdxsOrig(1:floor(length(topMetIdxsOrig)*.001));
excludeMets = {'ATP','ADP','Phosphate','NAD','NADH','CoA','PPi','NADP','NADPH','O2','AMP','NH3','CO2','UDP','S-Adenosyl-L-methionine','cpd00001','cpd00002','cpd00003','cpd00004','cpd00005','cpd00006','cpd00007','cpd00008','cpd00009','cpd00010','cpd00011','cpd00012','cpd00013','cpd00014','cpd00015','cpd00016','cpd00017','cpd00018','cpd00019','cpd00020','cpd00021'};
[~,intersectIdxs,~] = intersect(bigModel.metNames,excludeMets);
topMetIdxs = union(topMetIdxs,intersectIdxs);
for i=1:length(bigModel.rxns)
    metIdxs = find(bigModel.S(:,i)~=0);
    for j=1:length(metIdxs)
        if ~any(metIdxs(j)==topMetIdxs)
            connRxnIdxs = find(bigModel.S(metIdxs(j),:));
            if i<100
                disp(length(connRxnIdxs))
                disp(bigModel.metNames{metIdxs(j)})
            end
            for k=1:length(connRxnIdxs)
                rxnMatrix(i,connRxnIdxs(k))=1;
            end
        end
    end
end

centsRxnTable = betweenness_centrality(sparse(rxnMatrix));
ECsToCents = containers.Map;
for i=1:length(bigModel.rxnNames)
    rxnECs = bigModel.rxnECNums{i};
    for j=1:length(rxnECs)
        if isKey(ECsToCents,rxnECs{j})
            ECsToCents(rxnECs{j}) = ECsToCents(rxnECs{j}) + centsRxnTable(i);
        else
            ECsToCents(rxnECs{j}) = centsRxnTable(i);
        end
    end
end
ECsToCentsKeys = keys(ECsToCents);
correspondingCents = [];
for i=1:length(ECsToCentsKeys)
    correspondingCents(i) = mean(ECsToCents(ECsToCentsKeys{i}))*2/((length(ECsToCentsKeys)-1)*(length(ECsToCentsKeys)-2));
end
writeData({ECsToCentsKeys,correspondingCents},[baseDir filesep 'ECsToCents.txt'],'\t');
%centsTable = betweenness_centrality(sparse(connMatrixTable));
ask = ask+1;

titleString = 'Abs Log Diff Expression Vs. Centrality';
yvals = []; xvals = [];
count=0;
for i=1:length(bigModel.rxns)
    if isKey(rxnsToDiffExp,bigModel.rxns{i})
        count=count+1;
        if count<100
            disp(i)
            %disp(rxnsToDiffExp(bigModel.rxns{i}))
            disp(centsRxnTable(i))
            disp(sum(xvals~=0))
        end
        yvals(end+1) = rxnsToDiffExp(bigModel.rxns{i});
        xvals(end+1) = centsRxnTable(i);
    end
end
[~, sortIdxs] = sort(xvals); xvals=xvals(sortIdxs(1:end-1)); yvals=yvals(sortIdxs(1:end-1));
[rhoSpear pvalSpear] = corr(xvals',yvals','type','Spearman');
[rhoPear pvalPear] = corr(xvals',yvals','type','Pearson');
makeBar(xvals,yvals,titleString,outputDir1,'ylabelString', 'Abs Log Diff Expression','xlabelString','Centrality','isScatter',1);
writeForGGPlot(xvals,yvals,[outputDir filesep 'Abs_Log_Diff_Expression_Vs_Centrality.txt']);



