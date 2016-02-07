%configSEED;
outputDir1 = [outputDir filesep 'rxnCentVsDiffExp'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1);
end
connMatrixTable = makeConnMatrix(bigModel);
rxnMatrix = zeros(length(bigModel.rxns),length(bigModel.rxns));
[~, topMetIdxs] = sort(sum(bigModelTable.S~=0,2),'descend');
topMetIdxs = topMetIdxs(1:floor(length(topMetIdxs)*.001));
for i=1:length(bigModel.rxns)
    metIdxs = find(bigModel.S(:,i)~=0);
    for j=1:length(metIdxs)
        if ~any(metIdxs(j)==topMetIdxs)
            connRxnIdxs = find(bigModel.S(metIdxs(j),:));
            for k=1:length(connRxnIdxs)
                rxnMatrix(i,k)=1;
            end
        end
    end
end
centsRxnTable = betweenness_centrality(sparse(rxnMatrix));
centsTable = betweenness_centrality(sparse(connMatrixTable));

titleString = 'Abs Log Diff Expression Vs. Centrality';
yvals = []; xvals = [];
for i=1:length(bigModel.rxns)
    if isKey(rxnsToDiffExp,bigModel.rxns{i})
        yvals(end+1) = rxnsToDiffExp(bigModel.rxns{i});
        xvals(end+1) = centsRxnTable(i);
    end
end
[~, sortIdxs] = sort(xvals); xvals=xvals(sortIdxs(1:end-100)); yvals=yvals(sortIdxs(1:end-100));
[rhoSpear pvalSpear] = corr(xvals',yvals','type','Spearman');
[rhoPear pvalPear] = corr(xvals',yvals','type','Pearson');
makeBar(xvals,yvals,titleString,outputDir1,'ylabelString', 'Abs Log Diff Expression','xlabelString','Centrality','isScatter',1);