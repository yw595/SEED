configSEED;

writeOutModel(bigModelTableTest,[baseDir filesep 'bigModelTableTest.txt']);

bigModel = bigModelTable;

connMatrix = makeConnMatrix(bigModel);
cents = betweenness_centrality(sparse(connMatrix));

bigModel = addMustEx(bigModel);

[biomassMets biomassCoeffs] = makeMergedBiomass(values(modelNamesToModels),bigModel.mets);

modelNamesToModelsValues = values(modelNamesToModels);
bigModelAdded = mergeSmallest(bigModel,modelNamesToModelsValues);

sols = testBiomass(bigModelAdded,biomassMets,biomassCoeffs,0);

xvals = 1:length(sols); titleString = 'biomassMetProd'; yvals = []; outputDir= baseDir;
for i=1:length(sols)
    yvals(i) = sols{i}.f;
end
disp(sum(yvals>0)/length(yvals))
xlabels = {};
for i=1:length(biomassMets)
    xlabels{i} = bigModel.metNames{strcmp(bigModel.mets,biomassMets{i})};
    temp = xlabels{i};
    xlabels{i} = temp(1:min(10,length(temp)));
end
makeBar(xvals,yvals,titleString,outputDir,'ylabelString','Flux','xlabelString','Metabolite','xlabels',xlabels);

xvals = 1:length(bigModels); titleString = 'numRxnsAdded'; yvals = cellfun(@(x) length(x.rxns), bigModels); outputDir= baseDir;
xlabels = {};
for i=1:length(bigModels)
    xlabels{i} = num2str(i);
end
makeBar(xvals,yvals,titleString,outputDir,'ylabelString', 'numRxns','xlabelString','Model Number','xlabels',xlabels);

rxnsToExpress = mapExpToRxns(ECsToRxns,[baseDir filesep 'testEC3.txt']);
titleString = 'Total Expression Vs. Centrality';
yvals = []; xvals = [];
for i=1:length(bigModel.rxns)
    if isKey(rxnsToExpress,bigModel.rxns{i})
        yvals(end+1) = rxnsToExpress(bigModel.rxns{i});
        xvals(end+1) = cents(i);
end
makeBar(xvals,yvals,titleString,outputDir,'ylabelString', 'Total Expression','xlabelString','Centrality','isScatter',1);