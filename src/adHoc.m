if 0
    %configSEED;
modelNamesToModelsValues = values(modelNamesToModels);
writeOutModel(bigModelTable,[outputDir filesep 'bigModelTable.txt']);

bigModel = bigModelTable;
bigModel = addMustEx(bigModel);
testModel=modelNamesToModelsValues{1};
%[biomassMets biomassCoeffs] = makeMergedBiomass({testModel},bigModel.mets);
bigModel = mergeModels(bigModel,testModel);
%sols = testBiomass(bigModel,biomassMets,biomassCoeffs,0);
bigModel.c = zeros(length(bigModel.rxns),1);
bigModel = changeObjective(bigModel,testModel.rxns{strcmp(testModel.rxnNames,'Biomass')});
sol = optimizeCbModel(bigModel);
end

if 0
compares1={normFluxesNormal(:,2),normFluxesNormal(:,3),normFluxesNormal(:,1),normFluxesNormal(:,1)};
compares2={normFluxesObese(:,2),normFluxesObese(:,3),normFluxesNormal(:,2),normFluxesNormal(:,3)};
titles = {'Normal and Obese EFluxes','Normal and Obese FALCON Fluxes','Normal FBA and EFluxes','Normal FBA and FALCON Fluxes'};
legendLabelsArray = {{'Normal','Obese'},{'Normal','Obese'},{'FBA','EFlux'},{'FBA','FALCON'}};
for i=1:length(compares1)
    [~, ~, ~, ~, ~, ~, ~, ~, ~, SSAbsDiffs1TopTen, SSAbsDiffs2TopTen,SSAbsDiffsTopTenNames]=analyzeDiffFluxes(bigModel,compares1{i},compares2{i},1);
    titleString = titles{i};
    xvals = 1:10; yvals = [SSAbsDiffs1TopTen SSAbsDiffs2TopTen];
    legendLabels = legendLabelsArray{i};
    makeBar(xvals,yvals,titleString,outputDir,'ylabelString','Flux','xlabelString','Subsystem','xlabels',SSAbsDiffsTopTenNames,'legendLabels',legendLabels);
    writeForGGPlot(xvals,yvals,[outputDir filesep strrep(titleString,' ','_') '.txt'],SSAbsDiffsTopTenNames,legendLabels);
end

xvals = normFluxesNormal(:,2); yvals = normFluxesNormal(:,3)/max(normFluxesNormal(:,3)); titleString = 'Normal FALCON Vs EFluxes';
makeBar(xvals,yvals,titleString,outputDir,'ylabelString','FALCON Flux','xlabelString','EFlux','isScatter',1);
writeForGGPlot(xvals,yvals,[outputDir filesep 'Normal_FALCON_Vs_EFluxes.txt']);
end

if 0
configSEED;
bigModelsAppend = bigModelsAccum; bigModelsAppend{end+1}=bigModelTable;
xvals = 1:length(bigModelsAppend)+5; yvals = zeros(length(bigModelsAppend)+5,4); titleString = 'Reactions Added Per Model'; xlabels = {};
legendLabels = {'numNormal','numExch','numTrans','numUnclass'};
for i=1:length(bigModelsAppend)
    numUnclass = sum(strcmp(bigModelsAppend{i}.subSystems,''));
    numTrans = sum(cellfun(@(x) ~isempty(regexpi(x,'transport')), bigModelsAppend{i}.subSystems));
    numExch = sum(cellfun(@(x) ~isempty(regexpi(x,'exchange|MUST_EX')), bigModelsAppend{i}.subSystems));
    numNormal = length(bigModelsAppend{i}.subSystems)-numUnclass-numTrans-numExch;
    yvals(i,1) = numNormal; yvals(i,2) = numExch; yvals(i,3) = numTrans; yvals(i,4) = numUnclass;
    if i==length(bigModelsAppend)
        xlabels{i} = 'bigModelTable';
    elseif i==length(bigModelsAppend)-1
        xlabels{i} = 'bigModelAccum';
    else
        xlabels{i} = '';
    end
end
for i=1:5
    xlabels{end+1}='';
end
makeBar(xvals,yvals,titleString,outputDir,'ylabelString','Num Reactions','xlabelString','Model','xlabels',xlabels,'isStackBar',1,'legendLabels',legendLabels);
xlabels1 = keys(modelNamesToModels);
xlabels1{end+1} = 'Added All Model SEED Rxns';
for i=1:length(xlabels1)
    paddedIdx = num2str(i);
    while length(paddedIdx)<3
        paddedIdx = ['0' paddedIdx];
    end
    xlabels1{i} = [paddedIdx '_' xlabels1{i}];
end
writeForGGPlot(xvals(1:end-5),yvals(1:end-5,:),[outputDir filesep 'Reactions_Added_Per_Model.txt'],xlabels1,legendLabels);
end

if 1
%writeOutModel(bigModelTableTest,[outputDir filesep 'bigModelTableTest.txt']);

%bigModel = bigModelTable;

%connMatrix = makeConnMatrix(bigModel);
%cents = betweenness_centrality(sparse(connMatrix));

%bigModel = addMustEx(bigModel);

[biomassMets biomassCoeffs] = makeMergedBiomass(values(modelNamesToModels),bigModel.mets);

modelNamesToModelsValues = values(modelNamesToModels);
%bigModelAdded = mergeSmallest(bigModel,modelNamesToModelsValues);

sols = testBiomass(bigModelAdded,biomassMets,biomassCoeffs,0);

xvals = 1:length(sols); titleString = 'biomassMetProd'; yvals = [];
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
writeData({xvals,yvals,xlabels},[transferDir filesep 'biomassMetProd.txt'],'\t',{'xvals','yvals','xlabels'});
%makeBar(xvals,yvals,titleString,outputDir,'ylabelString','Flux','xlabelString','Metabolite','xlabels',xlabels);

% xvals = 1:length(bigModels); titleString = 'numRxnsAdded'; yvals = cellfun(@(x) length(x.rxns), bigModels);
% xlabels = {};
% for i=1:length(bigModels)
%     xlabels{i} = num2str(i);
% end
% makeBar(xvals,yvals,titleString,outputDir,'ylabelString', 'numRxns','xlabelString','Model Number','xlabels',xlabels);

% rxnsToExpress = mapExpToRxns(ECsToRxns,[baseDir filesep 'testEC3.txt']);
% titleString = 'Total Expression Vs. Centrality';
% yvals = []; xvals = [];
% for i=1:length(bigModel.rxns)
%     if isKey(rxnsToExpress,bigModel.rxns{i})
%         yvals(end+1) = rxnsToExpress(bigModel.rxns{i});
%         xvals(end+1) = cents(i);
%     end
% end
% makeBar(xvals,yvals,titleString,outputDir,'ylabelString', 'Total Expression','xlabelString','Centrality','isScatter',1);
end
