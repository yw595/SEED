configSEED;

if 1 
    %configSEED;
rxnsToExpressObeseKeys = keys(rxnsToExpressObese);

% bigModelTableFlux = addMustEx(bigModelTable);
% [bigModelTableFlux, ~, addedModel] = mergeSmallest(bigModelTableFlux,modelNamesToModelsValues,1);
% [biomassMets biomassCoeffs] = makeMergedBiomass({addedModel},bigModelTableFlux.mets);
% bigModelTableFlux = addReaction(bigModelTableFlux,'BIOMASS', {bigModelTableFlux.mets{end}},[-1],1);
% bigModelTableFlux.c = zeros(length(bigModelTableFlux.rxns),1);
% bigModelTableFlux = changeObjective(bigModelTableFlux,'BIOMASS');

bigModelTableFlux = bigModelAdded;
bigModelTableFlux.c = zeros(length(bigModelTableFlux.rxns),1);
for i=1:length(bigModelTableFlux.rxns)
    bigModelTableFlux.rxns{i} = strrep(strrep(bigModelTableFlux.rxns{i},')',''),'(','');
end
bigModelTableFlux.genes = bigModelTableFlux.rxns;
bigModelTableFlux.grRules = bigModelTableFlux.rxns;
bigModelTableFlux.rules = bigModelTableFlux.rxns;
bigModelTableFlux.rxnGeneMat = eye(length(bigModelTableFlux.rxns));
bigModelTableFlux = changeObjective(bigModelTableFlux,testModel.rxns{strcmp(testModel.rxnNames,'Biomass')});
bigModelTableFlux.description = 'bigModel';
bigModelTableFlux.metKEGGs = bigModelReconc.metKEGGs;
end

if 0
    %bigModelTableFlux = bigModelTableFluxSave;
bigModelTableFluxRestrict = bigModelTableFlux;
for i=1:length(bigModelTableFluxRestrict.rxns)
    if strcmp(bigModelTableFluxRestrict.subSystems{i},'Transport') || strcmp(bigModelTableFluxRestrict.subSystems{i},'Exchange')
        bigModelTableFluxRestrict.lb(i) = convertArrSave(5,i,3);
        bigModelTableFluxRestrict.ub(i) = convertArrSave(5,i,3);
    end
end
bigModelTableFlux = bigModelTableFluxRestrict;
end

bigModelTableFluxSave = bigModelTableFlux;
if 0
    bigModelTableFlux = restrictFecalMet(bigModelTableFlux,metabolomeData);
end

if 0
    
normFluxesNormalArr = {};
normFluxesObeseArr = {};
useERP = 1;
useXeno = 0;
useHadza = 0;
useScramble = 1;
zLim = 2;
if useXeno==1 or useHadza==1
    zLim = 1;
end
readERPHadzaXeno;

for z1=1:10
    for z=1:zLim
        expressionData = expressionDataArr{z};
        expressionSDs = expressionSDsArr{z};
        if useScramble
            zScramble = setdiff(1:zLim,z);
            expressionDataAlt = expressionDataArr{zScramble};
            expressionSDsAlt = expressionSDsArr{zScramble};
            randIdxs = randperm(length(expressionData),floor(length(expressionData)/10));
            expressionDataScramble = expressionData;
            expressionDataScramble(randIdxs) = expressionDataAlt(randIdxs);
            expressionSDsScramble = expressionSDs;
            expressionSDsScramble(randIdxs) = expressionSDsAlt(randIdxs);
        end

        if useScramble && z1~=1
            disp('CONTRADICT')
            normFluxes = runAllFluxMethods(bigModelTableFlux,expressionDataScramble,expressionSDsScramble,expressionIDs);
        else
            normFluxes = runAllFluxMethods(bigModelTableFlux,expressionData,expressionSDs,expressionIDs);
        end
            
        if z==1
            normFluxesNormal = normFluxes;
            normFluxesNormalArr{z1} = normFluxesNormal;
        else
            normFluxesObese = normFluxes;
            normFluxesObeseArr{z1} = normFluxesObese;
        end
    end
end

end

if 1
convertArr = [];
for i=1:length(normFluxesNormalArr)
    convertArr(i,:,:) = normFluxesNormalArr{i};
end
if useERP
convertArrObese = [];
for i=1:length(normFluxesObeseArr)
    convertArrObese(i,:,:) = normFluxesObeseArr{i};
end
end

if useScramble
FBAMeanNorm = mean(convertArr(1,:,1),1);
EFluxMeanNorm = mean(convertArr(1,:,2),1);
FALCONMeanNorm = mean(convertArr(1,:,3),1);
GXFBAMeanNorm = mean(convertArr(1,:,4),1);
if useERP
FBAMeanObese = mean(convertArrObese(1,:,1),1);
EFluxMeanObese = mean(convertArrObese(1,:,2),1);
FALCONMeanObese = mean(convertArrObese(1,:,3),1);
GXFBAMeanObese = mean(convertArrObese(1,:,4),1);
end
else
FBAMeanNorm = mean(convertArr(:,:,1),1);
EFluxMeanNorm = mean(convertArr(:,:,2),1);
FALCONMeanNorm = mean(convertArr(:,:,3),1);
GXFBAMeanNorm = mean(convertArr(:,:,4),1);
if useERP
FBAMeanObese = mean(convertArrObese(:,:,1),1);
EFluxMeanObese = mean(convertArrObese(:,:,2),1);
FALCONMeanObese = mean(convertArrObese(:,:,3),1);
GXFBAMeanObese = mean(convertArrObese(:,:,4),1);
end
end

if useScramble
FBAStds = std(convertArr(2:end,:,1),0,1)';
EFluxStds = std(convertArr(2:end,:,2),0,1)';
FALCONStds = std(convertArr(2:end,:,3),0,1)';
GXFBAStds = std(convertArr(2:end,:,4),0,1)';
else
FBAStds = std(convertArr(:,:,1),0,1)';
EFluxStds = std(convertArr(:,:,2),0,1)';
FALCONStds = std(convertArr(:,:,3),0,1)';
GXFBAStds = std(convertArr(:,:,4),0,1)';
end

uniqSubs = unique(bigModelTableFlux.subSystems);
methodsList = {'FBA','EFlux','FALCON','GXFBA'};
labels = {}; yvals = []; grouplabels = {}; labels1 = {}; yvals1 = [];
labelsDiffERP1 = {}; yvalsDiffERP1 = []; yvals2DiffERP1 = []; yvals3DiffERP1 = [];
labelsDiffERP2 = {}; yvalsDiffERP2 = []; yvals2DiffERP2 = []; yvals3DiffERP2 = [];
labelsDiffERP3 = {}; yvalsDiffERP3 = []; yvals2DiffERP3 = []; yvals3DiffERP3 = [];
labelsDiffERP4 = {}; yvalsDiffERP4 = []; yvals2DiffERP4 = []; yvals3DiffERP4 = [];
convertArr2 = std(convertArr,0,1);
for i=1:length(uniqSubs)
    labels1{end+1} = uniqSubs{i};
    yvals1(end+1) = sum(convertArr2(1,strcmp(bigModelTableFlux.subSystems,uniqSubs{i}),3)~=0);
    if useERP
    labelsDiffERP1{end+1} = uniqSubs{i};
    yvalsDiffERP1(end+1) = abs(sum(FBAMeanNorm(strcmp(bigModelTableFlux.subSystems,uniqSubs{i}))-FBAMeanObese(strcmp(bigModelTableFlux.subSystems,uniqSubs{i}))));
    yvals2DiffERP1(end+1) = abs(sum(FBAStds(strcmp(bigModelTableFlux.subSystems,uniqSubs{i}))));
    labelsDiffERP2{end+1} = uniqSubs{i};
    yvalsDiffERP2(end+1) = abs(sum(EFluxMeanNorm(strcmp(bigModelTableFlux.subSystems,uniqSubs{i}))-EFluxMeanObese(strcmp(bigModelTableFlux.subSystems,uniqSubs{i}))));
    yvals2DiffERP2(end+1) = abs(sum(EFluxStds(strcmp(bigModelTableFlux.subSystems,uniqSubs{i}))));
    labelsDiffERP3{end+1} = uniqSubs{i};
    yvalsDiffERP3(end+1) = abs(sum(FALCONMeanNorm(strcmp(bigModelTableFlux.subSystems,uniqSubs{i}))-FALCONMeanObese(strcmp(bigModelTableFlux.subSystems,uniqSubs{i}))));
    yvals2DiffERP3(end+1) = abs(sum(FALCONStds(strcmp(bigModelTableFlux.subSystems,uniqSubs{i}))));
    labelsDiffERP4{end+1} = uniqSubs{i};
    yvalsDiffERP4(end+1) = abs(sum(GXFBAMeanNorm(strcmp(bigModelTableFlux.subSystems,uniqSubs{i}))-GXFBAMeanObese(strcmp(bigModelTableFlux.subSystems,uniqSubs{i}))));
    yvals2DiffERP4(end+1) = abs(sum(GXFBAStds(strcmp(bigModelTableFlux.subSystems,uniqSubs{i}))));
    if useScramble
        yvals3DiffERP1(end+1) = yvalsDiffERP1(end)/yvals2DiffERP1(end);
        yvals3DiffERP2(end+1) = yvalsDiffERP2(end)/yvals2DiffERP2(end);
        yvals3DiffERP3(end+1) = yvalsDiffERP3(end)/yvals2DiffERP3(end);
        yvals3DiffERP4(end+1) = yvalsDiffERP4(end)/yvals2DiffERP4(end);
    end
    end
    %yvals1(end+1) = mean(convertArr2(1,strcmp(bigModelTableFlux.subSystems,uniqSubs{i}),3));
    for j=1:length(methodsList)
        labels{end+1} = [uniqSubs{i} '_' methodsList{j}];
        yvals(end+1) = mean(convertArr2(1,strcmp(bigModelTableFlux.subSystems,uniqSubs{i}),j));
        grouplabels{end+1} = methodsList{j};
    end
end
[yvalsDiffERP1, sortIdxs] = sort(yvalsDiffERP1,'descend');
labelsDiffERP1 = labelsDiffERP1(sortIdxs);
yvals2DiffERP1 = yvals2DiffERP1(sortIdxs);
yvals3DiffERP1 = yvals3DiffERP1(sortIdxs);
[yvalsDiffERP2, sortIdxs] = sort(yvalsDiffERP2,'descend');
labelsDiffERP2 = labelsDiffERP2(sortIdxs);
yvals2DiffERP2 = yvals2DiffERP2(sortIdxs);
yvals3DiffERP2 = yvals3DiffERP2(sortIdxs);
labelsTemp = labelsDiffERP2(~isnan(yvals3DiffERP2) & ~isinf(yvals3DiffERP2));
yvalsTemp = yvals3DiffERP2(~isnan(yvals3DiffERP2) & ~isinf(yvals3DiffERP2));
[yvalsTemp, sortIdxsTemp] = sort(yvalsTemp,'descend');
labelsTemp = labelsTemp(sortIdxsTemp);
[yvalsDiffERP3, sortIdxs] = sort(yvalsDiffERP3,'descend');
labelsDiffERP3 = labelsDiffERP3(sortIdxs);
yvals2DiffERP3 = yvals2DiffERP3(sortIdxs);
yvals3DiffERP3 = yvals3DiffERP3(sortIdxs);
[yvalsDiffERP4, sortIdxs] = sort(yvalsDiffERP4,'descend');
labelsDiffERP4 = labelsDiffERP4(sortIdxs);
yvals2DiffERP4 = yvals2DiffERP4(sortIdxs);
yvals3DiffERP4 = yvals3DiffERP4(sortIdxs);
[yvals1, sortIdxs] = sort(yvals1,'descend');
labels1 = labels1(sortIdxs);
for i=1:length(labels1)
    ithNum = num2str(i);
    while length(ithNum) < 3
        ithNum = ['0' ithNum];
    end
    labels1{i} = [ithNum '_' labels1{i}];
end
for i=1:length(labelsTemp)
    ithNum = num2str(i);
    while length(ithNum) < 3
        ithNum = ['0' ithNum];
    end
    labelsTemp{i} = [ithNum '_' labelsTemp{i}];
end
for i=1:length(labelsDiffERP1)
    ithNum = num2str(i);
    while length(ithNum) < 3
        ithNum = ['0' ithNum];
    end
    labelsDiffERP1{i} = [ithNum '_' labelsDiffERP1{i}];
end
for i=1:length(labelsDiffERP2)
    ithNum = num2str(i);
    while length(ithNum) < 3
        ithNum = ['0' ithNum];
    end
    labelsDiffERP2{i} = [ithNum '_' labelsDiffERP2{i}];
end
for i=1:length(labelsDiffERP3)
    ithNum = num2str(i);
    while length(ithNum) < 3
        ithNum = ['0' ithNum];
    end
    labelsDiffERP3{i} = [ithNum '_' labelsDiffERP3{i}];
end
for i=1:length(labelsDiffERP4)
    ithNum = num2str(i);
    while length(ithNum) < 3
        ithNum = ['0' ithNum];
    end
    labelsDiffERP4{i} = [ithNum '_' labelsDiffERP4{i}];
end

writeData({labelsTemp,yvalsTemp},'/home/fs01/yw595/FALCONOnlyDiffERPAdjusted.txt','\t',{'sub','avgdiff'});
writeData({labels,yvals,grouplabels},'/home/fs01/yw595/fourMethodsSubsStds.txt','\t',{'subsAndMethod','std','method'});
writeData({labels1,yvals1},'/home/fs01/yw595/fourMethodsSubsStdsFALCONOnly.txt','\t',{'sub','std'});
writeData({labelsDiffERP1,yvalsDiffERP1,yvalsDiffERP1-yvals2DiffERP1,yvalsDiffERP1+yvals2DiffERP1},'/home/fs01/yw595/FBAOnlyDiffERP.txt','\t',{'sub','avgdiff','lower','upper'});
writeData({labelsDiffERP2,yvalsDiffERP2,yvalsDiffERP2-yvals2DiffERP2,yvalsDiffERP2+yvals2DiffERP2},'/home/fs01/yw595/EFluxOnlyDiffERP.txt','\t',{'sub','avgdiff','lower','upper'});
writeData({labelsDiffERP3,yvalsDiffERP3,yvalsDiffERP3-yvals2DiffERP3,yvalsDiffERP3+yvals2DiffERP3},'/home/fs01/yw595/FALCONOnlyDiffERP.txt','\t',{'sub','avgdiff','lower','upper'});
writeData({labelsDiffERP4,yvalsDiffERP4,yvalsDiffERP4-yvals2DiffERP4,yvalsDiffERP4+yvals2DiffERP4},'/home/fs01/yw595/GXFBAOnlyDiffERP.txt','\t',{'sub','avgdiff','lower','upper'});
bigModelReconc.FALCONStds = FALCONStds;
printModel(bigModelReconc,convertArr(1,:,3),'/home/fs01/yw595/reducedFBABigModelReconcFALCONFlux.txt','/home/fs01/yw595/frequentMetsBigModelReconcFALCONFlux.txt');
printModel(bigModelReconc,bigModelReconc.FALCONStds,'/home/fs01/yw595/reducedFBABigModelReconcFALCONStds.txt','/home/fs01/yw595/frequentMetsBigModelReconcFALCONStds.txt');
printModel(bigModelReconc,bigModelReconc.express,'/home/fs01/yw595/reducedFBABigModelReconcExpress.txt','/home/fs01/yw595/frequentMetsBigModelReconcExpress.txt');
end

if 0
uniqSubs = unique(bigModelTableFlux.subSystems);
subPcts = [];
matchSignIdxs = find(sign(normFluxes(:,1)) == sign(normFluxes(:,2)) & sign(normFluxes(:,1)) == sign(normFluxes(:,3)));
for i=1:length(uniqSubs)
    subPcts(i) = sum(strcmp(bigModelTableFlux.subSystems(matchSignIdxs),uniqSubs{i}))/sum(strcmp(bigModelTableFlux.subSystems,uniqSubs{i}));
end
[sortSubPcts sortIdxs] = sort(subPcts,'descend');
sortUniqSubs = uniqSubs(sortIdxs);
for i=1:length(sortUniqSubs)
    ithIdx = num2str(i);
    while length(ithIdx) < 3
        ithIdx = ['0' ithIdx];
    end
    sortUniqSubs{i} = [ithIdx '_' sortUniqSubs{i}];
end
writeData({sortUniqSubs,sortSubPcts},'/home/fs01/yw595/pctMatchSubsystem.txt','\t',{'uniqSubs','pcts'});

matrixLabels = {'GX-FBA','E-Flux','FALCON','GX-FBA'};
corrMatrixLabelsX = {};
corrMatrixLabelsY = {};
corrVals = [];
for i=1:3
    for j=1:3
        corrMatrixLabelsX{end+1} = [num2str(i) '_' matrixLabels{i}];
        corrMatrixLabelsY{end+1} = [num2str(j) '_' matrixLabels{j}];
        if i > j
            corrVals(end+1) = corr(normFluxes(:,i),normFluxes(:,j),'type','Spearman');
        else
            corrVals(end+1) = 0;
        end
    end
end
writeData({corrMatrixLabelsX,corrMatrixLabelsY,corrVals},'/home/fs01/yw595/corrMatrixGeneExp.txt','\t',{'methodX','methodY','corr'});
end

if 0
uniqSubs = unique(bigModelTableFlux.subSystems);
uniqSubsDiffs = zeros(length(uniqSubs),3);

for i=1:4
    normFluxes1 = normFluxesNormal(:,i);
    normFluxes2 = normFluxesObese(:,i);
    for j=1:length(uniqSubs)
        uniqSubsDiffs(j,i) = sum(abs(normFluxes1(strcmp(bigModelTableFlux.subSystems,uniqSubs{j})) - normFluxes2(strcmp(bigModelTableFlux.subSystems,uniqSubs{j}))));
    end
end

[~,methodSortIdxs] = sort(uniqSubsDiffs(:,2),'descend');
uniqSubsNum = uniqSubs(methodSortIdxs);
for i=1:length(uniqSubsNum)
    ithLabel = num2str(i);
    while length(ithLabel) < 3
        ithLabel = ['0' ithLabel];
    end
    ithLabel = [ithLabel '_' uniqSubsNum{i}];
    uniqSubsNum{i} = ithLabel;
end
end

if 0
falconFluxesAbs = abs(normFluxesNormal(:,1)-normFluxesObese(:,2));
uniqueSubs = unique(bigModelAddedTemp.subSystems);
subsFluxDiffs = [];
for i=1:length(uniqueSubs)
    subsFluxDiffs(i) = sum(falconFluxesAbs(strcmp(bigModelTableFlux.subSystems,uniqueSubs{i})))/sum(strcmp(bigModelTableFlux.subSystems,uniqueSubs{i}));
end

subsFluxCorrArr = [];
subsConnCorrArr = [];
for i=1:length(uniqueSubs)
    for j=1:length(uniqueSubs)
        subsFluxCorrArr(end+1) = subsFluxDiffs(i)-subsFluxDiffs(j);
        subsConnCorrArr(end+1) = subsystemConnMatrixNorm(j,i);
    end
end

writeData({subsFluxCorrArr,subsConnCorrArr},'/home/fs01/yw595/subsConnFlux.txt','\t',{'subsFlux','subsConn'});
end
