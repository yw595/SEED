configSEED;

if 0
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
    bigModelTableFlux = restrictPrevTE(bigModelTableFlux,convertArrSave(5,:,3));
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

methodsList = {'FBA','EFlux','FALCON','GXFBA'};
uniqSubs = unique(bigModelTableFlux.subSystems);
convertArr2 = std(convertArr,0,1);
for i=3:3%length(methodsList)
    if useScramble
        meanNorm = mean(convertArr(1,:,i),1);
        if useERP
            meanObese = mean(convertArrObese(1,:,i),1);
        end
        methodstds = std(convertArr(2:end,:,1),0,1)';
    else
        meanNorm = mean(convertArr(:,:,i),1);
        if useERP
            meanObese = mean(convertArrObese(:,:,i),1);
        end
        methodstds = std(convertArr(:,:,1),0,1)';
    end

    labels = {}; yvals = []; grouplabels = {}; labels1 = {}; yvals1 = [];
    labelsDiffERP = {}; yvalsDiffERP = []; yvals2DiffERP = []; yvals3DiffERP = [];
    for j=1:length(uniqSubs)
        labels1{end+1} = uniqSubs{j};
        yvals1(end+1) = sum(convertArr2(1,strcmp(bigModelTableFlux.subSystems,uniqSubs{i}),3)~=0);
        if useERP
            labelsDiffERP{end+1} = uniqSubs{j};
            yvalsDiffERP(end+1) = abs(sum(meanNorm(strcmp(bigModelTableFlux.subSystems,uniqSubs{j}))-meanObese(strcmp(bigModelTableFlux.subSystems,uniqSubs{j}))));
            yvals2DiffERP(end+1) = abs(sum(methodstds(strcmp(bigModelTableFlux.subSystems,uniqSubs{j}))));
            if useScramble
                yvals3DiffERP(end+1) = yvalsDiffERP(end)/yvals2DiffERP(end);
            end
        end        
    end
    [yvalsDiffERP, sortIdxs] = sort(yvalsDiffERP,'descend');
    labelsDiffERP = labelsDiffERP(sortIdxs);
    yvals2DiffERP = yvals2DiffERP(sortIdxs);
    yvals3DiffERP = yvals3DiffERP(sortIdxs);
    labelsTemp = labelsDiffERP(~isnan(yvals3DiffERP) & ~isinf(yvals3DiffERP));
    yvalsTemp = yvals3DiffERP(~isnan(yvals3DiffERP) & ~isinf(yvals3DiffERP));
    [yvalsTemp, sortIdxsTemp] = sort(yvalsTemp,'descend');
    labelsTemp = labelsTemp(sortIdxsTemp);
    [yvals1, sortIdxs] = sort(yvals1,'descend');
    labels1 = labels1(sortIdxs);
    labels1 = addIdxStrings(labels1);
    labelsTemp = addIdxStrings(labelsTemp);
    labelsDiffERP = addIdxStrings(labelsDiffERP);

    writeData({labelsDiffERP,yvalsDiffERP,yvalsDiffERP-yvals2DiffERP,yvalsDiffERP+yvals2DiffERP},['/home/fs01/yw595/' methodsList{i} 'OnlyDiffERP.txt'],'\t',{'sub','avgdiff','lower','upper'});
    writeData({labelsTemp,yvalsTemp},['/home/fs01/yw595/' methodsList{i} 'OnlyDiffERPAdjusted.txt'],'\t',{'sub','avgdiff'});
    writeData({labels1,yvals1},['/home/fs01/yw595/fourMethodsSubsStds' methodsList{i} 'Only.txt'],'\t',{'sub','std'});
    printModel(bigModelReconc,convertArr(1,:,i),['/home/fs01/yw595/reducedFBABigModelReconc' methodsList{i} 'Flux.txt'],['/home/fs01/yw595/frequentMetsBigModelReconc' methodsList{i} 'Flux.txt']);
    printModel(bigModelReconc,abs(convertArr(1,:,i)-convertArrObese(1,:,i)),['/home/fs01/yw595/reducedFBABigModelReconc' methodsList{i} 'DiffFlux.txt'],['/home/fs01/yw595/frequentMetsBigModelReconc' methodsList{i} 'DiffFlux.txt']);
    printModel(bigModelReconc,methodstds,['/home/fs01/yw595/reducedFBABigModelReconc' methodsList{i} 'Stds.txt'],['/home/fs01/yw595/frequentMetsBigModelReconc' methodsList{i} 'Stds.txt']);
end
printModel(bigModelReconc,bigModelReconc.express,'/home/fs01/yw595/reducedFBABigModelReconcExpress.txt','/home/fs01/yw595/frequentMetsBigModelReconcExpress.txt');

for i=1:length(uniqSubs)
    for j=1:length(methodsList)
        labels{end+1} = [uniqSubs{i} '_' methodsList{j}];
        yvals(end+1) = mean(convertArr2(1,strcmp(bigModelTableFlux.subSystems,uniqSubs{i}),j));
        grouplabels{end+1} = methodsList{j};
    end
end
writeData({labels,yvals,grouplabels},'/home/fs01/yw595/fourMethodsSubsStds.txt','\t',{'subsAndMethod','std','method'});

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
sortUniqSubs = addIdxStrings(sortUniqSubs);
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
