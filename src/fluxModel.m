ERPData = textscan(fopen('/home/fs01/yw595/MATLAB/SEED/input/MGMData/ERPBowtieComplete.txt'),'%s%s%s%s','Delimiter','\t','HeaderLines',0);
xenoData = textscan(fopen('/home/fs01/yw595/MATLAB/SEED/input/xenobiotics/GSM935962_A1_EtOH_CDS.EC.txt'),'%s%s%s','Delimiter','\t','HeaderLines',0);
metabolomeData = textscan(fopen('/home/fs01/yw595/MATLAB/SEED/input/fecalmicrobiomeauto.txt'),'%s%s','Delimiter','\t','HeaderLines',0);
hadzaData = textscan(fopen('/home/fs01/yw595/output3.EC.txt'),'%s%s','Delimiter','\t','HeaderLines',0);

if 1
configSEED;
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
bigModelTableFlux = bigModelTableFluxSave;
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
if 1
fecalMet = metabolomeData{2};
fecalMet = fecalMet(~strcmp(fecalMet,''));
fecalMet = cellfun(@(x) ['C' x], fecalMet, 'UniformOutput', 0);
bigModelTableFlux = bigModelTableFluxSave;
bigModelTableFluxRestrict = bigModelTableFlux;
for i=1:length(bigModelTableFluxRestrict.rxns)
    if strcmp(bigModelTableFluxRestrict.subSystems{i},'Exchange')
        metIDs = bigModelTableFluxRestrict.metKEGGs(bigModelTableFluxRestrict.S(:,i)~=0);
        %metIDs = cellfun(@(x) ['C' x(4:8)], metIDs, 'UniformOutput',0);
        if length(intersect(metIDs,fecalMet)) ~= length(metIDs)
            disp(i)
            bigModelTableFluxRestrict.lb(i) = 0;
            bigModelTableFluxRestrict.ub(i) = 0;
        end
    end
end
bigModelTableFlux = bigModelTableFluxRestrict;
end
if 1
normFluxesNormalArr = {};
useERP = 0;
useXeno = 0;
useHadza = 1;
for z1=1:2
zLim = 2;
if useXeno==1 or useHadza==1
    zLim = 1;
end
for z=1:zLim

    expressionIDs = {};
    expressionData = [];
    expressionSDs = [];
    
    if useERP
        for i=1:length(ERPData{1})
            i
            ECNums = ERPData{1}{i};
                for k =1:length(bigModelTableFlux.rxnECNums)
                    bigModelECNums = bigModelTableFlux.rxnECNums{k};
                    if any(strcmp(ECNums,bigModelECNums))
                        if ~any(strcmp(expressionIDs,bigModelTableFlux.rxns{k}))
                            expressionIDs{end+1} = bigModelTableFlux.rxns{k};
                            if z==1
                                expressionData(end+1) = str2num(ERPData{2}{i});
                            else
                                expressionData(end+1) = str2num(ERPData{3}{i});
                            end
                            expressionSDs(end+1) = 1;
                        end
                    end
                end
        end
    elseif useXeno
        for i=1:length(xenoData{1})
            i
            ECNums = xenoData{1}{i};
                for k =1:length(bigModelTableFlux.rxnECNums)
                    bigModelECNums = bigModelTableFlux.rxnECNums{k};
                    if any(strcmp(ECNums,bigModelECNums))
                        if ~any(strcmp(expressionIDs,bigModelTableFlux.rxns{k}))
                            expressionIDs{end+1} = bigModelTableFlux.rxns{k};
                            expressionData(end+1) = str2num(xenoData{2}{i});
                            expressionSDs(end+1) = str2num(xenoData{3}{i});
                        end
                    end
                end
        end
    elseif useHadza
        for i=1:length(hadzaData{1})
            i
            ECNums = hadzaData{1}{i};
                for k =1:length(bigModelTableFlux.rxnECNums)
                    bigModelECNums = bigModelTableFlux.rxnECNums{k};
                    if any(strcmp(ECNums,bigModelECNums))
                        if ~any(strcmp(expressionIDs,bigModelTableFlux.rxns{k}))
                            expressionIDs{end+1} = bigModelTableFlux.rxns{k};
                            expressionData(end+1) = str2num(hadzaData{2}{i});
                            expressionSDs(end+1) = 1;%str2num(xenoData{3}{i});
                        end
                    end
                end
        end
    else
        if z==1
            expressionIDs = rxnsToExpressNormKeys;
        else
            expressionIDs = rxnsToExpressObeseKeys;
        end
        expressionData = [];
        for i=1:length(expressionIDs)
            if z==1
                expressionData(i) = rxnsToExpressNorm(expressionIDs{i});
            else
                expressionData(i) = rxnsToExpressObese(expressionIDs{i});
            end
        end
        expressionSDs = ones(size(expressionData));
    end

    normFluxes = runAllFluxMethods(bigModelTableFlux,expressionData,expressionSDs,expressionIDs);
    
    if z==1
        normFluxesNormal = normFluxes;
        normFluxesNormalArr{z1} = normFluxesNormal;
    else
        normFluxesObese = normFluxes;
    end
end
end
end

if 1
convertArr = [];
for i=1:length(normFluxesNormalArr)
    convertArr(i,:,:) = normFluxesNormalArr{i};
end
FBAStds = std(convertArr(:,:,1),0,1)';
EFluxStds = std(convertArr(:,:,2),0,1)';
FALCONStds = std(convertArr(:,:,3),0,1)';
GXFBAStds = std(convertArr(:,:,4),0,1)';
uniqSubs = unique(bigModelTableFlux.subSystems);
methodsList = {'FBA','EFlux','FALCON','GXFBA'};
labels = {}; yvals = []; grouplabels = {}; labels1 = {}; yvals1 = [];
convertArr2 = std(convertArr,0,1);
for i=1:length(uniqSubs)
    labels1{end+1} = uniqSubs{i};
    yvals1(end+1) = sum(convertArr2(1,strcmp(bigModelTableFlux.subSystems,uniqSubs{i}),3)~=0);
    %yvals1(end+1) = mean(convertArr2(1,strcmp(bigModelTableFlux.subSystems,uniqSubs{i}),3));
    for j=1:length(methodsList)
        labels{end+1} = [uniqSubs{i} '_' methodsList{j}];
        yvals(end+1) = mean(convertArr2(1,strcmp(bigModelTableFlux.subSystems,uniqSubs{i}),j));
        grouplabels{end+1} = methodsList{j};
    end
end
[yvals1, sortIdxs] = sort(yvals1,'descend');
labels1 = labels1(sortIdxs);
for i=1:length(labels1)
    ithNum = num2str(i);
    while length(ithNum) < 3
        ithNum = ['0' ithNum];
    end
    labels1{i} = [ithNum '_' labels1{i}];
end

writeData({labels,yvals,grouplabels},'/home/fs01/yw595/fourMethodsSubsStds.txt','\t',{'subsAndMethod','std','method'});
writeData({labels1,yvals1},'/home/fs01/yw595/fourMethodsSubsStdsFALCONOnly.txt','\t',{'sub','std'});
bigModelReconc.FALCONStds = FALCONStds;
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
