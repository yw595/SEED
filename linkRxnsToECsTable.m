if 0
configSEED;
FI = fopen([baseDir filesep 'reactionsCleaned.csv']);
dataFields = textscan(FI,repmat('%s',1,9),'Delimiter',',');
rxnData = [dataFields{:}]; fclose(FI);
FI = fopen([baseDir filesep 'compoundsCleaned.csv']);
dataFields = textscan(FI,repmat('%s',1,10),'Delimiter',',');
cpdData = [dataFields{:}]; fclose(FI);
rxnsToECsTable = containers.Map; ECsToRxnsTable = containers.Map;
for i=1:length(rxnData)
    rxnName = rxnData{i,1};
    ECCell = strrep(rxnData{i,3},'|',':');
    [regex1 regex2] = regexp(ECCell,'(\d|-)+\.(\d|-)+\.(\d|-)+\.(\d|-)+');
    if ~isempty(regex1)
        ECNums = arrayfun(@(x,y) ECCell(x:y), regex1,regex2, 'UniformOutput',0);
        [rxnsToECsTable ECsToRxnsTable] = updateTwoMaps(rxnsToECsTable, ECsToRxnsTable, rxnName, ECNums);
    end
end
end

if 0
warning('off')
bigModelTable = bigModelOrig;
cpdIDs = cpdData(:,1); cpdNames = cpdData(:,2); cpdAbbrvs = cpdData(:,3);
rxnIDs = rxnData(:,1); rxnNames = rxnData(:,2); equations = rxnData(:,7);
equations = cellfun(@(x) strrep(strrep(x,'(',''),')',''), equations, 'UniformOutput',0);
badIdxs = find(strcmp(equations,'<=> '));
for i=1:length(rxnIDs)
    if mod(i,100)==0
        disp(i)
        disp(length(bigModelTable.rxns));
    end
    if ~any(strcmp(rxnIDs{i},bigModel.rxns)) && isKey(rxnsToECsTable,rxnIDs{i}) && any(ismember(GreenblumEC,rxnsToECsTable(rxnIDs{i}))) && ~any(i==badIdxs)
        [metList coeffList revFlag] = parseRxnFormula(equations{i});
        bigModelTable = addReaction(bigModelTable,rxnIDs{i}, metList,coeffList,revFlag);
    end
end
end

rxnsToExpressObese = mapExpToRxns(ECsToRxns,[baseDir filesep 'testObese.txt']);
rxnsToExpressNorm = mapExpToRxns(ECsToRxns,[baseDir filesep 'testNorm.txt']);
rxnsToDiffExp = containers.Map;
rxnsToExpressNormKeys = keys(rxnsToExpressNorm);
for i=1:length(rxnsToExpressNormKeys)
    expKey = rxnsToExpressNormKeys{i};
    if isKey(rxnsToExpressObese,expKey)
        rxnsToDiffExp(expKey) = abs(log2( (rxnsToExpressObese(expKey)/37)/(rxnsToExpressNorm(expKey)/87) ));
    end
    i
end
if 0
connMatrixTable = makeConnMatrix(bigModelTable);
rxnMatrix = zeros(length(bigModelTable.rxns),length(bigModelTable.rxns));
[~, topMetIdxs] = sort(sum(bigModelTable.S~=0,2),'descend');
topMetIdxs = topMetIdxs(1:floor(length(topMetIdxs)*.001));
for i=1:length(bigModelTable.rxns)
    metIdxs = find(bigModelTable.S(:,i)~=0);
    for j=1:length(metIdxs)
        if ~any(metIdxs(j)==topMetIdxs)
            connRxnIdxs = find(bigModelTable.S(metIdxs(j),:));
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
outputDir= baseDir;
makeBar(xvals,yvals,titleString,outputDir,'ylabelString', 'Abs Log Diff Expression','xlabelString','Centrality','isScatter',1);
end

KEGGmets = bigModelTableTest.mets;
KEGGmets = cellfun(@(x) strrep(x,'cpd','C'),KEGGmets,'UniformOutput',0);
for i=1:length(KEGGmets)
    braceIdx = regexp(KEGGmets{i},'[');
    currKEGGmet = KEGGmets{i};
    if ~isempty(braceIdx)
        KEGGmets{i} = currKEGGmet(1:braceIdx-1);
    end
end
metField = {}; stoichField = []; rxnField = {}; subSystemField = {}; KEGGField = {}; longMetField = {}; longRxnField = {};
for i=1:length(bigModelTableTest.rxns)
    metIdxs = find(bigModelTableTest.S(:,i));
    for j=1:length(metIdxs)
        rxnField{end+1} = bigModelTableTest.rxns{i};
        stoichField(end+1) = full(bigModelTableTest.S(metIdxs(j),i));
        metField{end+1} = bigModelTableTest.mets{metIdxs(j)};
        subSystemField{end+1} = bigModelTableTest.subSystems{i};
        KEGGField{end+1} = KEGGmets{metIdxs(j)};
    end
end
writeData({rxnField,stoichField,metField,subSystemField,KEGGField},[baseDir filesep 'bigModelTableTest.txt'],'\t');

if 0
bigModelTableTest = addMustEx(bigModelTable);
[bigModelTableTest, ~, addedModel] = mergeSmallest(bigModelTableTest,modelNamesToModelsValues,1);
[biomassMets biomassCoeffs] = makeMergedBiomass({addedModel},bigModelTableTest.mets);
%[sols, bigModelTableTestBio] = testBiomass(bigModelTableTest,biomassMets,biomassCoeffs,0);
bigModelTableTest = addReaction(bigModelTableTest,'BIOMASS', {bigModelTableTest.mets{end}},[-1],1);
bigModelTableTest.c = zeros(length(bigModelTableTest.rxns),1);
bigModelTableTest = changeObjective(bigModelTableTest,'BIOMASS');
expressionIDs = rxnsToExpressNormKeys;
for i=1:length(expressionIDs)
    expressionData(i) = rxnsToExpressNorm(expressionIDs{i});
end
expressionSDs = ones(size(expressionData));
for i=1:length(bigModelTableTest.rxns)
    bigModelTableTest.rxns{i} = strrep(strrep(bigModelTableTest.rxns{i},')',''),'(','');
end
bigModelTableTest.genes = bigModelTableTest.rxns;
bigModelTableTest.grRules = bigModelTableTest.rxns;
bigModelTableTest.rules = bigModelTableTest.rxns;
bigModelTableTest.rxnGeneMat = eye(length(bigModelTableTest.rxns));
bigModelTableTest = changeObjective(bigModelTableTest,'BIOMASS');
bigModelTableTest.description = 'bigModel';
falconfluxes = runFluxMethod(expressionData,expressionIDs,'test',bigModelTableTest,'FALCON',expressionSDs);
end
if 0
bigModelTableTest = addMustEx(bigModelTable);
percentageProd = 0;
minRxnsAdded = 0;
while percentageProd < 1
    sols = testBiomass(bigModelTableTest,biomassMets,biomassCoeffs,0);
    producesAll = cellfun(@(x) x.f,sols);
    percentageProd = sum(producesAll > 0)/length(producesAll);
    disp(percentageProd)
    disp(minRxnsAdded)
    if percentageProd < 1
        sortNum = 1;
        [bigModelTableTestCand, minRxnsAddedCand, addedModel] = mergeSmallest(bigModelTableTest,modelNamesToModelsValues,sortNum);
        
        metToFind = biomassMets{find(producesAll==0,1)}; hasMet = 0;
        while ~hasMet
            biomassIdxs = find(cellfun(@(x) ~isempty(regexpi(x,'(biomass|grow)')),addedModel.rxnNames));
            if ~isempty(biomassIdxs)
                for j=1:length(biomassIdxs)
                    biomassMetIdxs = find(addedModel.S(:,biomassIdxs(j))~=0);
                    for k=1:length(biomassMetIdxs)
                        if any(strcmp(addedModel.mets{biomassMetIdxs(k)},metToFind))
                            hasMet=1;
                        end
                    end
                end
            end
            if ~hasMet
                sortNum=sortNum+1;
                [bigModelTableTestCand, minRxnsAddedCand, addedModel] = mergeSmallest(bigModelTableTest,modelNamesToModelsValues,sortNum);
            end        
        end
        
        percentageProdCand = percentageProd;
        while percentageProdCand==percentageProd
            solsCand = testBiomass(bigModelTableTestCand,biomassMets,biomassCoeffs,0);
            producesAllCand = cellfun(@(x) x.f,solsCand);
            percentageProdCand = sum(producesAllCand > 0)/length(producesAllCand);
            if percentageProdCand==percentageProd
                sortNum = sortNum+1;
                [bigModelTableTestCand, minRxnsAddedCand] = mergeSmallest(bigModelTableTest,modelNamesToModelsValues,sortNum);
            end
        end
        bigModelTableTest = bigModelTableTestCand;
        minRxnsAdded = minRxnsAddedCand;
    end
end
end

