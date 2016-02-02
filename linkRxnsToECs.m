configSEED;
if 1
if 1
filenames = dir(modelsDir);
rxnsToECs = containers.Map;
ECsToRxns = containers.Map;
bigModel = makeEmptyModel();
bigModels = {};
if ~exist('modelNamesToModels','var')
    modelNamesToModels = containers.Map;
end
end
for i=1:length(filenames)
    if ~isempty(regexp(filenames(i).name,'.tsv'))
        modelName = filenames(i).name; modelName = modelName(1:regexp(modelName,'.tsv')-1);
        if 1
        status=system(sprintf([baseDir filesep 'makeTemp.sh %s %s'],[modelsDir filesep modelName '.tsv'],[modelsDir filesep modelName '.temp']));
        disp(modelName);
        FI = fopen([modelsDir filesep modelName '.temp']);
        line = fgetl(FI);
        count = 0;
        while line ~= -1
            if mod(count,25)==0
                disp(modelName)
                disp(count)
                disp(line)
            end
            count = count+1;
            [regex1 regex2] = regexp(line,'(\d|-)+\.(\d|-)+\.(\d|-)+\.(\d|-)+');
            if ~isempty(regex1) 
                ECNums = arrayfun(@(x,y) line(x:y), regex1,regex2, 'UniformOutput',0);
                words = strsplit(line,',');
                rxnName = words{1};
                [rxnsToECs ECsToRxns] = updateTwoMaps(rxnsToECs, ECsToRxns, rxnName, ECNums);
            end
            line = fgetl(FI);
        end
        fclose(FI);
        end
        
        if count > 10 && ~strcmp(modelName,'iAbaylyiv4') && ~strcmp(modelName,'iSB619')
            if ~isKey(modelNamesToModels,strrep(modelName,'.','_'))
                if isempty(regexp(modelName,'^i.*$'))
                    modelTemp = readCbModel([modelsDir filesep modelName '.xml']);
                else
                    modelTemp = readCbModel([modelsDir filesep modelName '_4.xml']);
                end
                modelNamesToModels(strrep(modelName,'.','_')) = modelTemp;
            else
                modelTemp = modelNamesToModels(strrep(modelName,'.','_'));
            end
            for j=1:length(modelTemp.rxns)
                if ~any(strcmp(bigModel.rxns,modelTemp.rxns{j})) && isKey(rxnsToECs,modelTemp.rxns{j}) && any(ismember(GreenblumEC,rxnsToECs(modelTemp.rxns{j})))
                    bigModel = mergeModels(bigModel,modelTemp,modelTemp.rxns{j});
                    bigModel = checkModelDims(bigModel);
                end
            end
            bigModels{end+1} = bigModel;
        end
    end
end
end

bigModelOrig = bigModel;

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
outputDir= baseDir;
makeBar(xvals,yvals,titleString,outputDir,'ylabelString', 'Total Expression','xlabelString','Centrality','isScatter',1);