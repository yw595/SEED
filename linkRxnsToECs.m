if 0
if 1
outputDir = '/mnt/extra/SEED';
filenames = dir(outputDir);
rxnsToECs = containers.Map;
ECsToRxns = containers.Map;
allRxnNames = containers.Map;
rxnsToMets = containers.Map;
rxnsToCoeffs = containers.Map;
FI = fopen([outputDir filesep 'GreenblumECs.txt']);
GreenblumEC = textscan(FI,'%s\n');
fclose(FI);
GreenblumEC = GreenblumEC{1};
bigModel = struct(); bigModel.S = []; bigModel.rxns = {}; bigModel.mets = {}; bigModel.lb = []; bigModel.ub = []; bigModel.rxnNames = {}; bigModel.metNames = {};
bigModels = {};
if ~exist('modelNamesToModels','var')
    modelNamesToModels = containers.Map;
end
end
for i=1:length(filenames)
    if ~isempty(regexp(filenames(i).name,'.tsv'))
        modelName = filenames(i).name; modelName = modelName(1:regexp(modelName,'.tsv')-1);
        if 1
        status=system(sprintf(['/home/ubuntu/MATLAB/SEED/makeTemp.sh %s %s'],[outputDir filesep modelName '.tsv'],[outputDir filesep modelName '.temp']));
        disp(modelName);
        FI = fopen([outputDir filesep modelName '.temp']);
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
                if ~isKey(rxnsToECs,rxnName)
                    rxnsToECs(rxnName) = ECNums;
                else
                    temp = rxnsToECs(rxnName);
                    if ~iscell(temp)
                        temp = {temp};
                    end
                    for j=1:length(ECNums)
                        temp{end+1} = ECNums{j};
                    end
                    temp = unique(temp);
                    rxnsToECs(rxnName) = temp;
                end
                for j=1:length(ECNums)
                    if ~isKey(ECsToRxns,ECNums{j})
                        ECsToRxns(ECNums{j}) = rxnName;
                    else
                        temp = ECsToRxns(ECNums{j});
                        if ~iscell(temp) 
                            temp = {temp};
                        end
                        temp{end+1} = rxnName;
                        temp = unique(temp);
                        ECsToRxns(ECNums{j}) = temp;
                    end
                end
            end
            line = fgetl(FI);
        end
        fclose(FI);
        end
        
        if count > 10 && ~strcmp(modelName,'iAbaylyiv4') && ~strcmp(modelName,'iSB619')
            if ~isKey(modelNamesToModels,strrep(modelName,'.','_'))
                if isempty(regexp(modelName,'^i.*$'))
                    modelTemp = readCbModel([outputDir filesep modelName '.xml']);
                else
                    modelTemp = readCbModel([outputDir filesep modelName '_4.xml']);
                end
                modelNamesToModels(strrep(modelName,'.','_')) = modelTemp;
            else
                modelTemp = modelNamesToModels(strrep(modelName,'.','_'));
            end
            for j=1:length(modelTemp.rxns)
                if ~any(strcmp(bigModel.rxns,modelTemp.rxns{j})) && isKey(rxnsToECs,modelTemp.rxns{j}) && any(ismember(GreenblumEC,rxnsToECs(modelTemp.rxns{j})))
                    bigModel.rxns{end+1} = modelTemp.rxns{j};
                    bigModel.rxnNames{end+1} = modelTemp.rxnNames{j};
                    mets = modelTemp.mets(modelTemp.S(:,j)~=0);
                    metNames = modelTemp.metNames(modelTemp.S(:,j)~=0);
                    coeffs = modelTemp.S(modelTemp.S(:,j)~=0,j);
                    for k=1:length(mets)
                        if ~any(strcmp(bigModel.mets,mets{k}))
                            bigModel.mets{end+1} = mets{k};
                            bigModel.metNames{end+1} = metNames{k};
                        end
                        bigModel.S(strcmp(bigModel.mets,mets{k}),strcmp(bigModel.rxns,modelTemp.rxns{j})) = coeffs(k);
                        bigModel.lb(strcmp(bigModel.rxns,modelTemp.rxns{j})) = modelTemp.lb(j);
                        bigModel.ub(strcmp(bigModel.rxns,modelTemp.rxns{j})) = modelTemp.ub(j);
                    end
                end
            end
            bigModels{end+1} = bigModel;
        end
    end
end
end

connMatrix = makeConnMatrix(bigModel);
cents = betweenness_centrality(sparse(connMatrix));

bigModel = addMustEx(bigModel);

allModels = keys(modelNamesToModels);
biomassMets = {}; biomassCoeffs = [];
for i=1:length(allModels);
    modelTemp = modelNamesToModels(allModels{i});
    biomassIdxs = find(cellfun(@(x) ~isempty(regexpi(x,'(biomass|grow)')),modelTemp.rxnNames));
    if ~isempty(biomassIdxs)
        for j=1:length(biomassIdxs)
            biomassMetIdxs = find(modelTemp.S(:,biomassIdxs(j))~=0);
            for k=1:length(biomassMetIdxs);
                if any(strcmp(modelTemp.mets{biomassMetIdxs(k)},bigModel.mets))
                    biomassMets{end+1} = modelTemp.mets{biomassMetIdxs(k)};
                    biomassCoeffs(end+1) = modelTemp.S(biomassMetIdxs(k),biomassIdxs(j));
                end
            end
        end
    end
end
[biomassMets uniqIdxs ~] = unique(biomassMets);
biomassCoeffs = biomassCoeffs(uniqIdxs);

modelNamesToModelsKeys = keys(modelNamesToModels);
for i=1:length(modelNamesToModelsKeys)
    modelTemp = modelNamesToModels(modelNamesToModelsKeys{i});
    diffnums(i) = length(setdiff(modelTemp.rxns,bigModel.rxns));
end
[~, minIdx] = min(diffnums);
bigModelAdded = bigModel;
bigModelAdded=mergeModels(bigModelAdded,modelNamesToModels(modelNamesToModelsKeys{minIdx}));

sols = {};
cumulative = 0;
for i=1:length(biomassMets)
    disp(i)
    if i==1 | ~cumulative
        bigModelBio = bigModelAdded;
        bigModelBio.lb = bigModelBio.lb';
        bigModelBio.ub = bigModelBio.ub';
        bigModelBio.rxns{end+1} = 'BIOMASS';
        bigModelBio.rxnNames{end+1} = 'BIOMASS';
        biomassIdx = length(bigModelBio.rxns);
    end
    bigModelBio.S(ismember(bigModelBio.mets,biomassMets{i}),biomassIdx) = biomassCoeffs(i);
    bigModelBio.lb(biomassIdx)=-1000;bigModelBio.ub(biomassIdx)=1000;
    bigModelBio.c = zeros(length(bigModelBio.rxns),1);
    %bigModelBio.c(end) = 1;
    bigModelBio = changeObjective(bigModelBio,'BIOMASS');
    sols{end+1} = optimizeCbModel(bigModelBio);
end

xvals = 1:length(sols); titleString = 'biomassMetProd'; yvals = []; outputDir= '/home/ubuntu/MATLAB/SEED';
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

xvals = 1:length(bigModels); titleString = 'numRxnsAdded'; yvals = cellfun(@(x) length(x.rxns), bigModels); outputDir= '/home/ubuntu/MATLAB/SEED';
xlabels = {};
for i=1:length(bigModels)
    xlabels{i} = num2str(i);
end
makeBar(xvals,yvals,titleString,outputDir,'ylabelString', 'numRxns','xlabelString','Model Number','xlabels',xlabels);

FI = fopen('/home/ubuntu/Downloads/testEC2.txt');
dataFields = textscan(FI,'%s%s','Delimiter',' ');
fclose(FI);
dataFields = [dataFields{:}];
ECs = keys(ECsToRxns);
rxnsToExpress = containers.Map;
for i=1:length(ECs)
    matchIdx = strcmp(dataFields(:,1),ECs{i});
    matchIdx = find(matchIdx);
    if matchIdx > 0
        rxns = ECsToRxns(ECs{i});
        if ~iscell(rxns)
            rxns = {rxns};
        end
        for j=1:length(rxns)
            rxnsToExpress(rxns{j}) = str2num(dataFields{matchIdx,2});
        end
    end
end
xvals = cents(1:length(bigModel.rxns)); titleString = 'Total Expression Vs. Centrality';
yvals = [];
for i=1:length(xvals)
    if isKey(rxnsToExpress,bigModel.rxns{i})
        yvals(i) = rxnsToExpress(bigModel.rxns{i});
    else
        yvals(i) = 0;
    end
end
outputDir= '/home/ubuntu/MATLAB/SEED';
makeBar(xvals,yvals,titleString,outputDir,'ylabelString', 'Total Expression','xlabelString','Centrality','isScatter',1);

if 0
ECs = keys(ECsToRxns);
for i=1:length(ECs)
    if any(strcmp(ECs{i},GreenblumEC))
        rxns = ECsToRxns(ECs{i});
        if ~iscell(rxns)
            rxns = {rxns};
        end
        for j=1:length(rxns)
            if ~any(strcmp(bigModel.rxns,rxns{j})) && isKey(rxnsToMets,rxns{j})
                bigModel.rxns{end+1} = rxns{j};
                mets = rxnsToMets(rxns{j});
                coeffs = rxnsToCoeffs(rxns{j});
                for k=1:length(mets)
                    if ~any(strcmp(bigModel.mets,mets{k}))
                        bigModel.mets{end+1} = mets{k};
                    end
                    bigModel.S(strcmp(bigModel.mets,mets{k}),strcmp(bigModel.rxns,rxns{j})) = coeffs(k);
                end
            end
        end
    end
end
end