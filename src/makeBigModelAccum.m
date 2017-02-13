configSEED;
outputDir1 = [outputDir filesep 'makeBigModelAccum'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1);
end
filenames = dir(modelsDir);
filenames = filenames(randperm(length(filenames)));
rxnsToECsAccum = containers.Map;
ECsToRxnsAccum = containers.Map;
bigModelAccum = makeEmptyModel();
bigModelsAccum = {};
reloadModels = 0;
IGNOREDASH = 1;
if ~exist('modelNamesToModels','var') || reloadModels
    modelNamesToModels = containers.Map;
end
for i=1:length(filenames)
    if ~isempty(regexp(filenames(i).name,'.tsv'))
        modelName = filenames(i).name; modelName = modelName(1:regexp(modelName,'.tsv')-1);
        status=system(sprintf([baseDir filesep 'src' filesep 'makeTemp.sh %s %s'],[modelsDir filesep modelName '.tsv'],[modelsDir filesep modelName '.temp']));
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
            if ~IGNOREDASH
                [regex1 regex2] = regexp(line,'(\d|-)+\.(\d|-)+\.(\d|-)+\.(\d|-)+');
            else
                [regex1 regex2] = regexp(line,'(\d)+\.(\d)+\.(\d)+\.(\d)+');
            end
            if ~isempty(regex1) 
                ECNums = arrayfun(@(x,y) line(x:y), regex1,regex2, 'UniformOutput',0);
                words = strsplitYiping(line,',');
                rxnName = words{1};
                [rxnsToECsAccum ECsToRxnsAccum] = updateTwoMaps(rxnsToECsAccum, ECsToRxnsAccum, rxnName, ECNums);
            end
            line = fgetl(FI);
        end
        fclose(FI);
        
        if count > 10 && ~strcmp(modelName,'iAbaylyiv4') && ~strcmp(modelName,'iSB619')
            if ~isKey(modelNamesToModels,strrep(modelName,'.','_'))
                if isempty(regexp(modelName,'^i.*$'))
                    modelTemp = readCbModel([modelsDir filesep modelName '.xml']);
                else
                    modelTemp = readCbModel([modelsDir filesep modelName '_3.xml']);
                end
                modelNamesToModels(strrep(modelName,'.','_')) = modelTemp;
            else
                modelTemp = modelNamesToModels(strrep(modelName,'.','_'));
            end
            %strip off compartments for compatibility with bigModelTable
            for j=1:length(modelTemp.mets)
                bracketIdx = regexp(modelTemp.mets{j},'\[c\]');
                if ~isempty(bracketIdx)
                    modelTemp.mets{j} = modelTemp.mets{j}(1:bracketIdx-1);
                end
            end
            for j=1:length(modelTemp.rxns)
                if ~any(strcmp(bigModelAccum.rxns,modelTemp.rxns{j})) && (~isempty(regexp(modelTemp.rxnNames{j},'EX')) || ~isempty(regexpi(modelTemp.rxnNames{j},'transport')))
                    bigModelAccum = mergeModels(bigModelAccum,modelTemp,modelTemp.rxns{j});
                    bigModelAccum.rxnECNums{strcmp(bigModelAccum.rxns,modelTemp.rxns{j})} = {};
                    bigModelAccum = checkModelDims(bigModelAccum);
                end
                if ~any(strcmp(bigModelAccum.rxns,modelTemp.rxns{j})) && isKey(rxnsToECsAccum,modelTemp.rxns{j}) && any(ismember(GreenblumEC,rxnsToECsAccum(modelTemp.rxns{j})))
                    bigModelAccum = mergeModels(bigModelAccum,modelTemp,modelTemp.rxns{j});
                    bigModelAccum.rxnECNums{strcmp(bigModelAccum.rxns,modelTemp.rxns{j})} = rxnsToECsAccum(modelTemp.rxns{j});
                    bigModelAccum = checkModelDims(bigModelAccum);
                end
            end
            bigModelsAccum{end+1} = bigModelAccum;
        end
    end
end
bigModelsAccumSubsystems = cellfun(@(x) x.subSystems,bigModelsAccum,'UniformOutput',0);
save([outputDir1 filesep 'extra.mat'],'bigModelsAccumSubsystems');
save([outputDir1 filesep 'makeBigModelAccum.mat'],'bigModelAccum', ...
     'modelNamesToModels','rxnsToECsAccum','ECsToRxnsAccum','bigModelsAccum');

