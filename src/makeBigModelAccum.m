configSEED;
outputDir1 = [outputDir filesep 'makeBigModelAccum'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1);
end
load([outputDir1 filesep 'makeBigModelAccumSafe.mat'],'bigModelAccum','modelNamesToModels','rxnsToECsAccum','ECsToRxnsAccum');
filenames = dir(modelsDir);
filenames = filenames(randperm(length(filenames)));
rxnsToECsAccum = containers.Map;
ECsToRxnsAccum = containers.Map;
bigModelAccum = makeEmptyModel();
bigModelAccum.metKEGGs = {};
bigModelsAccum = {};
temporaryFix = 1;
reloadModels = 0;
if temporaryFix==1
    reloadModels = 0;
end
IGNOREDASH = 1;
if ~exist('modelNamesToModels','var') || reloadModels
    modelNamesToModels = containers.Map;
end

if temporaryFix==1
    modelNames = keys(modelNamesToModels);
    for i=1:length(modelNames)
        tempModel = modelNamesToModels(modelNames{i});
        tempModel.metKEGGs = {};
        for k=1:length(tempModel.mets)
            reconcmet = tempModel.mets{k};
            if ~isempty(regexp(reconcmet,'\['))
                reconcmet = reconcmet(1:end-3);
            end
            matchIdx = find(strcmp(reconcmet,cpdIDs));
            if ~isempty(matchIdx)
                matchIdx = matchIdx(1);
                tempModel.metKEGGs{k} = cpdKEGGs{matchIdx};
            else
                tempModel.metKEGGs{k} = '';
            end
        end
        modelNamesToModels(modelNames{i}) = tempModel;
    end
end

for i=1:length(filenames)
    if ~isempty(regexp(filenames(i).name,'.tsv'))
        disp('HERE')
        disp(i)
        modelName = filenames(i).name; modelName = modelName(1:regexp(modelName,'.tsv')-1);
        status=system(sprintf([baseDir filesep 'src' filesep 'Utils/makeTemp.sh %s %s'],[modelsDir filesep modelName '.tsv'],[modelsDir filesep modelName '.temp']));
        disp(modelName);
        FI = fopen([modelsDir filesep modelName '.temp']);
        line = fgetl(FI);
        count = 0;
        while line ~= -1
            if mod(count,25)==0
                %disp(modelName)
                %disp(count)
                %disp(line)
            end
            count = count+1;
            if count > 25
                break
            end
            
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
        
        disp('HERE1')
        disp(i)
        disp(length(keys(modelNamesToModels)));
        if count > 10 && ~strcmp(modelName,'iAbaylyiv4') && ~strcmp(modelName,'iSB619')
            if ~isKey(modelNamesToModels,strrep(modelName,'.','_'))
                if isempty(regexp(modelName,'^i.*$'))
                    modelTemp = readCbModel([modelsDir filesep modelName '.xml']);
                else
                    modelTemp = readCbModel([modelsDir filesep modelName '_3.xml']);
                end

                if temporaryFix~=1
                    modelTemp.metKEGGs = {};
                    for k=1:length(modelTemp.mets)
                        reconcmet = modelTemp.mets{k};
                        if ~isempty(regexp(reconcmet,'\['))
                            reconcmet = reconcmet(1:end-3);
                        end
                        matchIdx = find(strcmp(reconcmet,cpdIDs));
                        if ~isempty(matchIdx)
                            matchIdx = matchIdx(1);
                            modelTemp.metKEGGs{k} = cpdKEGGs{matchIdx};
                        else
                           modelTemp.metKEGGs{k} = '';
                        end
                    end
                end
                modelNamesToModels(strrep(modelName,'.','_')) = modelTemp;
            else
                modelTemp = modelNamesToModels(strrep(modelName,'.','_'));
            end
            disp('HERE2')
            disp(i)
            %strip off compartments for compatibility with bigModelTable
            for j=1:length(modelTemp.mets)
                bracketIdx = regexp(modelTemp.mets{j},'\[c\]');
                if ~isempty(bracketIdx)
                    modelTemp.mets{j} = modelTemp.mets{j}(1:bracketIdx-1);
                end
            end
            disp('HERE3')
            %disp(i)
            %disp(length(bigModelAccum.rxns))
            %bigModelAccum.rxnECNums
            %bigModelAccum.metKEGGs
            if temporaryFix~=1
                for j=1:length(modelTemp.rxns)
                    %disp('line1')
                    if ~any(strcmp(bigModelAccum.rxns,modelTemp.rxns{j}))
                        %disp('line2')
                        bigModelAccum = mergeModels(bigModelAccum,modelTemp,modelTemp.rxns{j});
            %disp('line3')
                        bigModelAccum.rxnECNums{strcmp(bigModelAccum.rxns,modelTemp.rxns{j})} = {};
                        %disp('line4')
                        if (~isempty(regexp(modelTemp.rxnNames{j},'EX')) || ~isempty(regexpi(modelTemp.rxnNames{j},'transport')))
                        end
                        %disp('line5')
                        if isKey(rxnsToECsAccum,modelTemp.rxns{j}) && any(ismember(GreenblumEC,rxnsToECsAccum(modelTemp.rxns{j})))
                            bigModelAccum.rxnECNums{strcmp(bigModelAccum.rxns,modelTemp.rxns{j})} = rxnsToECsAccum(modelTemp.rxns{j});
                            %disp('line6')
                        end
                        %disp('line7')
                        bigModelAccum = checkModelDims(bigModelAccum);
                        %disp('line8')
                    end
                end
                %disp('line9')
                bigModelsAccum{end+1} = bigModelAccum;
                %disp('line10')
            end
        end

        disp('HERE4')
        disp(i)
        if count > 10 && ~strcmp(modelName,'iAbaylyiv4') && ~strcmp(modelName,'iSB619')
            modelTemp = modelNamesToModels(strrep(modelName,'.','_'));
            if ~isfield(modelTemp,'modelName')
                if exist('/home/ubuntu','dir')
                    xmlFI = fopen(modelTemp.description);
                else
                    xmlFI = fopen(['/home/fs01/yw595/' modelTemp.description(14:end)]);
                end
                line = fgetl(xmlFI);
                while line ~= -1
                    if ~isempty(regexp(line,'model.*name'))
                        t = regexp(line,'name=\"(.*)\"','tokens');
                        modelTemp.modelName = t{1}{1};
                        t{1}{1}
                    end
                    line = fgetl(xmlFI);
                end
                fclose(xmlFI);
            end
            modelNamesToModels(strrep(modelName,'.','_')) = modelTemp;
        end
    end
end
bigModelsAccumSubsystems = cellfun(@(x) x.subSystems,bigModelsAccum,'UniformOutput',0);
modelNames = keys(modelNamesToModels);
modelNamesShortList = cellfun(@(x) strsplit(modelNamesToModels(x).modelName,'_'), modelNames, 'UniformOutput', 0);
modelNamesShort = {};
for i=1:length(modelNamesShortList)
    temp = modelNamesShortList{i};
    modelNamesShort{end+1} = [temp{1} ' ' temp{2}];
end

save([outputDir1 filesep 'makeBigModelAccum.mat'],'bigModelAccum','bigModelsAccumSubsystems','modelNames','modelNamesToModels','modelNamesShort','rxnsToECsAccum','ECsToRxnsAccum');

