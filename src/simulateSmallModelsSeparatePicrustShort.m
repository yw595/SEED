configSEED;
outputDir1 = [outputDir filesep 'simulateSmallModelsSeparatePicrust'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1);
end

pseudoFluxDistMat1 = {};
pseudoFluxDistMat2 = {};
pseudoFluxDistArr1 = {};
pseudoFluxDistArr2 = {};
modelNamesArr = {};
for i=1:10%length(modelNames)
    for j=1:10%length(modelNames)
	modelNamesArr{end+1,1} = modelNames{i};
        modelNamesArr{end,2} = modelNames{j};
    end
end
parfor z=1:length(modelNames)*length(modelNames)
    %parfor j=1:length(modelNames)
	%z = 2*(i*length(modelNames)+j);
        %count = count+1;
        initCobraToolbox;
        %i = 1;%floor(z/length(modelNames))+1;
        %j = 1;%mod(z,length(modelNames))+1;
        model1 = modelNamesToModels(modelNamesArr{z,1});
        model2 = modelNamesToModels(modelNamesArr{z,2});
	model1.description = model1.modelName;
	model2.description = model2.modelName;
        model1.rxnECNums = {}; model2.rxnECNums = {};
        for k=1:length(model1.rxns)
	    model1.rxnECNums{end+1} = {};
        end
        for k=1:length(model2.rxns)
	    model2.rxnECNums{end+1} = {};
        end
	model1.rxns = cellfun(@(x) [x '_1'],model1.rxns,'UniformOutput',0);
	model2.rxns = cellfun(@(x) [x '_2'],model2.rxns,'UniformOutput',0);
	mergedModel = mergeModels(model1,model2,'',1);
        mergedModel = fluxModelFunc(mergedModel);
        %mergedModel = makeBigModelTableFunc(pseudoModelFormatPicrust(mergedModel),rxnData,equations,rxnIDs,rxnNames,GreenblumEC,cpdIDs,cpdKEGGs);

        mergedModel = assignSortedBiom(mergedModel);
	if ~all(mergedModel.c==0)
	    cIdxs = find(mergedModel.c);
            if length(cIdxs)==1
	        pseudoFluxDist1 = pseudoMapPicrust(mergedModel,1,1,mergedModel.rxnNames{cIdxs});
                pseudoFluxDist2 = pseudoMapPicrust(mergedModel,0,1,mergedModel.rxnNames{cIdxs});
            else
	        pseudoFluxDist1 = pseudoMapPicrust(mergedModel,1,1,mergedModel.rxnNames{cIdxs(1)});
                pseudoFluxDist2 = pseudoMapPicrust(mergedModel,0,1,mergedModel.rxnNames{cIdxs(1)});
            end
            %pseudoFluxDistMat1{i,j} = pseudoFluxDist1;
            %pseudoFluxDistMat2{i,j} = pseudoFluxDist2;
            pseudoFluxDistArr1{z} = pseudoFluxDist1;
            pseudoFluxDistArr2{z} = pseudoFluxDist2;
            %disp(sum(cellfun(@(x) ~isempty(x), pseudoFluxDistArr1)))
        end
    %end
end
%end
