if 1
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
for i=1:length(modelNames)
    for j=1:length(modelNames)
	modelNamesArr{end+1,1} = modelNames{i};
        modelNamesArr{end,2} = modelNames{j};
    end
end

flexibleBiomass = 1;
    
for i=1:length(modelNames)
    model = modelNamesToModels(modelNames{i});
    model.description = model.modelName;
    model.rxnECNums = {};
    for k=1:length(model.rxns)
	model.rxnECNums{end+1} = {};
    end
    model = fluxModelFunc(model);
    model = assignSortedBiom(model);

    if flexibleBiomass
        cIdxs = find(model.c);
        for j=1:length(cIdxs)
	    jthBiomMetIdxs = find(model.S(:,cIdxs(j))<0);
            for k=1:length(jthBiomMetIdxs);
                model = addReaction(model,['FLEX_BIOM_' model.mets{jthBiomMetIdxs(k)}],{model.mets{jthBiomMetIdxs(k)}},[-1],0,-1000,0);
            end
	    model.lb(cIdxs(j)) = 0;
            model.ub(cIdxs(j)) = 0;
        end
    end

    if ~all(model.c==0)
	cIdxs = find(model.c);
        someMistake = 0;
        try
	    if length(cIdxs)==1
		pseudoFluxDist1 = pseudoMapPicrust(model,1,1,model.rxnNames{cIdxs});
		pseudoFluxDist2 = pseudoMapPicrust(model,0,1,model.rxnNames{cIdxs});
	    else
		pseudoFluxDist1 = pseudoMapPicrust(model,1,1,model.rxnNames{cIdxs(1)});
		pseudoFluxDist2 = pseudoMapPicrust(model,0,1,model.rxnNames{cIdxs(1)});
	    end
	catch
	    someMistake = 1;
	end
    end

    if someMistake==0 && ~all(isnan(pseudoFluxDist1)) && ~all(pseudoFluxDist1==0)
        outFI = fopen([outputDir1 filesep num2str(i) 'SpeciesAlone.txt'],'w');
	for k=1:length(pseudoFluxDist1)
	    fprintf(outFI,sprintf('%f %f\n',pseudoFluxDist1(k),pseudoFluxDist2(k)));
	end
	fclose(outFI);
    end
end
end

if 0
parfor z=1:length(modelNames)*length(modelNames)
    %parfor j=1:length(modelNames)
	%z = 2*(i*length(modelNames)+j);
        %count = count+1;
        if exist([outputDir1 filesep num2str(z) '.txt'],'file')
            continue;
        end
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
            %pseudoFluxDistArr1{z} = pseudoFluxDist1;
            %pseudoFluxDistArr2{z} = pseudoFluxDist2;
            %disp(sum(cellfun(@(x) ~isempty(x), pseudoFluxDistArr1)))
            outFI = fopen([outputDir1 filesep num2str(z) '.txt'],'w');
            for k=1:length(pseudoFluxDist1)
	        fprintf(outFI,sprintf('%f %f\n',pseudoFluxDist1(k),pseudoFluxDist2(k)));
            end
            fclose(outFI);
	    %save([outputDir1 filesep num2str(z) '.mat'],'pseudoFluxDist1','pseudoFluxDist2');
        end
    %end
end
%end
end
