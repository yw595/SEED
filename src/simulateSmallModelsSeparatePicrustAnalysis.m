if 1

outputDir1 = [outputDir filesep 'simulateSmallModelsSeparatePicrust'];
modelNamesArr = {};
pseudoBiomassMatNormal = -ones(length(modelNames),length(modelNames));
pseudoBiomassMatObese = -ones(length(modelNames),length(modelNames));
pseudoBiomassArrNormal = []; pseudoBiomassArrObese = [];
pseudoSubsMat = {};
for i=1:length(modelNames)
    for j=1:length(modelNames)
	modelNamesArr{end+1,1} = modelNames{i};
        modelNamesArr{end,2} = modelNames{j};
        pseudoSubsMat{i,j,1} = {};
        pseudoSubsMat{i,j,2} = {};
        pseudoSubsMat{i,j,3} = {};
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
    if exist([outputDir1 filesep num2str(i) 'SpeciesAlone.txt'],'file')
	inFI = fopen([outputDir1 filesep num2str(i) 'SpeciesAlone.txt']);
	dataFields = textscan(inFI,'%s%s','Delimiter',' ','HeaderLines',0);
	fclose(inFI);
	dataFields = [dataFields{:}];
	pseudoFluxDist1 = dataFields(:,1);
	pseudoFluxDist2 = dataFields(:,2);
	pseudoFluxDist1 = cellfun(@(x) str2num(x), pseudoFluxDist1);
	pseudoFluxDist2 = cellfun(@(x) str2num(x), pseudoFluxDist2);

        pseudoBiomassArrNormal(i) = picrustDetermineBiomass(model,pseudoFluxDist1,flexibleBiomass);
        pseudoBiomassArrObese(i) = picrustDetermineBiomass(model,pseudoFluxDist2,flexibleBiomass);
    else
        pseudoBiomassArrNormal(i) = nan;
        pseudoBiomassArrObese(i) = nan;
    end
end

oxsecreteIDs = {'C00003','C00004','C00005','C00006', ... %NADH and NADPH
                'C00051','C00127', ... %Glutathiones
		'C00072', ... %Ascorbate
		'C02477','C02483', ... %Tocopherols
		'C00386','C00884', ... %Carnosines
		'C00366', ... %Urate
		'C00486'}; %Bilirubin
oxsecreteMat = zeros(length(modelNames),length(modelNames),2,length(oxsecreteIDs));

z=0;
totalOxArr1 = [];
totalOxArr2 = [];
totalFermentArr1 = [];
totalFermentArr2 = [];
for i=1:length(modelNames)
    for j=1:length(modelNames)
	z = z+1;
	if ~exist([outputDir1 filesep num2str(z) '.txt'])
            continue;
        end
	z
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
	mergedModel = mergeModelsSEED(model1,model2,'',0);
        mergedModel = fluxModelFunc(mergedModel);
        %mergedModel = makeBigModelTableFunc(pseudoModelFormatPicrust(mergedModel),rxnData,equations,rxnIDs,rxnNames,GreenblumEC,cpdIDs,cpdKEGGs);

        mergedModel = assignSortedBiom(mergedModel);
		
	inFI = fopen([outputDir1 filesep num2str(z) '.txt']);
        dataFields = textscan(inFI,'%s%s','Delimiter',' ','HeaderLines',0);
	fclose(inFI);
	dataFields = [dataFields{:}];
	pseudoFluxDist1 = dataFields(:,1);
	pseudoFluxDist2 = dataFields(:,2);
	pseudoFluxDist1 = cellfun(@(x) str2num(x), pseudoFluxDist1);
	pseudoFluxDist2 = cellfun(@(x) str2num(x), pseudoFluxDist2);

        for l=1:length(oxsecreteIDs)
	    matchOxSecrete = find(strcmp(mergedModel.metKEGGs,oxsecreteIDs{l}));
	    if length(matchOxSecrete)==1
		matchOxSecreteRxns = find(mergedModel.S(matchOxSecrete,:)~=0);
		for k=1:length(matchOxSecreteRxns)
		    if (pseudoFluxDist1(matchOxSecreteRxns(k))>0 && mergedModel.S(matchOxSecrete,matchOxSecreteRxns(k))>0) || (pseudoFluxDist1(matchOxSecreteRxns(k))<0 && mergedModel.S(matchOxSecrete,matchOxSecreteRxns(k))<0)
			oxsecreteMat(i,j,1,l) = oxsecreteMat(i,j,1,l)+abs(pseudoFluxDist1(matchOxSecreteRxns(k)));
		    end
		    if (pseudoFluxDist2(matchOxSecreteRxns(k))>0 && mergedModel.S(matchOxSecrete,matchOxSecreteRxns(k))>0) || (pseudoFluxDist2(matchOxSecreteRxns(k))<0 && mergedModel.S(matchOxSecrete,matchOxSecreteRxns(k))<0)
			oxsecreteMat(i,j,2,l) = oxsecreteMat(i,j,2,l)+abs(pseudoFluxDist2(matchOxSecreteRxns(k)));
		    end
		end
	    end
        end

        [subLabels,subFluxSums,subFluxDiffSums] = segmentFluxBySubsystem(mergedModel,pseudoFluxDist1,1,pseudoFluxDist2,1);
        pseudoBiomassMatNormal(i,j) = picrustDetermineBiomass(mergedModel,pseudoFluxDist1,flexibleBiomass);
        pseudoBiomassMatObese(i,j) = picrustDetermineBiomass(mergedModel,pseudoFluxDist2,flexibleBiomass);
        pseudoSubsMat{i,j,1} = subLabels;
        pseudoSubsMat{i,j,2} = subFluxSums;
        pseudoSubsMat{i,j,3} = subFluxDiffSums;

        [totalOx1 totalFerment1] = measureOxFermFunc(mergedModel,pseudoFluxDist1);
        [totalOx2 totalFerment2] = measureOxFermFunc(mergedModel,pseudoFluxDist2);
        totalOxArr1(z) = totalOx1;
	totalOxArr2(z) = totalOx2;
	totalFermentArr1(z) = totalFerment1;
	totalFermentArr2(z) = totalFerment2;
    end
end
save([outputDir filesep 'simulateSmallModelsSeparatePicrustAnalysis' filesep 'simulateSmallModelsSeparatePicrustAnalysis.mat'],'modelNamesArr','pseudoBiomassMatNormal','pseudoBiomassMatObese','pseudoBiomassArrNormal','pseudoBiomassArrObese','pseudoSubsMat','totalOxArr1','totalOxArr2','totalFermentArr1','totalFermentArr2');
nonsense = nonsense+1;

end

justSEEDData = textscan(fopen('/mnt/vdb/home/ubuntu2/justSEEDDists.txt'),'%s%s%s','Delimiter','|','HeaderLines',0);
justSEEDMatrix = [];
for i=1:length(justSEEDData{1})
    justSEEDMatrix(strcmp(justSEEDData{1}{i},modelNamesShort),strcmp(justSEEDData{2}{i},modelNamesShort)) = str2num(justSEEDData{3}{i});
end
justSEEDArr = [];
normalArr = [];
species1Arr = {};
species2Arr = {};
for i=1:length(modelNames)
    for j=1:length(modelNames)
	%if pseudoBiomassMatNormal(i,j)~=-1 && pseudoBiomassMatObese(i,j)~=-1
	justSEEDArr(end+1) = justSEEDMatrix(i,j);
        normalArr(end+1) = pseudoBiomassMatNormal(i,j);
        species1Arr{end+1} = modelNamesShort{i};
        species2Arr{end+1} = modelNamesShort{j};
        %end
    end
end

totalSubsRanks = [];
totalSubsNames = {};
pseudoBiomassMatDiff = -ones(length(modelNames),length(modelNames));
for i=1:length(modelNames)
    for j=1:length(modelNames)
	if pseudoBiomassMatNormal(i,j)~=-1 && pseudoBiomassMatObese(i,j)~=-1
	    pseudoBiomassMatDiff(i,j) = pseudoBiomassMatNormal(i,j) - pseudoBiomassMatObese(i,j);
            subLabels = pseudoSubsMat{i,j,1};
            for k=1:length(subLabels)
	        occurIdx = strcmp(totalSubsNames,subLabels{k});
                if sum(occurIdx)==1
                    totalSubsRanks(occurIdx) = totalSubsRanks(occurIdx)+k;
                elseif sum(occurIdx)==0
                    totalSubsNames{end+1} = subLabels{k};
                    totalSubsRanks(end+1) = 1;
                end
	    end
        end
    end
end
[totalSubsRanks,sortIdxs] = sort(totalSubsRanks,'descend');
totalSubsNames = totalSubsNames(sortIdxs);
totalSubsNames = addIdxStrings(totalSubsNames);
writeData({totalSubsNames,totalSubsRanks},'/mnt/vdb/home/ubuntu2/mergedModelTotalSubs.txt','\t',{'subname','subrank'});
writeData({normalArr3,species1Arr,species2Arr},'/mnt/vdb/home/ubuntu2/mergedModelNormalBiomass.txt','\t',{'normalbiom','species1','species2'});
writeData({justSEEDArr,normalArr},'/mnt/vdb/home/ubuntu2/justSEEDVsSmallBiomass.txt','\t',{'justseeddist','smallbiomass'});
nonsense = nonsense+1;
		
pseudoBiomassMatNormal = [];
pseudoBiomassMatObese = [];
mergedSubsystemsArr = {};
totalExchangeMat = [];
if 1
mergedModelsArr = {};
for i=1:length(modelNames)
	for j=1:length(modelNames)
					model1 = modelNamesToModels(modelNames{i});
	model2 = modelNamesToModels(modelNames{j});
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
%mergedModelsArr{i,j} = mergeModels(model1,model2,'',1);

mergedModel = mergeModels(model1,model2,'',1);
pseudoFluxDist = pseudoFluxDistMatNormal{i,j};
totalExchangeMat(i,j) = sum(abs(pseudoFluxDist(strcmp(mergedModel.subSystems,'Exchange') | strcmp(mergedModel.subSystems,'Transport') | strcmp(mergedModel.subSystems,''))));
i
j
length(mergedModel.rxns)
end
end
end
if 0
for i=1:length(modelNames)
	for j=1:length(modelNames)
i
j
					model1 = modelNamesToModels(modelNames{i});
	model2 = modelNamesToModels(modelNames{j});
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
pseudoBiomassMatNormal(i,j) = picrustDetermineBiomass(mergedModel,pseudoFluxDistMatNormal{i,j},flexibleBiomass);
pseudoBiomassMatObese(i,j) = picrustDetermineBiomass(mergedModel,pseudoFluxDistMatObese{i,j},flexibleBiomass);
mergedSubsystemsArr{i,j} = mergedModel.subSystems;
end
end
end

if 0
for i=1:10%length(modelNames)
	for j=1:10%length(modelNames)
		%if all(size(pseudoFluxDistMatObese{i,j})==size(pseudoFluxDistMatNormal{i,j}))
		sum(abs(pseudoFluxDistMatObese{i,j}-pseudoFluxDistMatNormal{i,j}))
		%end
end
end
end
