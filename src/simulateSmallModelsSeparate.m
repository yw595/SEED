outputDir1 = [outputDir filesep 'simulateSmallModelsSeparate'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1);
end
if 0

FBAGapMatrix = [];
newFBAMatrix = [];
newFBAMatrixPureComp = [];
numIntersectMatrix = zeros(length(modelNames),length(modelNames));
numIntersectMatrixAdd = zeros(length(modelNames),length(modelNames));
numIntersectMatrixMets = {};
numIntersectMatrixMetsAdd = {};
blockOnlyPotentialShad = 1;
for i=1:length(modelNames)
    for j=1:length(modelNames)
	numIntersectMatrixMets{i,j} = {};
        numIntersectMatrixMetsAdd{i,j} = {};
    end
end
for i=1:length(modelNames)
    for j=1:length(modelNames)
        disp([num2str(i) ' ' num2str(j)]);
        FBAGapMatrix(i,j) = allBiomassRates(i)-allBiomassRates(j);

        model1 = modelNamesToModels(modelNames{i});
        model2 = modelNamesToModels(modelNames{j});
        dist1 = allBiomassDists{i};
        dist2 = allBiomassDists{j};
        involvedMetKEGGs1 = {};
        involvedIdxs1 = [];
        involvedMetKEGGs1Add = {};
        involvedIdxs1Add = [];
        for k=1:length(dist1)
            if strcmp(model1.subSystems{k},'') && dist1(k)<0
                metKEGGToAdd = model1.metKEGGs(model1.S(:,k)~=0);
                involvedMetKEGGs1{end+1} = metKEGGToAdd;
                involvedIdxs1(end+1) = k;
            end
            if strcmp(model1.subSystems{k},'') && dist1(k)>0
                metKEGGToAdd = model1.metKEGGs(model1.S(:,k)~=0);
                involvedMetKEGGs1Add{end+1} = metKEGGToAdd;
                involvedIdxs1Add(end+1) = k;
            end
        end
        involvedMetKEGGs2 = {};
        involvedIdxs2 = [];
        involvedMetKEGGs2Add = {};
        involvedIdxs2Add = [];
        for k=1:length(dist2)
            if strcmp(model2.subSystems{k},'') && dist2(k)<0
                metKEGGToAdd = model2.metKEGGs(model2.S(:,k)~=0);
                involvedMetKEGGs2{end+1} = metKEGGToAdd;
                involvedIdxs2(end+1) = k;
            end
            if strcmp(model2.subSystems{k},'') && dist2(k)>0
                metKEGGToAdd = model2.metKEGGs(model2.S(:,k)~=0);
                involvedMetKEGGs2Add{end+1} = metKEGGToAdd;
                involvedIdxs2Add(end+1) = k;
            end
        end
        for k=1:length(involvedMetKEGGs1)
            if length(involvedMetKEGGs1{k})==1
                intersectIdx = find(cellfun(@(x) length(x)==1 && strcmp(involvedMetKEGGs1{k}{1},x{1}), involvedMetKEGGs2));
                if length(intersectIdx)==1
                    %length(intersectIdx)
                    %disp([num2str(i) ' ' num2str(j) ' ' num2str(k)]);
                    totalFlux = dist1(involvedIdxs1(k))+dist2(involvedIdxs2(intersectIdx));
                    if totalFlux > (model1.lb(involvedIdxs1(k))+model2.lb(involvedIdxs2(intersectIdx)))/2 || blockOnlyPotentialShad==0
                        model1.lb(k) = totalFlux/2;
                        model2.lb(intersectIdx) = totalFlux/2;
                    end
                    numIntersectMatrix(i,j) = numIntersectMatrix(i,j)+1;
		    tempArr = numIntersectMatrixMets{i,j};
                    metKEGGTemp = involvedMetKEGGs1{k}{1};
                    metTemp = model1.metNames{strcmp(model1.metKEGGs,metKEGGTemp)};
                    tempArr{end+1} = metTemp;
		    numIntersectMatrixMets{i,j} = tempArr;
                end
            end
        end
        model1PureComp = model1;
        model2PureComp = model2;
        for k=1:length(involvedMetKEGGs1Add)
            if length(involvedMetKEGGs1Add{k})==1
                intersectIdx = find(cellfun(@(x) length(x)==1 && strcmp(involvedMetKEGGs1Add{k}{1},x{1}), involvedMetKEGGs2));
                if length(intersectIdx)==1
                    totalFlux = dist1(involvedIdxs1Add(k));
                    model2.lb(intersectIdx) = model2.lb(intersectIdx)-totalFlux;
                    numIntersectMatrixAdd(i,j) = numIntersectMatrixAdd(i,j)+1;
		    tempArr = numIntersectMatrixMetsAdd{i,j};
                    metKEGGTemp = involvedMetKEGGs1Add{k}{1};
                    metTemp = model1.metNames{strcmp(model1.metKEGGs,metKEGGTemp)};
                    tempArr{end+1} = metTemp;
		    numIntersectMatrixMetsAdd{i,j} = tempArr;
                end
            end
        end
        for k=1:length(involvedMetKEGGs2Add)
            if length(involvedMetKEGGs2Add{k})==1
                intersectIdx = find(cellfun(@(x) length(x)==1 && strcmp(involvedMetKEGGs2Add{k}{1},x{1}), involvedMetKEGGs1));
                if length(intersectIdx)==1
                    totalFlux = dist2(involvedIdxs2Add(k));
                    model1.lb(intersectIdx) = model1.lb(intersectIdx)-totalFlux;
                    numIntersectMatrixAdd(i,j) = numIntersectMatrixAdd(i,j)+1;
		    tempArr = numIntersectMatrixMetsAdd{i,j};
                    metKEGGTemp = involvedMetKEGGs2Add{k}{1};
                    metTemp = model2.metNames{strcmp(model2.metKEGGs,metKEGGTemp)};
                    tempArr{end+1} = metTemp;
		    numIntersectMatrixMetsAdd{i,j} = tempArr;
                end
            end
        end

        newBiom1 = examineBiomassFactored(model1);
        newBiom2 = examineBiomassFactored(model2);
        newBiom1PureComp = examineBiomassFactored(model1PureComp);
        newBiom2PureComp = examineBiomassFactored(model2PureComp);
        newFBAMatrix(i,j,1) = newBiom1; newFBAMatrix(i,j,2) = newBiom2;
        newFBAMatrixPureComp(i,j,1) = newBiom1PureComp; newFBAMatrixPureComp(i,j,2) = newBiom2PureComp;
    end
end
end

save([outputDir1 filesep 'simulateSmallModelsSeparate.mat'],'FBAGapMatrix','newFBAMatrix','newFBAMatrixPureComp','numIntersectMatrix','numIntersectMatrixAdd','numIntersectMatrixMets','numIntersectMatrixMetsAdd');

if 1
modelsAddMet = {};
for i=1:length(modelNamesShort)
    modelsAddMet{i} = addMetSubsystems(modelNamesToModels(modelNames{i}));
end
end

if 0
numIntersectMatrixMetsArr1 = {};
numIntersectMatrixMetsArr2 = [];
numIntersectMatrixMetsArr1Add = {};
numIntersectMatrixMetsArr2Add = [];
numIntersectMatrixMetsSubArr1 = {};
numIntersectMatrixMetsSubArr2 = [];
numIntersectMatrixMetsSubArr1Add = {};
numIntersectMatrixMetsSubArr2Add = [];
for i=1:length(modelNamesShort)
    model1 = modelsAddMet{i};
    for j=1:length(modelNamesShort)
	model2 = modelsAddMet{j};
	tempArr = numIntersectMatrixMets{i,j};
        for k=1:length(tempArr)
	    if sum(strcmp(numIntersectMatrixMetsArr1,tempArr{k}))==0
	        numIntersectMatrixMetsArr1{end+1} = tempArr{k};
                numIntersectMatrixMetsArr2(end+1) = 0;
            end
            numIntersectMatrixMetsArr2(strcmp(numIntersectMatrixMetsArr1,tempArr{k})) = numIntersectMatrixMetsArr2(strcmp(numIntersectMatrixMetsArr1,tempArr{k}))+1;
            if sum(strcmp(model1.metNames,tempArr{k}))~=0 && sum(strcmp(model2.metNames,tempArr{k}))~=0
	       idx1 = find(strcmp(model1.metNames,tempArr{k}));
               idx1 = idx1(end);
               idx2 = find(strcmp(model2.metNames,tempArr{k}));
               idx2 = idx2(end);
	       sub1 = model1.metSubsystems{idx1};
               if sum(strcmp(sub1,numIntersectMatrixMetsSubArr1))==0
	            numIntersectMatrixMetsSubArr1{end+1} = sub1;
		    numIntersectMatrixMetsSubArr2(end+1) = 0;
		end
		numIntersectMatrixMetsSubArr2(strcmp(numIntersectMatrixMetsSubArr1,sub1)) = numIntersectMatrixMetsSubArr2(strcmp(numIntersectMatrixMetsSubArr1,sub1))+1;
		sub2 = model2.metSubsystems{idx2};
		if sum(strcmp(sub2,numIntersectMatrixMetsSubArr1))==0
		    numIntersectMatrixMetsSubArr1{end+1} = sub2;
		    numIntersectMatrixMetsSubArr2(end+1) = 0;
		end
		numIntersectMatrixMetsSubArr2(strcmp(numIntersectMatrixMetsSubArr1,sub2)) = numIntersectMatrixMetsSubArr2(strcmp(numIntersectMatrixMetsSubArr1,sub2))+1;
            end
        end
	tempArrAdd = numIntersectMatrixMetsAdd{i,j};
        for k=1:length(tempArrAdd)
	    if sum(strcmp(numIntersectMatrixMetsArr1Add,tempArrAdd{k}))==0
	        numIntersectMatrixMetsArr1Add{end+1} = tempArrAdd{k};
                numIntersectMatrixMetsArr2Add(end+1) = 0;
            end
            numIntersectMatrixMetsArr2Add(strcmp(numIntersectMatrixMetsArr1Add,tempArrAdd{k})) = numIntersectMatrixMetsArr2Add(strcmp(numIntersectMatrixMetsArr1Add,tempArrAdd{k}))+1;
            if sum(strcmp(model1.metNames,tempArrAdd{k}))~=0 && sum(strcmp(model2.metNames,tempArrAdd{k}))~=0
	       idx1 = find(strcmp(model1.metNames,tempArrAdd{k}));
               idx1 = idx1(end);
               idx2 = find(strcmp(model2.metNames,tempArrAdd{k}));
               idx2 = idx2(end);
	       sub1 = model1.metSubsystems{idx1};
               if sum(strcmp(sub1,numIntersectMatrixMetsSubArr1Add))==0
	            numIntersectMatrixMetsSubArr1Add{end+1} = sub1;
		    numIntersectMatrixMetsSubArr2Add(end+1) = 0;
		end
		numIntersectMatrixMetsSubArr2Add(strcmp(numIntersectMatrixMetsSubArr1Add,sub1)) = numIntersectMatrixMetsSubArr2Add(strcmp(numIntersectMatrixMetsSubArr1Add,sub1))+1;
		sub2 = model2.metSubsystems{idx2};
		if sum(strcmp(sub2,numIntersectMatrixMetsSubArr1Add))==0
		    numIntersectMatrixMetsSubArr1Add{end+1} = sub2;
		    numIntersectMatrixMetsSubArr2Add(end+1) = 0;
		end
		numIntersectMatrixMetsSubArr2Add(strcmp(numIntersectMatrixMetsSubArr1Add,sub2)) = numIntersectMatrixMetsSubArr2Add(strcmp(numIntersectMatrixMetsSubArr1Add,sub2))+1;
            end
        end
    end
    %nonsense = nonsense+1;
end
[numIntersectMatrixMetsArr2,sortIdxs] = sort(numIntersectMatrixMetsArr2,'descend');
numIntersectMatrixMetsArr1 = addIdxStrings(numIntersectMatrixMetsArr1(sortIdxs));
[numIntersectMatrixMetsArr2Add,sortIdxs] = sort(numIntersectMatrixMetsArr2Add,'descend');
numIntersectMatrixMetsArr1Add = addIdxStrings(numIntersectMatrixMetsArr1Add(sortIdxs));
[numIntersectMatrixMetsSubArr2,sortIdxs] = sort(numIntersectMatrixMetsSubArr2,'descend');
numIntersectMatrixMetsSubArr1 = addIdxStrings(numIntersectMatrixMetsSubArr1(sortIdxs));
[numIntersectMatrixMetsSubArr2Add,sortIdxs] = sort(numIntersectMatrixMetsSubArr2Add,'descend');
numIntersectMatrixMetsSubArr1Add = addIdxStrings(numIntersectMatrixMetsSubArr1Add(sortIdxs));
writeData({numIntersectMatrixMetsArr1,numIntersectMatrixMetsArr2},[transferDir filesep 'numCompetingMets.txt'],'\t',{'met','numintersect'});
writeData({numIntersectMatrixMetsArr1Add,numIntersectMatrixMetsArr2Add},[transferDir filesep 'numCooperatingMets.txt'],'\t',{'met','numintersect'});
writeData({numIntersectMatrixMetsSubArr1,numIntersectMatrixMetsSubArr2},[transferDir filesep 'numCompetingSubs.txt'],'\t',{'sub','numintersect'});
writeData({numIntersectMatrixMetsSubArr1Add,numIntersectMatrixMetsSubArr2Add},[transferDir filesep 'numCooperatingSubs.txt'],'\t',{'sub','numintersect'});
end

justSEEDData = textscan(fopen('/mnt/vdb/home/ubuntu2/justSEEDDists.txt'),'%s%s%s','Delimiter','|','HeaderLines',0);
justSEEDMatrix = [];
for i=1:length(justSEEDData{1})
    justSEEDMatrix(strcmp(justSEEDData{1}{i},modelNamesShort),strcmp(justSEEDData{2}{i},modelNamesShort)) = str2num(justSEEDData{3}{i});
end

justSEEDArr = [];
FBAGapArr = [];
for i=1:length(modelNames)
    for j=1:length(modelNames)
	justSEEDArr(end+1) = justSEEDMatrix(i,j);
	FBAGapArr(end+1) = FBAGapMatrix(i,j);
    end
end
		 
modelArr1 = {}; modelArr2 = {};
diffFromOneMatrix = [];
for i=1:length(modelNames)
    for j=1:length(modelNames)
        diffFromOneMatrix(i,j) = max([allBiomassRates(j),allBiomassRates(i)])-newFBAMatrix(i,j,1);
    end
end
diffFromOneMatrix(diffFromOneMatrix<0) = - diffFromOneMatrix(diffFromOneMatrix<0);
pureCompDiffMatrix = newFBAMatrix-newFBAMatrixPureComp;

diffArr = [];
compDiffArr = [];
FBAGapArr = [];
interactionArr = [];
for i=1:length(modelNames)
    for j=1:length(modelNames)
        modelArr1{end+1} = modelNamesShort{i};
        modelArr2{end+1} = modelNamesShort{j};
        diffArr(end+1) = diffFromOneMatrix(i,j);
        compDiffArr(end+1) = (pureCompDiffMatrix(i,j)+pureCompDiffMatrix(i,j))/2;
        FBAGapArr(end+1) = FBAGapMatrix(i,j);
        interactionArr(end+1) = numIntersectMatrixAdd(i,j);
    end
end
writeData({modelArr1,modelArr2,-diffArr},[transferDir filesep 'diffFromOneMatrixWithCoop.txt'],'\t',{'species1','species2','diffFromOne'});
writeData({modelArr1,modelArr2,compDiffArr},[transferDir filesep 'pureCompDiffMatrix.txt'],'\t',{'species1','species2','pureCompDiff'});
writeData({modelArr1,modelArr2,interactionArr},[transferDir filesep 'numInteractionsAddMatrix.txt'],'\t',{'species1','species2','numInteractionsAdd'});
[FBAGapArrSorted sortIdxsGap] = sort(FBAGapArr,'descend');
writeData({modelArr1,modelArr2,FBAGapArr},[transferDir filesep 'FBAGapHeatmap.txt'],'\t',{'species1','species2','fbagap'});
writeData({justSEEDArr,FBAGapArr},[transferDir filesep 'justSEEDDistVsFBAGap.txt'],'\t',{'justseeddist','fbagap'});






