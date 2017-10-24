function [FBAGapMatrix newFBAMatrix newFBAMatrixPureComp numIntersectMatrix numIntersectMatrixAdd numIntersectMatrixMets numIntersectMatrixMetsAdd] = simulateSmallModelsSeparateFactored(modelLimit,modelNames,modelNamesToModels,allBiomassRates,allBiomassDists)

FBAGapMatrix = [];
newFBAMatrix = [];
newFBAMatrixPureComp = [];
numIntersectMatrix = zeros(length(modelNames),length(modelNames));
numIntersectMatrixAdd = zeros(length(modelNames),length(modelNames));
numIntersectMatrixMets = {};
numIntersectMatrixMetsAdd = {};
for i=1:modelLimit%length(modelNames)
    for j=1:modelLimit%length(modelNames)
	numIntersectMatrixMets{i,j} = {};
        numIntersectMatrixMetsAdd{i,j} = {};
    end
end
for i=1:modelLimit%length(modelNames)
    for j=1:modelLimit%length(modelNames)
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
	%splitRatio = rand(1);
        splitRatio = allBiomassRates(i)/(allBiomassRates(i)+allBiomassRates(j));
	if isnan(splitRatio) || isinf(splitRatio)
	    splitRatio = 0.5;
	end
        for k=1:length(involvedMetKEGGs1)
            if length(involvedMetKEGGs1{k})==1
                intersectIdx = find(cellfun(@(x) length(x)==1 && strcmp(involvedMetKEGGs1{k}{1},x{1}), involvedMetKEGGs2));
                if length(intersectIdx)==1
                    %length(intersectIdx)
                    %disp([num2str(i) ' ' num2str(j) ' ' num2str(k)]);
                    totalFlux = dist1(involvedIdxs1(k))+dist2(involvedIdxs2(intersectIdx));
                    model1.lb(k) = totalFlux*splitRatio;
                    model2.lb(intersectIdx) = totalFlux*(1-splitRatio);
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
