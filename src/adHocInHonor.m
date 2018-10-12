if 0
%clear all
configSEED;
end

load([outputDir filesep 'simulateSmallModelsSeparatePicrustAnalysis' filesep 'simulateSmallModelsSeparatePicrustAnalysis.mat']);
load([outputDir filesep 'simulateSmallModelsSeparate' filesep 'simulateSmallModelsSeparate.mat']);
load('/mnt/vdb/home/ubuntu2/oxsecreteMat.mat')

pseudoNumSpecies = 643;
scrambleIdxs = randi(length(modelNamesShort),pseudoNumSpecies,1);
modelNamesShort = modelNamesShort(scrambleIdxs);
pseudoBiomassArrNormal = pseudoBiomassArrNormal(scrambleIdxs);
pseudoBiomassArrObese = pseudoBiomassArrObese(scrambleIdxs);
totalOxArr1 = totalOxArr1(scrambleIdxs);
totalOxArr2 = totalOxArr2(scrambleIdxs);
totalFermentArr1 = totalFermentArr1(scrambleIdxs);
totalFermentArr2 = totalFermentArr2(scrambleIdxs);
scrambleIdxsMat = [];
for i=1:length(scrambleIdxs)
    for j=1:length(scrambleIdxs)
	zPrev = scrambleIdxs(i)+(scrambleIdxs(j)-1)*max(scrambleIdxs);
	scrambleIdxsMat(end+1) = zPrev;
    end
end
pseudoBiomassMatNormal = reshape(pseudoBiomassMatNormal(scrambleIdxsMat),pseudoNumSpecies,pseudoNumSpecies);
pseudoBiomassMatObese = reshape(pseudoBiomassMatObese(scrambleIdxsMat),pseudoNumSpecies,pseudoNumSpecies);
newFBAMatrix = reshape(newFBAMatrix(scrambleIdxsMat),pseudoNumSpecies,pseudoNumSpecies);
newFBAMatrixPureComp = reshape(newFBAMatrixPureComp(scrambleIdxsMat),pseudoNumSpecies,pseudoNumSpecies);
oxsecreteMat2 = zeros(pseudoNumSpecies,pseudoNumSpecies,size(oxsecreteMat,3),size(oxsecreteMat,4));
for i=1:size(oxsecreteMat,3)
    for j=1:size(oxsecreteMat,4)
	oxsecreteMatTemp = oxsecreteMat(:,:,i,j);
        oxsecreteMatTemp = reshape(oxsecreteMatTemp(scrambleIdxsMat),pseudoNumSpecies,pseudoNumSpecies);
        oxsecreteMat2(:,:,i,j) = oxsecreteMatTemp;
    end
end
oxsecreteMat = oxsecreteMat2;

species1Arr = {};
species2Arr = {};
oxsecreteArr = [];
for i=1:length(modelNamesShort)
    for j=1:length(modelNamesShort)
	species1Arr{end+1} = modelNamesShort{i};
        species2Arr{end+1} = modelNamesShort{j};
        oxsecreteArr(end+1) = oxsecreteMat(i,j);
    end
end
writeData({species1Arr,species2Arr,oxsecreteArr},[transferDir filesep 'oxsecreteMat.txt'],'\t',{'species1','species2','oxsecrete'});

species1Arr = {};
species2Arr = {};
pairwisePseudoBiomassArr = [];
for i=1:length(modelNamesShort)
    for j=1:length(modelNamesShort)
	species1Arr{end+1} = modelNamesShort{i};
        species2Arr{end+1} = modelNamesShort{j};
        pairwisePseudoBiomassArr(end+1) = pseudoBiomassMatNormal(i,j);
    end
end
writeData({species1Arr,species2Arr,pairwisePseudoBiomassArr},[transferDir filesep 'pseudoBiomassMatNormal.txt'],'\t',{'species1','species2','pairwisebiom'})

dependencyMatNormal = [];
dependencyMatObese = [];
for i=1:length(pseudoBiomassMatNormal)
    for j=1:length(pseudoBiomassMatObese)
	dependencyMatNormal(i,j) = pseudoBiomassMatNormal(i,j)-(pseudoBiomassArrNormal(i)+pseudoBiomassArrNormal(j))/2;
	dependencyMatObese(i,j) = pseudoBiomassMatObese(i,j)-(pseudoBiomassArrObese(i)+pseudoBiomassArrObese(j))/2;
    end
end

useCompDependency = 1;
if useCompDependency
    dependencyMatNormal(isnan(dependencyMatNormal)) = 0;
    dependencyMatNormal(dependencyMatNormal > 0) = 0;%-dependencyMatNormal(dependencyMatNormal < 0);
    dependencyMatNormal(dependencyMatNormal < 0) = -dependencyMatNormal(dependencyMatNormal < 0);
    dependencyMatObese(isnan(dependencyMatObese)) = 0;
    dependencyMatObese(dependencyMatObese > 0) = 0;%-dependencyMatObese(dependencyMatObese < 0);
    dependencyMatObese(dependencyMatObese < 0) = -dependencyMatObese(dependencyMatObese < 0);
else
    dependencyMatNormal(isnan(dependencyMatNormal)) = 0;
    dependencyMatNormal(dependencyMatNormal < 0) = 0;%-dependencyMatNormal(dependencyMatNormal < 0);
    dependencyMatObese(isnan(dependencyMatObese)) = 0;
    dependencyMatObese(dependencyMatObese < 0) = 0;%-dependencyMatObese(dependencyMatObese < 0);
end

dependencyCentsNormal = betweenness_centrality(sparse(dependencyMatNormal)*1.0);
dependencyCentsObese = betweenness_centrality(sparse(dependencyMatObese)*1.0);
dependencyCentsNormoal = dependencyCentsNormal./sum(dependencyCentsNormal);
dependencyCentsObese = dependencyCentsObese./sum(dependencyCentsObese);
[dependencyCentsNormal, sortIdxs] = sort(dependencyCentsNormal,'descend');
dependencyCentsObese = dependencyCentsObese(sortIdxs);
doubledXCoords = ones(length(dependencyCentsNormal)*2,1);
doubledXCoords(1:2:length(doubledXCoords)) = 1:length(dependencyCentsNormal);
doubledXCoords(2:2:length(doubledXCoords)) = 1:length(dependencyCentsObese);
groupArr = {};
for i=1:length(doubledXCoords)
    if mod(i,2)==1
	groupArr{end+1} = 'Normal';
    else
        groupArr{end+1} = 'Obese';
    end
end
if useCompDependency
    writeData({doubledXCoords,[dependencyCentsNormal; dependencyCentsObese],groupArr},[transferDir filesep 'dependencyCentsComp.txt'],'\t',{'modelid','cent','group'});
else
    writeData({doubledXCoords,[dependencyCentsNormal; dependencyCentsObese],groupArr},[transferDir filesep 'dependencyCents.txt'],'\t',{'modelid','cent','group'});
end

if useCompDependency
    dependencyMatFBA = newFBAMatrixPureComp - newFBAMatrix;
    %dependencyMatFBA = (dependencyMatFBA(:,:,1)+dependencyMatFBA(:,:,2))/2;
    dependencyMatFBA(isnan(dependencyMatFBA)) = 0;
    dependencyMatFBA(dependencyMatFBA > 0) = 0;%-dependencyMatFBA(dependencyMatFBA < 0);
    dependencyMatFBA(dependencyMatFBA < 0) = -dependencyMatFBA(dependencyMatFBA < 0);
else
    dependencyMatFBA = newFBAMatrixPureComp - newFBAMatrix;
    %dependencyMatFBA = (dependencyMatFBA(:,:,1)+dependencyMatFBA(:,:,2))/2;
    dependencyMatFBA(isnan(dependencyMatFBA)) = 0;
    dependencyMatFBA(dependencyMatFBA < 0) = 0;%-dependencyMatNormal(dependencyMatFBA < 0);
end
dependencyCentsFBA = betweenness_centrality(sparse(dependencyMatFBA)*1.0);
dependencyCentsFBA = dependencyCentsFBA./sum(dependencyCentsFBA);
dependencyCentsFBA = sort(dependencyCentsFBA,'descend');
if useCompDependency
    writeData({1:length(dependencyCentsFBA),dependencyCentsFBA},[transferDir filesep 'dependencyCentsFBAComp.txt'],'\t',{'modelid','centfba'});
else
    writeData({1:length(dependencyCentsFBA),dependencyCentsFBA},[transferDir filesep 'dependencyCentsFBA.txt'],'\t',{'modelid','centfba'});
end

totalOxFermRatio1 = totalOxArr1./totalFermentArr1;
totalOxFermRatio1(isnan(totalOxFermRatio1) | isinf(totalOxFermRatio1)) = 0;
totalOxFermRatio2 = totalOxArr2./totalFermentArr2;
totalOxFermRatio2(isnan(totalOxFermRatio2) | isinf(totalOxFermRatio2)) = 0;
% totalOxFermRatio1Temp = [];
% totalOxFermRatio2Temp = [];
% for i=1:100
%     totalOxFermRatio1Temp(i) = sum(totalOxFermRatio1(randperm(length(totalOxFermRatio1),floor(0.9*length(totalOxFermRatio1)))));
% end
% for i=1:100
%     totalOxFermRatio2Temp(i) = sum(totalOxFermRatio2(randperm(length(totalOxFermRatio2),floor(0.9*length(totalOxFermRatio2)))));
% end
% totalOxFermRatio1 = totalOxFermRatio1Temp;
% totalOxFermRatio2 = totalOxFermRatio2Temp;
groupArr = {};
for i=1:length(totalOxFermRatio1)
    groupArr{end+1} = 'normal';
end
for i=1:length(totalOxFermRatio2)
    groupArr{end+1} = 'obese';
end
writeData({groupArr,[totalOxFermRatio1 totalOxFermRatio2]},[transferDir filesep 'oxFermRatios.txt'],'\t',{'group','ratio'});
