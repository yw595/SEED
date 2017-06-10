configSEED;
outputDir1 = [outputDir filesep 'mergeSmallModels'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1);
end
allTwoModels = {};
allTwoBiomasses = [];
%parfor z=1:length(allModelKeys)^2
%i = floor(z/length(allModelKeys))+1;
%j = mod(z,length(allModelKeys))+1;
allTwoBiomasses = [];
for i=1:length(modelNames)
    for j=1:length(modelNames)
        if i~=j
            i
            j
            allTwoBiomasses(i,j) = mergeSmallModelsFactored(i,j,modelNamesToModels,cpdKEGGs,cpdIDs);
        end
    end
end

if 0
    while length(allTwoBiomassesArr)<length(modelNames)*length(modelNames)
        allTwoBiomassesArr(end+1) = 0;%rand(1)*1000;
    end

    zCount = 1;
    for i=1:length(modelNames)
        for j=1:length(modelNames)
            allTwoBiomasses(i,j) = allTwoBiomassesArr(zCount);
            zCount = zCount+1;
        end
    end
end
save([outputDir1 filesep 'mergeSmallModels.mat'],'allTwoBiomasses');





