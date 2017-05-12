modelNames = keys(modelNamesToModels);
a = cellfun(@(x) strsplit(modelNamesToModels(x).modelName,'_'), modelNames, 'UniformOutput', 0);
b = {};
for i=1:length(a)
    temp = a{i};
    b{end+1} = [temp{1} ' ' temp{2}];
end

allTwoBiomasses = [];
while length(allTwoBiomassesArr)<length(allModelKeys)*length(allModelKeys)
    allTwoBiomassesArr(end+1) = 0;
end
zCount = 1;
for i=1:length(allModelKeys)
    for j=1:length(allModelKeys)
        allTwoBiomasses(i,j) = allTwoBiomassesArr(zCount);
        zCount = zCount+1;
    end
end

extraBiomMat = [];

for i=1:length(allModelKeys)
    for j=1:length(allModelKeys)
        if i~=j
            extraBiomMat(i,j) = allTwoBiomasses(i,j) - (allBiomassRates(i)+allBiomassRates(j))/2;
        end
    end
end

avgExtraBioms = [];
for i=1:length(allModelKeys)
    avgExtraBioms(i) = mean(extraBiomMat(i,:));
end
[avgExtraBioms, sortIdxs] = sort(avgExtraBioms,'descend');
sortModelNames = b(sortIdxs);
for i=1:length(sortModelNames)
    ithIdx = num2str(i);
    while length(ithIdx)<3
        ithIdx = ['0' ithIdx];
    end
    sortModelNames{i} = [ithIdx '_' sortModelNames{i}];
end
writeData({sortModelNames, avgExtraBioms},'/home/fs01/yw595/avgExtraBioms.txt','\t',{'species','avgextrabiom'});

modelArr1 = {}; modelArr2 = {}; extraBiomArr = [];
for i=1:length(allModelKeys)
    for j=1:length(allModelKeys)
        if i~=j
            modelArr1{end+1} = b{i};
            modelArr2{end+1} = b{j};
            extraBiomArr(end+1) = extraBiomMat(i,j);
        end
    end
end
[extraBiomArrSorted sortIdxsExtra] = sort(extraBiomArr,'descend');
writeData({modelArr1,modelArr2,extraBiomArr},'/home/fs01/yw595/coopBiomassHeatmap.txt','\t',{'species1','species2','extrabiom'});

allAbundsData = textscan(fopen('/home/fs01/yw595/MATLAB/SEED/src/HMPAllAbunds.txt'),repmat('%s',1,3840),'Delimiter','\t','HeaderLines',0);
closestFamiliesRevData = textscan(fopen('/home/fs01/yw595/closestFamiliesRev.txt'),'%s%s','Delimiter','|','HeaderLines',0);
nonsense = nonsense+1;

HMPFamilies = allAbundsData{1};
for i=1:length(HMPFamilies)
    matchIdxs1 = find(strcmp(closestFamiliesRevData{1},HMPFamilies{i}));
    for j=1:length(HMPFamilies)
        if i~=j
            matchIdxs2 = find(strcmp(closestFamiliesRevData{1},HMPFamilies{j}));
            if length(matchIdxs1) == 1 && length(matchIdxs2) == 1
                abunds1 = cellfun(@(x) str2num(x{i}),allAbundsData(2:end));
                abunds2 = cellfun(@(x) str2num(x{j}), allAbundsData(2:end));
                corrCoabundRhoArr(end+1) = corr(abunds1,abunds2,'type','Spearman');
                closestModel1 = closestFamiliesRevData{2}{matchIdxs1};
                closestModel2 = closestFamiliesRevData{2}{matchIdxs2};
                matchIdxs3 = find(strcmp(b,closestModel1));
                matchIdxs4 = find(strcmp(b,closestModel2));
                corrBiomArr(end+1) = extraBiomMat(matchIdxs3,matchIdxs4);
            end
        end
    end
end








