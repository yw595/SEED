configSEED;

extraBiomMat = [];

for i=1:length(modelNames)
    for j=1:length(modelNames)
        if i~=j
            minBiom = min([allBiomassRates(i),allBiomassRates(j)]);
            extraBiomMat(i,j) = allTwoBiomasses(i,j) - minBiom;%(allBiomassRates(i)+allBiomassRates(j))/2;
            if extraBiomMat(i,j)<0
                extraBiomMat(i,j) = mergeSmallModelsFactored(i,j,modelNamesToModels,cpdKEGGs,cpdIDs) - minBiom;
            end
        end
    end
end

avgExtraBioms = [];
for i=1:length(modelNames)
    avgExtraBioms(i) = mean(extraBiomMat(i,:));
end
[avgExtraBioms, sortIdxs] = sort(avgExtraBioms,'descend');
sortModelNames = modelNamesShort(sortIdxs);
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
            modelArr1{end+1} = modelNamesShort{i};
            modelArr2{end+1} = modelNamesShort{j};
            extraBiomArr(end+1) = extraBiomMat(i,j);
        end
    end
end
[extraBiomArrSorted sortIdxsExtra] = sort(extraBiomArr,'descend');
%modelArr1 = modelArr1(sortIdxsExtra);
%modelArr2 = modelArr2(sortIdxsExtra);
writeData({modelArr1,modelArr2,extraBiomArr},'/home/fs01/yw595/coopBiomassHeatmap.txt','\t',{'species1','species2','extrabiom'});








