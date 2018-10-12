if 0

excludeModels = 0;
%modelsToExclude = lessThanTwentyModels;

ucrTotalArr = [];
percentageLessThanTwentyArr = [];
fbaTotalArr = [];
subsystemExpressionArr = {};
parfor z=1:length(ucrFolders)
    if excludeModels
        [ucrTotal, percentageLessThanTwenty, fbaTotal] = countOxPresentFunc(z,ucrFolders,inputDir,AGORAMat,fbaOxMat,modelsToExclude);
        ucrTotalArr(z) = ucrTotal;
        percentageLessThanTwentyArr(z) = percentageLessThanTwenty;
        fbaTotalArr(z) = fbaTotal;
    else
        [ucrTotal, percentageLessThanTwenty, fbaTotal, subsystemExpression] = countOxPresentFunc(z,ucrFolders,inputDir,AGORAMat,fbaOxMat);
        ucrTotalArr(z) = ucrTotal;
        percentageLessThanTwentyArr(z) = percentageLessThanTwenty;
        fbaTotalArr(z) = fbaTotal;
        subsystemExpressionArr{z} = subsystemExpression;
    end
end
end

allSubsystems = {};
for i=1:length(subsystemExpressionArr)
    allSubsystems = union(allSubsystems,keys(subsystemExpressionArr{i}));
end
subsystemExpressionMat = zeros(length(ucrNames),length(allSubsystems));
for i=1:length(ucrNames)
    subsystemExpression = subsystemExpressionArr{i};
    for j=1:length(allSubsystems)
	if isKey(subsystemExpression,allSubsystems{j})
	    subsystemExpressionMat(i,j) = subsystemExpression(allSubsystems{j});
        end
    end
end
obeseids = importdata('/mnt/vdb/home/ubuntu2/MHobese.txt');
normalids = importdata('/mnt/vdb/home/ubuntu2/MHnormal.txt');
    
[~,~,intersectIdxs] = intersect(obeseids,ucrNames);
normalOx = ucrTotalArr(intersectIdxs);
normalPct = percentageLessThanTwentyArr(intersectIdxs);
normalFBA = fbaTotalArr(intersectIdxs);
normalSubsystemExpressionMat = subsystemExpressionMat(intersectIdxs,:);
[~,~,intersectIdxs] = intersect(normalids,ucrNames);
obeseOx = ucrTotalArr(intersectIdxs);
obesePct = percentageLessThanTwentyArr(intersectIdxs);
obeseFBA = fbaTotalArr(intersectIdxs);
obeseSubsystemExpressionMat = subsystemExpressionMat(intersectIdxs,:);
for i=1:length(allSubsystems)
    [h p] = ttest2(normalSubsystemExpressionMat(:,i),obeseSubsystemExpressionMat(:,i))
end
obeseOx = [obeseOx obeseOx];
normalOx = [normalOx normalOx];
obesePct = [obesePct obesePct];
normalPct = [normalPct normalPct];
obeseFBA = [obeseFBA obeseFBA];
normalFBA = [normalFBA normalFBA];
oxAmounts = [obeseOx normalOx];
pctAmounts = [obesePct normalPct];
fbaAmounts = [obeseFBA normalFBA];
grouplabels = {};
for i=1:length(obeseOx);
    grouplabels{end+1} = 'obexe';
end
for i=1:length(normalOx)
    grouplabels{end+1} = 'normal';
end
writeData({oxAmounts,grouplabels},['/mnt/vdb/home/ubuntu2' filesep 'oxAmounts16s.txt'],'\t',{'oxamount','group'});
writeData({oxAmounts,grouplabels},['/mnt/vdb/home/ubuntu2' filesep 'oxAmounts16sLessThanTwenty.txt'],'\t',{'oxamount','group'});
writeData({pctAmounts,grouplabels},['/mnt/vdb/home/ubuntu2' filesep 'pcts16sLessThanTwenty.txt'],'\t',{'pct','group'});
writeData({1-pctAmounts,grouplabels},['/mnt/vdb/home/ubuntu2' filesep 'pcts16sMoreThanTwenty.txt'],'\t',{'pct','group'});
writeData({fbaAmounts,grouplabels},['/mnt/vdb/home/ubuntu2' filesep 'fbaOxAmounts16s.txt'],'\t',{'oxamount','group'});
