function [totalSubsNames totalSubsRanks] = makeTotalSubsRanks(subLabelsArr)

totalSubsRanks = [];
totalSubsNames = {};
for i=1:length(subLabelsArr)
    subLabels = subLabelsArr{i};
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
