if 0
for i=1:length(AGORAMat)
	if ~isempty(AGORAMat{i})
	load(AGORAMat{i});
i
AGORAModelArr{i} = AGORAModel;
end
end
end

if 0
uniqueOxECs = containers.Map;
uniqueFermentECs = containers.Map;
for i=1:length(AGORAModelArr)
	[~,~,oxIdxs,fermentIdxs] = measureOxFermFunc(AGORAModelArr{i});
i
for j=1:length(oxIdxs)
	if ~isKey(uniqueOxECs,AGORAModelArr{i}.rxnECNumbers{oxIdxs(j)})
	uniqueOxECs(AGORAModelArr{i}.rxnECNumbers{oxIdxs(j)}) = '';
end
end
end
end

for i=1:length(AGORAModelArr)
	i
		[~,~,oxIdxs,fermentIdxs] = measureOxFermFunc(AGORAModelArr{i});
for j=1:length(fermentIdxs)
	if ~isKey(uniqueFermentECs,AGORAModelArr{i}.rxnECNumbers{fermentIdxs(j)})
	uniqueFermentECs(AGORAModelArr{i}.rxnECNumbers{fermentIdxs(j)}) = '';
end
end
end
