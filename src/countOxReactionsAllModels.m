AGORAFiles = dir([inputDir filesep 'AGORAModels/Western-Diet-Paper']);
AGORAMat = {};
count = 0;
for i = 1:length(AGORAFiles)
    dotIdx = regexp(AGORAFiles(i).name,'\.');
    if ~strcmp(AGORAFiles(i).name,'.') && ~strcmp(AGORAFiles(i).name,'..') && strcmp(AGORAFiles(i).name(dotIdx+1:end),'mat')
        AGORAMat{end+1} = [inputDir filesep 'AGORAModels/Western-Diet-Paper' filesep AGORAFiles(i).name];
    end
end

if 1
AGORAModelsArr = {};
for i=1:length(AGORAMat)
    load(AGORAMat{i});
i
    AGORAModelsArr{i} = AGORAModel;
end
end
oxFermNums = [];
for i=1:length(AGORAModelsArr)
	i
    oxFermNums(i) = measureOxFermFunc(AGORAModelsArr{i});
end
