if ~exist('AGORAMat','var')
AGORAFiles = dir([inputDir filesep 'AGORAModels/Western-Diet-Paper']);
AGORAMat = {};
count = 0;
for i = 1:length(AGORAFiles)
    dotIdx = regexp(AGORAFiles(i).name,'\.');
    if ~strcmp(AGORAFiles(i).name,'.') && ~strcmp(AGORAFiles(i).name,'..') && strcmp(AGORAFiles(i).name(dotIdx+1:end),'mat')
        AGORAMat{end+1} = [inputDir filesep 'AGORAModels/Western-Diet-Paper' filesep AGORAFiles(i).name];
    end
end
end
