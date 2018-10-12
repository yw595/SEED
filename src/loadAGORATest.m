AGORAFiles = dir([inputDir filesep 'AGORAModels/Western-Diet-Paper']);
count = 0;
AGORAModelArr = {};
for i = 1:length(AGORAFiles)
    dotIdx = regexp(AGORAFiles(i).name,'\.');
    if ~strcmp(AGORAFiles(i).name,'.') && ~strcmp(AGORAFiles(i).name,'..') && strcmp(AGORAFiles(i).name(dotIdx+1:end),'mat')
        count = count+1;
        load([inputDir filesep 'AGORAModels/Western-Diet-Paper' filesep AGORAFiles(i).name]);
        disp(count);
        AGORAModelArr{end+1} = AGORAModel;
        disp(length(AGORAModel.rxns));
    end
end
