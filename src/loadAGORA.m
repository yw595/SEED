AGORAFiles = dir([inputDir filesep 'AGORAModels/Western-Diet-Paper']);
totalCount = 0;
while totalCount*20 < length(AGORAFiles)
    AGORAModelsArr = {};
    parfor i=totalCount*20+1:min((totalCount+1)*20,length(AGORAFiles))
        dotIdx = regexp(AGORAFiles(i).name,'\.');
        baseFilename = AGORAFiles(i).name(1:dotIdx-1);
        if ~strcmp(AGORAFiles(i).name,'.') && ~strcmp(AGORAFiles(i).name,'..') && strcmp(AGORAFiles(i).name(dotIdx+1:end),'xml') && ~exist([inputDir filesep 'AGORAModels/Western-Diet-Paper' filesep baseFilename '.mat'],'file')
	disp(baseFilename)
	    AGORAModel = readCbModel([inputDir filesep 'AGORAModels/Western-Diet-Paper' filesep AGORAFiles(i).name]);
            AGORADesc = AGORAModel.description;
	    AGORADesc = AGORADesc(1:length(AGORADesc)-4);
	    AGORAModelsArr{i} = AGORAModel;
	else
	    AGORAModelsArr{i} = '';
	end
	disp(i);
    end

    for i=1:length(AGORAModelsArr)
	if ~isempty(AGORAModelsArr{i})
	    AGORAModel = AGORAModelsArr{i};
	    AGORADesc = AGORAModel.description;
	    AGORADesc = AGORADesc(1:length(AGORADesc)-4);
	    save([inputDir filesep 'AGORAModels/Western-Diet-Paper' filesep AGORADesc '.mat'],'AGORAModel');
        end
    end
    totalCount = totalCount+1;
end
