configSEED;
outputDir1 = [outputDir filesep 'readTIGRFams'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1);
end

inputDir1 = [inputDir filesep 'TIGRFAMS'];
TIGRFiles = dir(inputDir1);
TIGRIDsToECNums = containers.Map;
TIGRIDsToDEs = containers.Map;
DEsToTIGRIDs = containers.Map;
for i=1:length(TIGRFiles)
    if ~strcmp(TIGRFiles(i).name,'.') && ~strcmp(TIGRFiles(i).name,'..')
    if mod(i,100)==0
        disp(i);
    end
    TIGRFI = fopen([inputDir1 filesep TIGRFiles(i).name]);
    dataFields = textscan(TIGRFI,'%s','Delimiter','\n', 'HeaderLines',0);
    dataFields = [dataFields{:}];
    TIGRIdx = regexp(dataFields,'^AC'); TIGRIdx = find(cellfun(@(x) ~isempty(x),TIGRIdx));
    TIGRID = dataFields{TIGRIdx}; TIGRID = TIGRID(5:end);
    ECIdx = regexp(dataFields,'^EC'); ECIdx = find(cellfun(@(x) ~isempty(x),ECIdx));
    if ~isempty(ECIdx)
        ECNum = dataFields{ECIdx}; ECNum = strsplitYiping(ECNum(5:end),' ');
        TIGRIDsToECNums(TIGRID) = ECNum;
        DEIdx = regexp(dataFields,'^DE'); DEIdx = find(cellfun(@(x) ~isempty(x),DEIdx));
        DE = strrep(dataFields{DEIdx},' ','_'); DE = DE(5:end);
        updateTwoMaps(TIGRIDsToDEs,DEsToTIGRIDs,TIGRID,DE);
        %TIGRIDsToDEs(TIGRID) = DE;
        %DEsToTIGRIDs(DE) = TIGRID;
    end
    fclose(TIGRFI);
    end
end
flatDEs = cellfun(@(x) x{1}, values(TIGRIDsToDEs), 'UniformOutput',0);

save([outputDir1 filesep 'readTIGRFams.mat'],'TIGRIDsToECNums','TIGRIDsToDEs','DEsToTIGRIDs');