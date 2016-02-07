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
    disp(i)
    TIGRFI = fopen([inputDir1 filesep TIGRFiles(i).name]);
    dataFields = textscan(TIGRFI,'%s','Delimiter','\n', 'HeaderLines',0);
    dataFields = [dataFields{:}];
    TIGRIdx = regexp(dataFields,'^AC'); TIGRIdx = find(cellfun(@(x) ~isempty(x),TIGRIdx));
    TIGRID = dataFields{TIGRIdx}; TIGRID = TIGRID(5:end);
    ECIdx = regexp(dataFields,'^EC'); ECIdx = find(cellfun(@(x) ~isempty(x),ECIdx));
    if ~isempty(ECIdx)
        ECNum = dataFields{ECIdx}; ECNum = strsplit(ECNum(5:end),' ');
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

allECs = {};
inputDir2 = [inputDir filesep 'HMPRef'];
HMPRefFiles = dir(inputDir2);
count = 0;
for i=1:length(HMPRefFiles)
    if ~strcmp(HMPRefFiles(i).name,'.') && ~ strcmp(HMPRefFiles(i).name,'..')
        disp(i)
        disp(HMPRefFiles(i).name)
        count = count+1;
        HMPFI = fopen([inputDir2 filesep HMPRefFiles(i).name]);
        dataFields = textscan(HMPFI,'%s','Delimiter','\n', 'HeaderLines',0,'BufSize',50000);
        dataFields = [dataFields{:}];
        geneNames = cellfun(@(x) x(regexp(x,'\|')+1:end),dataFields,'UniformOutput',0);
        geneNames = geneNames(cellfun(@(x) ~isempty(x),geneNames));
        covered = cellfun(@(x) any(strcmp(x,strrep(flatDEs,' ','_'))),geneNames);
        geneNames = geneNames(covered);
        coveredTIGRs = cellfun(@(x) DEsToTIGRIDs(x),geneNames,'UniformOutput',0);
        coveredTIGRs = [coveredTIGRs{:}];
        coveredECs = cellfun(@(x) TIGRIDsToECNums(x),coveredTIGRs,'UniformOutput',0);
        coveredECs = unique([coveredECs{:}]);
        allECs = union(allECs,coveredECs);
        %break;
        fclose(HMPFI);
    end
end
save([outputDir1 filesep 'readTIGRFams'],'TIGRIDsToECNums','TIGRIDsToDEs','DEsToTIGRIDs','allECs');