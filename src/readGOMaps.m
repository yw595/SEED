configSEED;
outputDir1 = [outputDir filesep 'readGOMaps'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1);
end

inputFI1 = fopen([inputDir filesep 'ShoaieRefInfo' filesep 'GOmappings' filesep 'pfam2go']);
dataFields = textscan(inputFI1,'%s%s','Delimiter',{' > '}, 'HeaderLines',6);
fclose(inputFI1);
dataFields = [dataFields{:}];
%dataFields(:,1) = cellfun(@(x) x(1:end-1),dataFields(:,1),'UniformOutput',0);
%dataFields(:,2) = cellfun(@(x) x(1:end),dataFields(:,2),'UniformOutput',0); 
pfams = cellfun(@(x) strsplitYiping(x,' '), dataFields(:,1),'UniformOutput',0); pfams= cellfun(@(x) x{1}, pfams, 'UniformOutput',0); 
pfams = cellfun(@(x) strsplitYiping(x,':'), pfams,'UniformOutput',0); pfams= cellfun(@(x) x{2}, pfams, 'UniformOutput',0); 
goTerms = cellfun(@(x) strsplitYiping(x,' ; '),dataFields(:,2),'UniformOutput',0);
pfamToGo = containers.Map;
for i=1:length(pfams)
    if isKey(pfamToGo,pfams{i})
        currentGO = pfamToGo(pfams{i});
        currentGO = unique([currentGO goTerms{i}]);
        pfamToGo(pfams{i}) = currentGO;
    else
        pfamToGo(pfams{i}) = goTerms{i};
    end
end

inputFI2 = fopen([inputDir filesep 'ShoaieRefInfo' filesep 'GOmappings' filesep 'ec2go']);
dataFields = textscan(inputFI2,'%s%s','Delimiter','>', 'HeaderLines',2);
fclose(inputFI2);
dataFields = [dataFields{:}];
dataFields(:,1) = cellfun(@(x) x(1:end-1),dataFields(:,1), 'UniformOutput',0);
dataFields(:,2) = cellfun(@(x) x(1:end),dataFields(:,2), 'UniformOutput',0);
ecs = dataFields(:,1); ecs = cellfun(@(x) x(4:end),ecs,'UniformOutput',0);
for i=1:length(ecs)
    %use strsplit, my version gives error on both . and \.
    words = strsplit(ecs{i},'.');
    for j=length(words)+1:4
        words{end+1} = '-';
    end
    ecs{i} = strjoin(words,'.');
end
goTerms = cellfun(@(x) strsplitYiping(x,' ; '),dataFields(:,2), 'UniformOutput',0);
goToEC = containers.Map;
for i=1:length(goTerms)
    ithGoTerms = goTerms{i};
    for j=1:length(ithGoTerms)
        if isKey(goToEC,ithGoTerms{j})
            currentEC = goToEC(ithGoTerms{j});
            currentEC = unique([currentEC {ecs{i}}]);
            goToEC(ithGoTerms{j}) = currentEC;
        else
            goToEC(ithGoTerms{j}) = {ecs{i}};
        end
    end
end

pfamToEC = containers.Map;
for i=1:length(pfams)
    currentGO = pfamToGo(pfams{i});
    for j=1:length(currentGO)
        if isKey(goToEC,currentGO{j})
            if isKey(pfamToEC,pfams{i})
                currentEC = pfamToEC(pfams{i});
                currentEC = unique([currentEC goToEC(currentGO{j})]);
                pfamToEC(pfams{i}) = currentEC;
            else
                pfamToEC(pfams{i}) = goToEC(currentGO{j});
            end
        end
    end
end

inputFI3 = fopen([inputDir filesep 'ShoaieRefInfo' filesep 'Pfam-A.clans.tsv']);
dataFields = textscan(inputFI3,'%s%s%s%s%s','Delimiter','\t','HeaderLines',0);
dataFields = [dataFields{:}];
pfamToPfamLong = containers.Map; pfamLongToPfam = containers.Map;
pfamToECKeys = keys(pfamToEC);
for i=1:length(pfamToECKeys)
    dataFieldIdx = find(strcmp(dataFields(:,1),pfamToECKeys{i}));
    if ~isempty(dataFieldIdx)
        pfamLongName = dataFields{dataFieldIdx,5};
        pfamLongName = strrep(pfamLongName,' ','_');
        updateTwoMaps(pfamToPfamLong,pfamLongToPfam,pfamToECKeys{i},pfamLongName);
    end
end

flatPfamLong = cellfun(@(x) x{1}, values(pfamToPfamLong), 'UniformOutput',0);

save([outputDir1 filesep 'readGOMaps.mat'],'pfamToEC','pfamToPfamLong','pfamLongToPfam','flatPfamLong');