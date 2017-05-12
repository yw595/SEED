inputFI1 = fopen('/home/fs01/yw595/ko2ec.xl.txt');
dataFields = textscan(inputFI1,'%s%s','Delimiter','\t','HeaderLines',1);
fclose(inputFI1);
dataFields = [dataFields{:}];
ko2ec = containers.Map;
for i=1:length(dataFields)
    ko = dataFields{i,1};
    ecs = strsplit(dataFields{i,2}(5:end-1),' ');
    disp(ko)
    ko2ec(ko) = ecs;
end

useHadza = 0; useMatsumoto = 1;

if useHadza
    
inputFI2 = fopen('/home/fs01/yw595/manualmetmap.txt');
dataFields = textscan(inputFI2,'%s%s','Delimiter','|','HeaderLines',0);
fclose(inputFI2);
dataFields = [dataFields{:}];
metsToIDs = containers.Map;
for i=1:length(dataFields)
    metsToIDs(dataFields{i,1}) = dataFields{i,2};
end

inputFI3 = fopen('/home/fs01/yw595/hadzametabolomesupp.txt');
dataFields = textscan(inputFI3,'%s%s%s%s','Delimiter','|','HeaderLines',0);
fclose(inputFI3);
dataFields = [dataFields{:}];
metsToAbunds = containers.Map;
for i=1:length(dataFields)
    metsToAbunds(dataFields{i,1}) = [str2num(dataFields{i,3}), str2num(dataFields{i,4})];
end
end

if useMatsumoto
    inputFI3 = fopen('/home/fs01/yw595/matsumotometabolomesupp.txt');
    dataFields = textscan(inputFI3,'%s%s%s%s','Delimiter','|','HeaderLines',0);
    fclose(inputFI3);
    dataFields = [dataFields{:}];
    metsToAbunds = containers.Map;
    metsToIDs = containers.Map;
    for i=1:length(dataFields)
        metsToIDs(dataFields{i,1}) = dataFields{i,4};
        metsToAbunds(dataFields{i,1}) = [str2num(dataFields{i,2}), str2num(dataFields{i,3})];
    end
end

ecs = {}; ecAbunds = [];
inputFI4 = fopen('/home/fs01/yw595/output2.res');
dataFields = textscan(inputFI4,'%s%s%s%s%s%s%s','Delimiter','\t','HeaderLines',0);
fclose(inputFI4);
dataFields = [dataFields{:}];
for i=1:length(dataFields)
    ko = dataFields{i,4};
    if isKey(ko2ec,ko)
        matchECs = ko2ec(ko);
        for j=1:length(matchECs)
            alreadyIn = find(strcmp(matchECs{j},ecs));
            if ~isempty(alreadyIn)
                ecAbunds(alreadyIn) = ecAbunds(alreadyIn)+1;
            else
                ecs{end+1} = matchECs{j};
                ecAbunds(end+1) = 1;
            end
        end
    end
end
writeData({ecs,ecAbunds},'/home/fs01/yw595/output2.EC.txt','\t',{'ec','ecAbund'});

ecs = {}; ecAbunds = [];
inputFI4 = fopen('/home/fs01/yw595/output3.res');
dataFields = textscan(inputFI4,'%s%s%s%s%s%s%s','Delimiter','\t','HeaderLines',0);
fclose(inputFI4);
dataFields = [dataFields{:}];
for i=1:length(dataFields)
    ko = dataFields{i,4};
    if isKey(ko2ec,ko)
        matchECs = ko2ec(ko);
        for j=1:length(matchECs)
            alreadyIn = find(strcmp(matchECs{j},ecs));
            if ~isempty(alreadyIn)
                ecAbunds(alreadyIn) = ecAbunds(alreadyIn)+1;
            else
                ecs{end+1} = matchECs{j};
                ecAbunds(end+1) = 1;
            end
        end
    end
end
writeData({ecs,ecAbunds},'/home/fs01/yw595/output3.EC.txt','\t',{'ec','ecAbund'});






















    






