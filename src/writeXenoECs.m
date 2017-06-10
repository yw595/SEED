configSEED;
outputDir1 = [outputDir filesep 'writeXenoEC'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1);
end
inputDir1 = [inputDir filesep 'xenobiotics'];
xenoFI = fopen([inputDir1 filesep 'GSM935962_A1_EtOH_CDS.counts.txt']);
dataFields = textscan(xenoFI,'%s%s%s','Delimiter','\t', 'HeaderLines',0);
fclose(xenoFI);
dataFields = [dataFields{:}];
xenoHMPIDs = {};
xenoMeans = dataFields(:,2);
xenoStds = dataFields(:,3);
xenoMeans = cellfun(@(x) str2num(x), xenoMeans);
xenoStds = cellfun(@(x) str2num(x), xenoStds);
for i=1:length(dataFields)
    xenoHMPIDs{i} = dataFields{i,1};
    regIdx = regexp(xenoHMPIDs{i},'_');
    regIdx = regIdx(2);
    xenoHMPIDs{i} = xenoHMPIDs{i}(1:regIdx-1);
end
%xenoHMPIDs = cellfun(@(x) x(1:regexp(x,'_')(2)-1),dataFields(:,1),'UniformOutput',0);
count = 0;
count1 = 0;
xenoExpECs = {};
xenoExps = [];
xenoExpStds = [];
for i=1:length(xenoHMPIDs)
    disp(i)
    xenoHMPID = xenoHMPIDs{i};
    if isKey(HMPIDsToGeneNames,xenoHMPID)
        candGeneName = HMPIDsToGeneNames(xenoHMPID);
        count1 = count1+1;
        if isKey(geneNameToEC,candGeneName)
            candECs = geneNameToEC(candGeneName);
            for j=1:length(candECs)
                candEC = candECs{j};
                count = count+1;
                disp([num2str(count) ' ' num2str(count1)]);
                disp(length(xenoExpECs))
                existIdx = find(strcmp(xenoExpECs,candEC));
                if ~isempty(existIdx)
                    xenoExps(existIdx) = xenoExps(existIdx) + xenoMeans(i);
                    temp = xenoExpStds(existIdx);
                    temp = temp^2;
                    temp = temp + xenoStds(i)^2;
                    xenoExpStds(existIdx) = sqrt(temp);
                else
                    xenoExpECs{end+1} = candEC;
                    xenoExps(end+1) = xenoMeans(i);
                    xenoExpStds(end+1) = xenoStds(i);
                end
            end
        end
    end
end

save([outputDir1 filesep 'writeXenoECs.mat'],'xenoExpECs','xenoExps','xenoExpStds');

writeData({xenoExpECs,xenoExps,xenoExpStds},[outputDir1 filesep 'GSM935962_A1_EtOH_CDS.EC.txt'],'\t');