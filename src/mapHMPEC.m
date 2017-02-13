configSEED;
inputDir2 = [inputDir filesep 'HMPRef'];
HMPRefFiles = dir(inputDir2);
count = 0;
allECs = {};
percentageGeneCovs = [];
percentageECCovs = [];
speciesNames = {};
allECLengths = [];
flatTableECNums = unique([bigModelTable.rxnECNums{:}]);
geneNameToEC = containers.Map;
HMPIDsToGeneNames = containers.Map;
for i=1:length(HMPRefFiles)
    if ~strcmp(HMPRefFiles(i).name,'.') && ~ strcmp(HMPRefFiles(i).name,'..') && count<1000
        disp(i)
        disp(HMPRefFiles(i).name)
        speciesName = HMPRefFiles(i).name;
        speciesName = speciesName(1:regexp(speciesName,'2')-1);
        speciesNames{end+1} = speciesName;
        count = count+1;
        HMPFI = fopen([inputDir2 filesep HMPRefFiles(i).name]);
        dataFields = textscan(HMPFI,'%s','Delimiter','\n', 'HeaderLines',0);
        dataFields = [dataFields{:}];
        HMPIDs = cellfun(@(x) x(2:regexp(x,'\|')-1),dataFields,'UniformOutput',0);
        geneNames = cellfun(@(x) x(regexp(x,'\|')+1:end),dataFields,'UniformOutput',0);
        for z=1:length(HMPIDs)
            if ~isKey(HMPIDsToGeneNames,HMPIDs{z})
                HMPIDsToGeneNames(HMPIDs{z}) = geneNames{z};
            end
        end
        geneNames = geneNames(cellfun(@(x) ~isempty(x),geneNames));
        %covered = cellfun(@(x) any(strcmp(x,strrep(flatDEs,' ','_'))),geneNames);
        geneNames = geneNames(cellfun(@(x) isempty(regexp(x,'uncharacterized|putative|hypothetical|predicted')), geneNames));
        origGeneNames = geneNames;
        percentageCov = 0;
        covered = cellfun(@(x) any(strcmp(x,flatDEs)),geneNames);
        coveredGeneNames = geneNames(covered);
        percentageCov = length(coveredGeneNames)/length(origGeneNames);
        percentageGeneCovs(end+1,1) = percentageCov;
        disp(percentageCov)
        
        coveredTIGRs = cellfun(@(x) DEsToTIGRIDs(x),coveredGeneNames,'UniformOutput',0);
        coveredTIGRs = [coveredTIGRs{:}];
        coveredECs = cellfun(@(x) TIGRIDsToECNums(x),coveredTIGRs,'UniformOutput',0);
        coveredECs = unique([coveredECs{:}]);
        allECs = union(allECs,coveredECs);
        percentageECCovs(end+1,1) = length(intersect(allECs,flatTableECNums))/length(allECs);
        disp(length(allECs))

        for z=1:length(coveredGeneNames)
            if ~isKey(geneNameToEC,coveredGeneNames{z})
                if isKey(DEsToTIGRIDs,coveredGeneNames{z})
                    candTIGR = DEsToTIGRIDs(coveredGeneNames{z});
                end
                if isKey(TIGRIDsToECNums,candTIGR)
                    candEC = TIGRIDsToECNums(candTIGR{1});
                end
                geneNameToEC(coveredGeneNames{z}) = candEC;
            end
        end
        disp(length(keys(geneNameToEC)))        
        
        geneNames = setdiff(geneNames,coveredGeneNames);
        
        covered = cellfun(@(x) any(strcmp(x,flatPfamLong)),geneNames);
        coveredGeneNames = geneNames(covered);
        percentageCov = percentageCov+length(coveredGeneNames)/length(origGeneNames);
        percentageGeneCovs(end,2) = percentageCov;
        disp(percentageCov)
        coveredPfams = cellfun(@(x) pfamLongToPfam(x),coveredGeneNames,'UniformOutput',0);
        coveredPfams = [coveredPfams{:}];
        if isempty(coveredPfams)
            coveredPfams = {};
        end
        coveredECs = cellfun(@(x) pfamToEC(x),coveredPfams,'UniformOutput',0);
        coveredECs = unique([coveredECs{:}]);
        if isempty(coveredECs)
            coveredECs = {};
        end
        allECs = union(allECs,coveredECs);
        percentageECCovs(end,2) = length(intersect(allECs,flatTableECNums))/length(allECs);
        disp(length(allECs))

        for z=1:length(coveredGeneNames)
            if ~isKey(geneNameToEC,coveredGeneNames{z})
                if isKey(pfamLongToPfam,coveredGeneNames{z})
                    candPfam = pfamLongToPfam(coveredGeneNames{z});
                end
                if isKey(pfamToEC,candPfam)
                    candEC = pfamToEC(candPfam{1});
                end
                geneNameToEC(coveredGeneNames{z}) = candEC;
            end
        end
        disp(length(keys(geneNameToEC)))

        geneNames = setdiff(geneNames,coveredGeneNames);

        allECLengths(end+1) = length(allECs);
        %break;
        fclose(HMPFI);
    end
end