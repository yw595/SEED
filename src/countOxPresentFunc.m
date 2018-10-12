function [ucrTotal percentageLessThanTwenty fbaTotal subsystemExpression] = countOxPresentFunc(z,ucrFolders,inputDir,AGORAMat,fbaOxMat,modelsToExclude)
        ucrTotal = 0;
        fbaTotal = 0;
        subsystemExpression = containers.Map;
        ucrFolder = ucrFolders{z};
        system(['python /mnt/vdb/home/ubuntu2/pickPresent.py ' inputDir filesep 'MGMData' filesep ucrFolder filesep 'normalized_otus.tsv ' num2str(z)]);
        presentSpecies = textscan(fopen(['/mnt/vdb/home/ubuntu2/pickPresent' num2str(z) '.txt']),'%s','Delimiter','\n','HeaderLines',0);
	presentSpecies = presentSpecies{1};
	bigPresent = '';
        totalAll = 0;
        totalLessThanTwenty = 0;
	for i=1:length(presentSpecies)-1
	    bigPresent = [bigPresent presentSpecies{i} ','];
	end
	bigPresent = [bigPresent presentSpecies{end}];
	for i=1:length(presentSpecies)
	    [status, result] = system(['python /mnt/vdb/home/ubuntu2/MATLAB/SEED/src/simulateSmallModelsSeparatePicrustAGORA2.py ' ucrFolder ' ' presentSpecies{i}]);
	    result = strsplit(result,'Result ');
            result = strsplit(result{2},' ');
            totalresult = str2num(result{2});
            result = str2num(result{1});
            disp(result/totalresult)
	    for j=1:length(AGORAMat)
		if ~isempty(regexp(AGORAMat{j},presentSpecies{i}))
		    load(AGORAMat{j})
		    oxFermNum = measureOxFermFunc(AGORAModel);
                    totalAll = totalAll+result;
                    if oxFermNum <= 20
                        totalLessThanTwenty = totalLessThanTwenty+result/totalresult;
                    end
		    ucrTotal = ucrTotal + oxFermNum*result/totalresult;
                    fbaTotal = fbaTotal + fbaOxMat(j,2)*result/totalresult;

                    uniqSubsystems = unique(cellfun(@(x) x{:},AGORAModel.subSystems,'UniformOutput',0));
                    AGORAModel.subSystems = cellfun(@(x) x{:},AGORAModel.subSystems,'UniformOutput',0);
                    for k=1:length(uniqSubsystems)
			disp(sum(strcmp(AGORAModel.subSystems,uniqSubsystems{k})))    
	                if isKey(subsystemExpression,uniqSubsystems{k})
	                    subsystemExpression(uniqSubsystems{k}) = subsystemExpression(uniqSubsystems{k}) + sum(strcmp(AGORAModel.subSystems,uniqSubsystems{k}))*result/totalresult;
	                else
	                    subsystemExpression(uniqSubsystems{k}) = sum(strcmp(AGORAModel.subSystems,uniqSubsystems{k}))*result/totalresult;
	                end
		    end
		end
	    end
	end
        percentageLessThanTwenty = totalLessThanTwenty/totalAll;
	%ucrTotalArr(z) = ucrTotal;
