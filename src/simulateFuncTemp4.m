function [AGORAModelWithExpr] = simulateFunc(mergeModelsAGORAFlag,useFBA,singleSpeciesOnly,inputDir,outputDir1,ucrFolder,i,j1,blocksize,AGORAMat,AGORAModelArr,presentModelIdxs,realOneSpecies,doSimulation,justTwoSpecies,IBDFlag,justAGORAIBDMap)

if ~exist('IBDFlag','var')
    IBDFlag = 0
end
AGORAModelWithExpr = [];
disp(['i ' num2str(i)])
useCommon = 1;
alreadyRunPicrust = 0;
for j=(j1-1)*blocksize+1:min(length(AGORAMat),j1*blocksize)
    if i~=j && ~isempty(AGORAMat{i}) && ~isempty(AGORAMat{j}) && ~alreadyRunPicrust
	isPresent1 = 0;
	isPresent2 = 0;
        if presentModelIdxs(i)==1
	    isPresent1 = 1;
        end
	if presentModelIdxs(j)==1
	    isPresent2 = 1;
        end
	if 0%IBDFlag
	    AGORAName1 = strsplit(AGORAMat{i},'/');
	    AGORAName1 = strsplit(AGORAName1{end},'_');
	    AGORAName1 = [AGORAName1{1} ' ' AGORAName1{2}];
	    AGORAName2 = strsplit(AGORAMat{j},'/');
	    AGORAName2 = strsplit(AGORAName2{end},'_');
	    AGORAName2 = [AGORAName2{1} ' ' AGORAName2{2}];
	    if isKey(justAGORAIBDMap,AGORAName1) && isKey(justAGORAIBDMap,AGORAName2)
		for k=1:length(presentSpecies)
		    if ~isempty(regexp(presentSpecies{k},justAGORAIBDMap(AGORAName1)))
			isPresent1 = 1;
		    end
		    if ~isempty(regexp(presentSpecies{k},justAGORAIBDMap(AGORAName2)))
			isPresent2 = 1;
		    end
		end
	    end
	elseif 0
	    for k=1:length(presentSpecies)
		if ~isempty(regexp(AGORAMat{i},presentSpecies{k}))
		    isPresent1 = 1;
		end
		if ~isempty(regexp(AGORAMat{j},presentSpecies{k}))
		    isPresent2 = 1;
		end
	    end
	end

	if isPresent1==1 && isPresent2==1   
	    disp(['found ' num2str(i)])
	    disp(['found ' num2str(j)])
	    AGORAModel1 = AGORAModelArr{i};
	    AGORAModel2 = AGORAModelArr{j};
	    if mergeModelsAGORAFlag
		AGORAModelMerged = mergeModelsAGORA(AGORAModel1,AGORAModel2);
		status = 1;
		if useFBA==0
		    writeData({AGORAModelMerged.rxnECNumbers},[inputDir filesep 'MGMData/' ucrFolder '/speciesMerged.modelec'],'\t');
		    system(['python /mnt/vdb/home/ubuntu2/MATLAB/SEED/src/simulateSmallModelsSeparatePicrustAGORA4.py' ' --ucrFolder ' ucrFolder ' --modelnameslist ' AGORAModel1.description ',' AGORAModel2.description ' --modelnumberslist ' num2str(i) ',' num2str(j) ' --isSpeciesMerged True']);
		    status = importdata(['/mnt/vdb/home/ubuntu2/' ucrFolder AGORAModel1.description AGORAModel2.description 'True.txt']);
		end

		if status==0 || useFBA==1
		    simulateFuncFlux(useFBA,[inputDir filesep 'MGMData/' ucrFolder '/speciesMerged.modelexpr'],AGORAModelMerged,'_merged_');
		end

	    else
		status = 1;
		if useFBA==0 && realOneSpecies==0
		    writeData({AGORAModel1.rxnECNumbers},[inputDir filesep 'MGMData/' ucrFolder filesep 'speciesSep' num2str(i) '.modelec'],'\t');
		    writeData({AGORAModel2.rxnECNumbers},[inputDir filesep 'MGMData/' ucrFolder filesep 'speciesSep' num2str(j) '.modelec'],'\t');
		    if ~doSimulation
		    disp('where')
		        if useCommon
			    modelExprFile = [inputDir filesep 'MGMData' filesep 'IBDCommon' filesep 'speciesSep' num2str(i) '.modelexpr']
			else
		            modelExprFile = [inputDir filesep 'MGMData' filesep ucrFolder filesep 'speciesSep' num2str(i) '.modelexpr'];
                        end
                        if exist(modelExprFile,'file')
			    status = 0;
                            alreadyRunPicrust = 1;
			else
			    command = ['python /mnt/vdb/home/ubuntu2/MATLAB/SEED/src/simulateSmallModelsSeparatePicrustAGORA4.py' ' --ucrFolder ' ucrFolder ' --modelnameslist ' AGORAModel1.description ' --modelnumberslist ' num2str(i)]; 
			    if justTwoSpecies
			        command = [command ' --justTwoSpecies True'];
			    else
			        command = [command ' --justTwoSpecies False'];
			    end
			    if useCommon
			        command = [command ' --useCommon True'];
			    else
			        command = [command ' --useCommon False'];
			    end
			    system(command);
			    status = importdata(['/mnt/vdb/home/ubuntu2/' ucrFolder AGORAModel1.description 'False.txt']);
                            if status==0
                                alreadyRunPicrust = 1;
                            end
			end
		    else
			if IBDFlag
			    %system(['python /mnt/vdb/home/ubuntu2/MATLAB/SEED/src/simulateSmallModelsSeparatePicrustAGORA4.py' ' --ucrFolder ' ucrFolder ' --modelnameslist ' AGORAModel1.description ',' AGORAModel2.description ' --modelnumberslist ' num2str(i) ',' num2str(j) ' --isSpeciesMerged False' ' --IBDFlag True']);
			else
			    %system(['python /mnt/vdb/home/ubuntu2/MATLAB/SEED/src/simulateSmallModelsSeparatePicrustAGORA4.py' ' --ucrFolder ' ucrFolder ' --modelnameslist ' AGORAModel1.description ',' AGORAModel2.description ' --modelnumberslist ' num2str(i) ',' num2str(j) ' --isSpeciesMerged False' ' --IBDFlag False']);
			end
			%status = importdata(['/mnt/vdb/home/ubuntu2/' ucrFolder AGORAModel1.description AGORAModel2.description 'False.txt']);
		    end
		end

		if doSimulation
		    if exist(['/mnt/vdb/home/ubuntu2/' ucrFolder AGORAModel1.description 'False.txt'],'file')
			status = importdata(['/mnt/vdb/home/ubuntu2/' ucrFolder AGORAModel1.description 'False.txt']);
		    end
		    if exist([inputDir filesep 'MGMData/IBDCommon/' 'speciesSep' num2str(i) '.modelexpr'],'file') && exist([inputDir filesep 'MGMData/IBDCommon/' 'speciesSep' num2str(j) '.modelexpr'],'file')
			status = 0;
		    end
		    if status==0 || useFBA==1 || realOneSpecies==1
			if realOneSpecies==0
                            if useCommon
                                species1expr = simulateFuncFlux(useFBA,[inputDir filesep 'MGMData/IBDCommon/' 'speciesSep' num2str(i) '.modelexpr'],outputDir1,i,j,ucrFolder,AGORAModel1,'_1_',doSimulation);
                                simulateFuncFlux(useFBA,[inputDir filesep 'MGMData/IBDCommon/' 'speciesSep' num2str(j) '.modelexpr'],outputDir1,i,j,ucrFolder,AGORAModel2,'_2_',doSimulation);
                            else
                                species1expr = simulateFuncFlux(useFBA,[inputDir filesep 'MGMData/' ucrFolder '/' 'speciesSep' num2str(i) '.modelexpr'],outputDir1,i,j,ucrFolder,AGORAModel1,'_1_',doSimulation);
                                simulateFuncFlux(useFBA,[inputDir filesep 'MGMData/' ucrFolder '/' 'speciesSep' num2str(j) '.modelexpr'],outputDir1,i,j,ucrFolder,AGORAModel2,'_2_',doSimulation);
                            end
			    if ~doSimulation
				AGORAModelWithExpr = AGORAModel1;
				AGORAModelWithExpr.expression = species1expr;
			    end
			    if singleSpeciesOnly
				AGORAModelWithExpr.expression = species1expr;
				break;
			    end
			    else
				species1expr = simulateFuncFlux(useFBA,[inputDir filesep 'MGMData/' ucrFolder '/' 'speciesSep' num2str(i) '.modelexpr'],outputDir1,i,j,ucrFolder,AGORAModel1,'_1_',doSimulation);
                                simulateFuncFlux(useFBA,[inputDir filesep 'MGMData/' ucrFolder '/' 'speciesSep' num2str(j) '.modelexpr'],outputDir1,i,j,ucrFolder,AGORAModel2,'_2_',doSimulation);
			end
		    end
		end
	    end
	end
    end
end
