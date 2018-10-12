function [AGORAModelWithExpr] = simulateFunc(mergeModelsAGORAFlag,useFBA,singleSpeciesOnly,inputDir,outputDir1,ucrFolder,i,j1,blocksize,AGORAMat,AGORAModelArr,presentSpecies,realOneSpecies,doSimulation,justTwoSpecies,IBDFlag,justAGORAIBDMap)

if ~exist('IBDFlag','var')
    IBDFlag = 0
end
AGORAModelWithExpr = [];
disp(['i ' num2str(i)])
    for j=(j1-1)*blocksize+1:min(length(AGORAMat),j1*blocksize)
	if i~=j && ~isempty(AGORAMat{i}) && ~isempty(AGORAMat{j})
	    isPresent1 = 0;
	    isPresent2 = 0;
            if IBDFlag
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
            else
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
		    AGORAModel3 = mergeModelsAGORA(AGORAModel1,AGORAModel2);
		    status = 1;
		    if useFBA==0
			writeData({AGORAModel3.rxnECNumbers},[inputDir filesep 'MGMData/' ucrFolder '/species3.modelec'],'\t');
			system(['python /mnt/vdb/home/ubuntu2/MATLAB/SEED/src/simulateSmallModelsSeparatePicrustAGORA3.py' ' --isTwoSpecies True' ' --ucrFolder ' ucrFolder ' --modelname1 ' AGORAModel1.description ' --modelname2 ' AGORAModel2.description ' --isSpecies3 True']);
			status = importdata(['/mnt/vdb/home/ubuntu2/' ucrFolder AGORAModel1.description AGORAModel2.description 'True.txt']);
		    end

		    if status==0 || useFBA==1
		        simulateFuncFlux(useFBA,[inputDir filesep 'MGMData/' ucrFolder '/species3.modelexpr'],AGORAModel3,'_3_');
		    end

		else
		    status = 1;
		    if useFBA==0 && realOneSpecies==0
			writeData({AGORAModel1.rxnECNumbers},[inputDir filesep 'MGMData/' ucrFolder filesep AGORAModel1.description AGORAModel2.description 'species1.modelec'],'\t');
			writeData({AGORAModel2.rxnECNumbers},[inputDir filesep 'MGMData/' ucrFolder filesep AGORAModel1.description AGORAModel2.description 'species2.modelec'],'\t');
                        if ~doSimulation
			    if justTwoSpecies
			        system(['python /mnt/vdb/home/ubuntu2/MATLAB/SEED/src/simulateSmallModelsSeparatePicrustAGORA3.py' ' --isTwoSpecies False' ' --ucrFolder ' ucrFolder ' --modelname ' AGORAModel1.description ' --modelnumber ' num2str(i) ' --justTwoSpecies True']);
			    else
			        system(['python /mnt/vdb/home/ubuntu2/MATLAB/SEED/src/simulateSmallModelsSeparatePicrustAGORA3.py' ' --isTwoSpecies False' ' --ucrFolder ' ucrFolder ' --modelname ' AGORAModel1.description ' --modelnumber ' num2str(i) ' --justTwoSpecies False']);
			    end
                        else
			    if IBDFlag
				system(['python /mnt/vdb/home/ubuntu2/MATLAB/SEED/src/simulateSmallModelsSeparatePicrustAGORA3.py' ' --isTwoSpecies True' ' --ucrFolder ' ucrFolder ' --modelname1 ' AGORAModel1.description ' --modelname2 ' AGORAModel2.description ' --isSpecies3 False' ' --IBDFlag True']);
			    else
				system(['python /mnt/vdb/home/ubuntu2/MATLAB/SEED/src/simulateSmallModelsSeparatePicrustAGORA3.py' ' --isTwoSpecies True' ' --ucrFolder ' ucrFolder ' --modelname1 ' AGORAModel1.description ' --modelname2 ' AGORAModel2.description ' --isSpecies3 False' ' --IBDFlag False']);
			    end
			end
			status = importdata(['/mnt/vdb/home/ubuntu2/' ucrFolder AGORAModel1.description AGORAModel2.description 'False.txt']);
		    end

		    if status==0 || useFBA==1 || realOneSpecies==1 || ~doSimulation
			if realOneSpecies==0
		            species1expr = simulateFuncFlux(useFBA,[inputDir filesep 'MGMData/' ucrFolder '/' AGORAModel1.description AGORAModel2.description 'species1.modelexpr'],outputDir1,i,j,ucrFolder,AGORAModel1,'_1_',doSimulation);
                            simulateFuncFlux(useFBA,[inputDir filesep 'MGMData/' ucrFolder '/' AGORAModel1.description AGORAModel2.description 'species2.modelexpr'],outputDir1,i,j,ucrFolder,AGORAModel2,'_2_',doSimulation);
                            if ~doSimulation
                                AGORAModelWithExpr = AGORAModel1;
                                AGORAModelWithExpr.expression = species1expr;
                            end
			    if singleSpeciesOnly
				AGORAModelWithExpr.expression = species1expr;
				break;
			    end
			else
			    simulateFuncFlux(useFBA,[inputDir filesep 'MGMData/' ucrFolder filesep 'speciesSep' num2str(i) '.modelexpr'],outputDir1,i,j,ucrFolder,AGORAModel1,'_1_',doSimulation);
                            simulateFuncFlux(useFBA,[inputDir filesep 'MGMData/' ucrFolder filesep 'speciesSep' num2str(j) '.modelexpr'],outputDir1,i,j,ucrFolder,AGORAModel1,'_2_',doSimulation);    
			end
		    end
		end
	    end
	end
    end
