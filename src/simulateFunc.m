function [AGORAModelWithExpr] = simulateFunc(mergeModelsAGORAFlag,useFBA,singleSpeciesOnly,inputDir,outputDir1,ucrFolder,i,j1,blocksize,AGORAMat,AGORAModelArr,presentSpecies,realOneSpecies,doSimulation)

AGORAModelWithExpr = [];
disp(['i ' num2str(i)])
if realOneSpecies
    if ~doSimulation
	if ~isempty(AGORAMat{i})
	    isPresent = 0;
	    for k=1:length(presentSpecies)
		if ~isempty(regexp(AGORAMat{i},presentSpecies{k}))
		    isPresent = 1;
		end
	    end
	    if isPresent==1
		disp(['found ' num2str(i)])
		AGORAModel = AGORAModelArr{i};
		AGORAModelWithExpr = AGORAModel;
		status = 1;
		if useFBA==0
		    writeData({AGORAModel.rxnECNumbers},[inputDir filesep 'MGMData/' ucrFolder '/speciesSep' num2str(i) '.modelec'],'\t');
                    system(['python /mnt/vdb/home/ubuntu2/MATLAB/SEED/src/simulateSmallModelsSeparatePicrustAGORA.py False ' ucrFolder ' ' AGORAModel.description ' ' num2str(i)]);
		    status = importdata(['/mnt/vdb/home/ubuntu2/' ucrFolder AGORAModel.description '.txt']);
		end
            end
	end
    else
	for j=(j1-1)*blocksize+1:min(length(AGORAMat),j1*blocksize)
	    if i~=j && ~isempty(AGORAMat{i}) && ~isempty(AGORAMat{j})
		isPresent1 = 0;
		isPresent2 = 0;
		for k=1:length(presentSpecies)
		    if ~isempty(regexp(AGORAMat{i},presentSpecies{k}))
			isPresent1 = 1;
		    end
		    if ~isempty(regexp(AGORAMat{j},presentSpecies{k}))
			isPresent2 = 1;
		    end
		end
		if isPresent1==1 && isPresent2==1   
		disp(['found ' num2str(i)])
		disp(['found ' num2str(j)])
		%load(AGORAMat{i});
		AGORAModel1 = AGORAModelArr{i};
		AGORAModelWithExpr = AGORAModel1;
		%load(AGORAMat{j});
		AGORAModel2 = AGORAModelArr{j};
		if status==0 || useFBA==1
		    if useFBA==0
			speciesexpr = importdata([inputDir filesep 'MGMData/' ucrFolder '/species.modelexpr']);
			expressionData = speciesexpr;
		    else
			expressionData = ones(size(AGORAModel.rxns));
		    end			    
		    expressionSDs = ones(size(AGORAModel.rxns)).*min(expressionData,1);
		    expressionIDs = AGORAModel.rxns;
		    AGORAModel.rxnECNums = {};
		    for k=1:length(AGORAModel.rxns)
			AGORAModel.rxnECNums{end+1} = {};
		    end
		    AGORAModel = fluxModelFunc(AGORAModel);
		    AGORAModel = assignSortedBiom(AGORAModel);
		    try
			if useFBA
			    %picrustFluxes = runFluxMethod(expressionData,expressionIDs,'testfalcon',AGORAModel,'FBA',expressionSDs);
			    %writeData({picrustFluxes},[outputDir1 filesep num2str(i) '_' num2str(j) '_' ucrFolder '_FBA.flux'],'\t');
			else
		% 		picrustFluxes = runFluxMethod(expressionData,expressionIDs,'testfalcon',AGORAModel,'FALCON',expressionSDs);
		% 		writeData({picrustFluxes3},[outputDir1 filesep num2str(i) '_' num2str(j) '_' ucrFolder '.flux'],'\t');
			end
		    catch
		% 	    disp('PICRUST FAIL')
		    end
		end
		end
	    end
	end
    end
else
    for j=(j1-1)*blocksize+1:min(length(AGORAMat),j1*blocksize)
	if i~=j && ~isempty(AGORAMat{i}) && ~isempty(AGORAMat{j})
	    isPresent1 = 0;
	    isPresent2 = 0;
	    for k=1:length(presentSpecies)
		if ~isempty(regexp(AGORAMat{i},presentSpecies{k}))
		    isPresent1 = 1;
		end
		if ~isempty(regexp(AGORAMat{j},presentSpecies{k}))
		    isPresent2 = 1;
		end
	    end

	    if isPresent1==1 && isPresent2==1   
		disp(['found ' num2str(i)])
		disp(['found ' num2str(j)])
		%load(AGORAMat{i});
		AGORAModel1 = AGORAModelArr{i};
		AGORAModelWithExpr = AGORAModel1;
		%load(AGORAMat{j});
		AGORAModel2 = AGORAModelArr{j};
		if mergeModelsAGORAFlag
		    AGORAModel3 = mergeModelsAGORA(AGORAModel1,AGORAModel2);
		    status = 1;
		    if useFBA==0
			writeData({AGORAModel3.rxnECNumbers},[inputDir filesep 'MGMData/' ucrFolder '/species3.modelec'],'\t');
			system(['python /mnt/vdb/home/ubuntu2/MATLAB/SEED/src/simulateSmallModelsSeparatePicrustAGORA.py True ' ucrFolder ' ' AGORAModel1.description ' ' AGORAModel2.description ' True']);
			status = importdata(['/mnt/vdb/home/ubuntu2/' ucrFolder AGORAModel1.description AGORAModel2.description 'True.txt']);
		    end

		    if status==0 || useFBA==1
			if useFBA==0
			    species3expr = importdata([inputDir filesep 'MGMData/' ucrFolder '/species3.modelexpr']);
			    expressionData = species3expr;
			else
			    expressionData = ones(size(AGORAModel3.rxns));
			end			    
			expressionSDs = ones(size(AGORAModel3.rxns)).*min(expressionData,1);
			expressionIDs = AGORAModel3.rxns;
			AGORAModel3.rxnECNums = {};
			for k=1:length(AGORAModel3.rxns)
			    AGORAModel3.rxnECNums{end+1} = {};
			end
			AGORAModel3 = fluxModelFunc(AGORAModel3);
			AGORAModel3 = assignSortedBiom(AGORAModel3);
			try
			    if useFBA
				picrustFluxes3 = runFluxMethod(expressionData,expressionIDs,'testfalcon',AGORAModel3,'FBA',expressionSDs);
				writeData({picrustFluxes3},[outputDir1 filesep num2str(i) '_' num2str(j) '_3_' ucrFolder '_FBA.flux'],'\t');
			    else
				picrustFluxes3 = runFluxMethod(expressionData,expressionIDs,'testfalcon',AGORAModel3,'FALCON',expressionSDs);
				writeData({picrustFluxes3},[outputDir1 filesep num2str(i) '_' num2str(j) '_3_' ucrFolder '.flux'],'\t');
			    end
			catch
			    disp('PICRUST 3 FAIL')
			end
		    end

		else
		    status = 1;
		    if useFBA==0 && realOneSpecies==0
			writeData({AGORAModel1.rxnECNumbers},[inputDir filesep 'MGMData/' ucrFolder filesep AGORAModel1.description AGORAModel2.description 'species1.modelec'],'\t');
			writeData({AGORAModel2.rxnECNumbers},[inputDir filesep 'MGMData/' ucrFolder filesep AGORAModel1.description AGORAModel2.description 'species2.modelec'],'\t');
			system(['python /mnt/vdb/home/ubuntu2/MATLAB/SEED/src/simulateSmallModelsSeparatePicrustAGORA.py True ' ucrFolder ' ' AGORAModel1.description ' ' AGORAModel2.description ' False']);
			status = importdata(['/mnt/vdb/home/ubuntu2/' ucrFolder AGORAModel1.description AGORAModel2.description 'False.txt']);
		    end

		    if status==0 || useFBA==1 || realOneSpecies==1
			if useFBA==0
			    species1expr = importdata([inputDir filesep 'MGMData/' ucrFolder filesep AGORAModel1.description AGORAModel2.description 'species1.modelexpr']);
			    species2expr = importdata([inputDir filesep 'MGMData/' ucrFolder filesep AGORAModel1.description AGORAModel2.description 'species2.modelexpr']);
			    if singleSpeciesOnly
				AGORAModelWithExpr.expression = species1expr;
				break;
			    end
			    expressionData = species1expr;
			else
			    expressionData = ones(size(AGORAModel1.rxns));
			end

			AGORAModel1.description
			AGORAModel2.description
			size(expressionData)
			size(AGORAModel1.rxns)
			expressionSDs = ones(size(AGORAModel1.rxns)).*min(expressionData,1);
			expressionIDs = AGORAModel1.rxns;
			AGORAModel1.rxnECNums = {};
			for k=1:length(AGORAModel1.rxns)
			    AGORAModel1.rxnECNums{end+1} = {};
			end
			AGORAModel1 = fluxModelFunc(AGORAModel1);
			AGORAModel1 = assignSortedBiom(AGORAModel1);
			try
			    if useFBA
				picrustFluxes1 = runFluxMethod(expressionData,expressionIDs,'testfalcon',AGORAModel1,'FBA',expressionSDs);
				writeData({picrustFluxes1},[outputDir1 filesep num2str(i) '_' num2str(j) '_1_' ucrFolder '_FBA.flux'],'\t');
			    else
				picrustFluxes1 = runFluxMethod(expressionData,expressionIDs,'testfalcon',AGORAModel1,'FALCON',expressionSDs);
				writeData({picrustFluxes1},[outputDir1 filesep num2str(i) '_' num2str(j) '_1_' ucrFolder '.flux'],'\t');
			    end
			catch
			    disp('PICRUST 1 FAIL')
			    %nonsense = nonsense+1;
			end

			if useFBA==0
			    expressionData = species2expr;
			else                                               expressionData = ones(size(AGORAModel2.rxns));
			end
			expressionSDs = ones(size(AGORAModel2.rxns)).*min(expressionData,1);
			expressionIDs = AGORAModel2.rxns;
			AGORAModel2.rxnECNums = {};
			for k=1:length(AGORAModel2.rxns)
			    AGORAModel2.rxnECNums{end+1} = {};
			end
			AGORAModel2 = fluxModelFunc(AGORAModel2);
			AGORAModel2 = assignSortedBiom(AGORAModel2);

			try
			    if useFBA
				picrustFluxes2 = runFluxMethod(expressionData,expressionIDs,'testfalcon',AGORAModel2,'FBA',expressionSDs);
				writeData({picrustFluxes2},[outputDir1 filesep num2str(i) '_' num2str(j) '_2_' ucrFolder '_FBA.flux'],'\t');
			    else
				picrustFluxes2 = runFluxMethod(expressionData,expressionIDs,'testfalcon',AGORAModel2,'FALCON',expressionSDs);
				writeData({picrustFluxes2},[outputDir1 filesep num2str(i) '_' num2str(j) '_2_' ucrFolder '.flux'],'\t');
			    end
			catch
			    disp('PICRUST 2 FAIL')
			    %nonsense = nonsense+1;
			end
		    end
		end
	    end
	end
    end
end
