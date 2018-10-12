useNorm = 1;
outputDir1 = [outputDir filesep 'simulateSmallModelsSeparatePicrustAGORA'];
if useNorm
    outputDir1 = [outputDir filesep 'simulateSmallModelsSeparatePicrustAGORANormObeseNorm'];
end
if ~exist(outputDir1,'dir')
    system(['mkdir ' outputDir1]);
end

useFBA = 0;
blocksize = 800;
singleSpeciesOnly = 0;
AGORAModelWithExprArr = {};
count = 0;
realOneSpecies = 1;
doSimulation = 0;
justTwoSpecies = 1;
useMergedModel = 1;
picrustFluxesArrNormObese = {};
picrustFluxesArrNormObeseEFlux = {};
useEFlux = 0;
for z=1:length(ucrFolders)
    ucrFolder = ucrFolders{z};
    if useMergedModel
        biomassrxns = tobemerged.rxns(cellfun(@(x) ~isempty(regexp(x,'biomass')), tobemerged.rxns));
        biomassrxnNames = tobemerged.rxnNames(cellfun(@(x) ~isempty(regexp(x,'biomass')), tobemerged.rxns));
        tobemerged = changeObjective(tobemerged,biomassrxns,ones(length(biomassrxns),1));
        writeData({tobemerged.rxnECNumbers},[inputDir filesep 'MGMData/' ucrFolder filesep 'normalized_otus' '.modelec'],'\t');
        if ~exist([inputDir '/MGMData/' ucrFolder '/normalized_otus.modelexpr'],'file')
            system(['python /mnt/vdb/home/ubuntu2/MATLAB/SEED/src/simulateSmallModelsSeparatePicrustAGORA5.py --ucrFolder ' ucrFolder])
        end
        ac = importdata([inputDir '/MGMData/' ucrFolder '/normalized_otus.modelexpr']);
        if length(ac)~=length(tobemerged.rxns)
            system(['python /mnt/vdb/home/ubuntu2/MATLAB/SEED/src/simulateSmallModelsSeparatePicrustAGORA5.py --ucrFolder ' ucrFolder])
            ac = importdata([inputDir '/MGMData/' ucrFolder '/normalized_otus.modelexpr']);			end
	if useNorm
            ac = ac/sum(ac);
        end
	if useEFlux
 	    picrustFluxesArrNormObeseEFlux{z} = runFluxMethod(ac,tobemerged.rxns,'testeflux',tobemerged,'EFlux',ones(1,length(tobemerged.rxns)),biomassrxnNames{1});
            writeData({picrustFluxesArrNormObeseEFlux{z}},['/mnt/vdb/home/ubuntu2/' ucrFolder 'efluxmerged.flux'],'\t');
	else
            picrustFluxesArrNormObese{z} = runFluxMethod(ac,tobemerged.rxns,'testfalcon',tobemerged,'FALCON',ones(1,length(tobemerged.rxns)));
            if useNorm
	        writeData({picrustFluxesArrNormObese{z}},['/mnt/vdb/home/ubuntu2/' ucrFolder 'normobesenormmerged.flux'],'\t');
            else
                writeData({picrustFluxesArrNormObese{z}},['/mnt/vdb/home/ubuntu2/' ucrFolder 'merged.flux'],'\t');
            end
        end
    else
	system(['python /mnt/vdb/home/ubuntu2/pickPresent.py ' inputDir filesep 'MGMData' filesep ucrFolder filesep 'normalized_otus.tsv']);
	presentSpecies = textscan(fopen('/mnt/vdb/home/ubuntu2/pickPresent.txt'),'%s','Delimiter','\n','HeaderLines',0);
	presentSpecies = presentSpecies{1};
	mergeModelsAGORAFlag = 0;
	AGORAModelArr = {};
	AGORAFoundIdxs = [];
	for i=1:length(AGORAMat)
	    for k=1:length(presentSpecies)
		if ~isempty(regexp(AGORAMat{i},presentSpecies{k}))
		    load(AGORAMat{i});
		    i
		    AGORAModelArr{i} = AGORAModel;
		    AGORAFoundIdxs(end+1) = i;
		end
	    end
	end

	for i1=1:ceil(length(AGORAMat)/blocksize)
	    for j1=1:ceil(length(AGORAMat)/blocksize)
		for i=(i1-1)*blocksize+1:min(length(AGORAMat),i1*blocksize)
		    if sum(AGORAFoundIdxs==i)~=0
			changeCobraSolver('glpk');
			simulateFuncTemp2(mergeModelsAGORAFlag,useFBA,singleSpeciesOnly,inputDir,outputDir1,ucrFolder,i,j1,blocksize,AGORAMat,AGORAModelArr,presentSpecies,realOneSpecies,doSimulation,justTwoSpecies);
		    end
		    if 0
			foundFlag = 0;
			if ~isempty(AGORAModelWithExpr) && isfield(AGORAModelWithExpr,'expression') && singleSpeciesOnly
			    AGORAModelWithExprArr{z,i} = AGORAModelWithExpr;
			end
		    end
		end
	    end
	end
    end
end
