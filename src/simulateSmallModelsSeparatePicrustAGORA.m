outputDir1 = [outputDir filesep 'simulateSmallModelsSeparatePicrustAGORA'];

useFBA = 0;
blocksize = 150;
singleSpeciesOnly = 0;
AGORAModelWithExprArr = {};
count = 0;
realOneSpecies = 1;
doSimulation = 0;
for z=1:length(ucrFolders)
    ucrFolder = ucrFolders{z};
    system(['python /mnt/xvdf/home/ubuntu2/pickPresent.py ' inputDir filesep 'MGMData' filesep ucrFolder filesep 'normalized_otus.tsv']);
    presentSpecies = textscan(fopen('/mnt/xvdf/home/ubuntu2/pickPresent.txt'),'%s','Delimiter','\n','HeaderLines',0);
    presentSpecies = presentSpecies{1};
    mergeModelsAGORAFlag = 0;
    for i1=1:ceil(length(AGORAMat)/blocksize)
	for j1=1:ceil(length(AGORAMat)/blocksize)
	    if 1
	        AGORAModelArr = {};
	        for i=1:length(AGORAMat)
	    	    if ~isempty(AGORAMat{i}) && ( (i>=(i1-1)*blocksize && i<=i1*blocksize) || (i>=(j1-1)*blocksize && i<=j1*blocksize) ) 
	    	        load(AGORAMat{i});
	    	        i
	    	        AGORAModelArr{i} = AGORAModel;
	    	    end
	        end
	    end
	    parfor i=(i1-1)*blocksize+1:min(length(AGORAMat),i1*blocksize)
		%changeCobraSolver('glpk');
AGORAModelWithExpr = simulateFunc(mergeModelsAGORAFlag,useFBA,singleSpeciesOnly,inputDir,outputDir1,ucrFolder,i,j1,blocksize,AGORAMat,AGORAModelArr,presentSpecies,realOneSpecies,doSimulation);
                if i>50
                    %nonsense = nonsense+1;
                end
                foundFlag = 0;
		if ~isempty(AGORAModelWithExpr) && isfield(AGORAModelWithExpr,'expression') && singleSpeciesOnly
		    AGORAModelWithExprArr{z,i} = AGORAModelWithExpr;
		end
	    end
	end
    end
end
