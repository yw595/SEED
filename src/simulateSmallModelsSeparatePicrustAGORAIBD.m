useDiabetes = 1;
useNorm = 1;
if useDiabetes
    outputDir1 = [outputDir filesep 'simulateSmallModelsSeparatePicrustAGORADiabetes'];
    if useNorm
        outputDir1 = [outputDir filesep 'simulateSmallModelsSeparatePicrustAGORADiabetesNorm'];
    end
else
    outputDir1 = [outputDir filesep 'simulateSmallModelsSeparatePicrustAGORAIBD'];
end
if ~exist(outputDir1,'dir')
    system(['mkdir ' outputDir1]);
end


runSim = 1;
extractDiffAbund = 0;
writeDiffAbund = 0;
if runSim==1
    useMergedModel = 1;
    useFBA = 0;
    blocksize = 4;
    singleSpeciesOnly = 0;
    AGORAModelWithExprArr = {};
    count = 0;
    realOneSpecies = 0;
    doSimulation = 1;
    justTwoSpecies = 0;
    if useDiabetes
        taxonomy_profiles = importdata('/mnt/vdb/home/ubuntu2/taxonomyDiabetes.tsv');
        taxonomy_profilestemp = {};
        taxonomy_profilestemp{1} = taxonomy_profiles.textdata{1,1};
        for j=2:size(taxonomy_profiles.textdata,2)
	    taxonomy_profilestemp{1} = [taxonomy_profilestemp{1} sprintf(['\t' taxonomy_profiles.textdata{1,j}])];
        end
	for i=2:size(taxonomy_profiles.textdata,1)
	    taxonomy_profilestemp{i} = taxonomy_profiles.textdata{i,1};
        end
        for i=1:size(taxonomy_profiles.data,1)
	    for j=1:size(taxonomy_profiles.data,2)
		taxonomy_profilestemp{i+1} = [taxonomy_profilestemp{i+1} sprintf(['\t' num2str(taxonomy_profiles.data(i,j))])];
            end
        end
	taxonomy_profiles = taxonomy_profilestemp;
    else
        taxonomy_profiles = importdata([inputDir '/taxonomic_profiles.tsv']);
    end
    sampleNames = taxonomy_profiles{1};
    sampleNames = strsplit(sampleNames);
    if useDiabetes
        sampleNames = sampleNames(3:end);
    else
        sampleNames = sampleNames(3:end-1);
    end
    speciesNames = {};
    speciesAbunds = [];
    justAGORAIBDMap = containers.Map;
    if useDiabetes
        justAGORAIBD = textscan(fopen('/mnt/vdb/home/ubuntu2/justAGORADiabetesDists.txt'),'%s%s','Delimiter','|');
    else
        justAGORAIBD = textscan(fopen('/mnt/vdb/home/ubuntu2/justAGORAIBDDists.txt'),'%s%s','Delimiter','|');
    end
    for i=1:length(justAGORAIBD{1})
	justAGORAIBDMap(justAGORAIBD{1}{i}) = justAGORAIBD{2}{i};
    end
    for i=2:length(taxonomy_profiles)
	row = taxonomy_profiles{i};
	row = strsplit(row,'\t');
        if useDiabetes
	    speciesAbunds(i-1,:) = cellfun(@(x) str2num(x), row(2:end));
	    speciesNames{i-1} = row{1};
        else
	    speciesAbunds(i-1,:) = cellfun(@(x) str2num(x), row(2:end-1));
	    speciesNames{i-1} = row{end};
        end
    end
    AGORAModelAbundMat = zeros(length(sampleNames),length(AGORAMat));
    for z=1:length(sampleNames)
	sampleName = sampleNames{z};
	z
	if useDiabetes
	    if ~exist([inputDir filesep 'MGMData/Diabetes/' sampleName],'dir')
		mkdir([inputDir filesep 'MGMData/Diabetes/' sampleName]);
	    end
	else
	    if ~exist([inputDir filesep 'MGMData/' sampleName],'dir')
		mkdir([inputDir filesep 'MGMData/' sampleName]);
	    end
	end
	presentSpecies = speciesNames(speciesAbunds(:,z)~=0);
	presentIdxs = find(speciesAbunds(:,z)~=0);
	mergeModelsAGORAFlag = 0;
	IBDFlag = 1;
	%AGORAModelArr = {};
	for i=1:length(AGORAMat)
	    %if ~isempty(AGORAMat{i}) && ( (i>=(i1-1)*blocksize && i<=i1*blocksize) || (i>=(j1-1)*blocksize && i<=j1*blocksize) )
	        %load(AGORAMat{i});
	        %i
	        %AGORAModelArr{i} = AGORAModel;
            %end
	end
	AGORAFoundIdxs = [];
        for i=1:length(AGORAMat)
	    AGORAName = strsplit(AGORAMat{i},'/');
	    AGORAName = strsplit(AGORAName{end},'_');
	    AGORAName = [AGORAName{1} ' ' AGORAName{2}];
	    for k=1:length(presentSpecies)
		if isKey(justAGORAIBDMap,AGORAName) && ~isempty(regexp(presentSpecies{k},justAGORAIBDMap(AGORAName)))
		    AGORAModelAbundMat(z,i) = AGORAModelAbundMat(z,i) + speciesAbunds(presentIdxs(k),z);
		    AGORAFoundIdxs(end+1) = i;
		end
	    end
	end

	if useDiabetes
	    normFI = fopen([inputDir filesep 'MGMData/Diabetes/' sampleName filesep 'normalized_otus.tsv'],'w');
	else
	    normFI = fopen([inputDir filesep 'MGMData/' sampleName filesep 'normalized_otus.tsv'],'w');
        end
        fprintf(normFI,'# Constructed from excel file\n')
        fprintf(normFI,'#OTU ID DUMMY\n')
        AGORADistsManualTemp = textscan(fopen('/mnt/vdb/home/ubuntu2/justAGORADistsManual.txt'),'%s%s','Delimiter','|');
        ggIDsToTaxTemp = textscan(fopen('/mnt/vdb/home/ubuntu2/ggIDsToTax.txt'),'%s%s','Delimiter','\t');
        AGORADistsManual = containers.Map;
        for i=1:length(AGORADistsManualTemp{1})
	    AGORADistsManual(AGORADistsManualTemp{1}{i}) = AGORADistsManualTemp{2}{i};
        end
        ggIDsToTax = containers.Map;
        for i=1:length(ggIDsToTaxTemp{2})
	    ggIDsToTax(ggIDsToTaxTemp{2}{i}) = ggIDsToTaxTemp{1}{i};
        end
	ggAbundsMap = containers.Map;
	for i=1:length(AGORAMat)
	    AGORAName = strsplit(AGORAMat{i},'/');
	    AGORAName = strsplit(AGORAName{end},'_');
	    AGORAName = [AGORAName{1} ' ' AGORAName{2}];
            if isKey(AGORADistsManual,AGORAName)
                taxid = AGORADistsManual(AGORAName);
                if isKey(ggIDsToTax,taxid)
                    ggid = ggIDsToTax(taxid);
                    if AGORAModelAbundMat(z,i)~=0
		        if isKey(ggAbundsMap,ggid)
		            ggAbundsMap(ggid) = ggAbundsMap(ggid) + AGORAModelAbundMat(z,i);
		        else
		            ggAbundsMap(ggid) = AGORAModelAbundMat(z,i);
                        end
                    end
                end
            end
        end
	ggids = keys(ggAbundsMap);
        for i=1:length(ggids)
	    ggid = ggids{i};
	    fprintf(normFI,'%s\t%s\n',ggid,num2str(ggAbundsMap(ggids{i})));
        end
	fclose(normFI);


	if useMergedModel
	    parentFolder = '';
            if useDiabetes
                parentFolder = 'Diabetes/';
            end
	    writeData({tobemerged.rxnECNumbers},[inputDir filesep 'MGMData/' parentFolder sampleNames{z} filesep 'normalized_otus' '.modelec'],'\t');
	    if ~exist([inputDir '/MGMData/' parentFolder sampleNames{z} '/normalized_otus.modelexpr'],'file')
		system(['python /mnt/vdb/home/ubuntu2/MATLAB/SEED/src/simulateSmallModelsSeparatePicrustAGORA5.py --ucrFolder ' parentFolder sampleNames{z}])
	    end
	    ac = importdata([inputDir '/MGMData/' parentFolder sampleNames{z} '/normalized_otus.modelexpr']);
            if useNorm
                ac = ac/sum(ac);
            end
	    picrustFluxesArr{z} = runFluxMethod(ac,tobemerged.rxns,'testfalcon',tobemerged,'FALCON',ones(1,length(tobemerged.rxns)));
            if useDiabetes
	        if useNorm
	            writeData({picrustFluxesArr{z}},['/mnt/vdb/home/ubuntu2/' sampleNames{z} 'diabetesnormmerged.flux'],'\t');
	        else
		    writeData({picrustFluxesArr{z}},['/mnt/vdb/home/ubuntu2/' sampleNames{z} 'diabetesmerged.flux'],'\t');
	        end
            else
	        writeData({picrustFluxesArr{z}},['/mnt/vdb/home/ubuntu2/' sampleNames{z} 'merged.flux'],'\t');
            end
	else
	    presentModelIdxs = [];
	    for i=1:length(AGORAMat)
		if ~any(AGORAFoundIdxs==i)
		    presentModelIdxs(i) = 0;
		else
		    presentModelIdxs(i) = 1;
		end
	    end
	    for i1=1:ceil(length(AGORAMat)/blocksize)
		for j1=1:ceil(length(AGORAMat)/blocksize)
		    if 1
			AGORAModelArr = {};
			for i=1:length(AGORAMat)
			    if ~isempty(AGORAMat{i}) && ( (i>=(i1-1)*blocksize && i<=i1*blocksize) || (i>=(j1-1)*blocksize && i<=j1*blocksize) ) 
				%load(AGORAMat{i});
				%i
				AGORAModelArr{i} = AGORAModelArrFull{i};
			    end
			end
		    end
		    for i=(i1-1)*blocksize+1:min(length(AGORAMat),i1*blocksize)
			if ~exist([inputDir filesep 'MGMData/IBDCommon' filesep 'speciesSep' num2str(i) '.modelexpr'],'file')
			    disp('where2')
			    disp(exist([inputDir filesep 'MGMData/IBDCommon' filesep 'speciesSep' num2str(i) '.modelexpr'],'file'))
			    simulateFuncTemp4(mergeModelsAGORAFlag,useFBA,singleSpeciesOnly,inputDir,outputDir1,sampleName,i,j1,blocksize,AGORAMat,AGORAModelArr,presentModelIdxs,realOneSpecies,0,justTwoSpecies,IBDFlag,justAGORAIBDMap);
                        end
		    end
		end
	    end

	    presentModelIdxsTrue = find(presentModelIdxs);
	    ithRandPerm = randperm(length(presentModelIdxsTrue),1);
	    jthRandPerm = randperm(length(AGORAMat),100);
	    for i1=1:ceil(length(AGORAMat)/blocksize)
		for j1=1:ceil(length(AGORAMat)/blocksize)
		    parfor i=(i1-1)*blocksize+1:min(length(AGORAMat),i1*blocksize)
			    %if sum(i==ithRandPerm)~=0
			    %changeCobraSolver('glpk');
			    simulateFuncTemp4(mergeModelsAGORAFlag,useFBA,singleSpeciesOnly,inputDir,outputDir1,sampleName,i,j1,blocksize,AGORAMat,AGORAModelArr,presentModelIdxs,realOneSpecies,doSimulation,justTwoSpecies,IBDFlag,justAGORAIBDMap);
			    %end
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
elseif extractDiffAbund==0
    uniqSubsystems = containers.Map;
    for i=1:length(AGORAMat)
	i
	load(AGORAMat{i})
	for j=1:length(AGORAModel.subSystems)
	    uniqSubsystems(AGORAModel.subSystems{j}{1}) = '';
        end
    end
    uniqSubsystems = sort(keys(uniqSubsystems));
    diffAbunds = zeros(length(uniqSubsystems),2);
    for j=1:size(AGORAModelAbundMat,2)
	j
	load(AGORAMat{j})
	for i=1:size(AGORAModelAbundMat,1)
	    for k=1:length(AGORAModel.subSystems)
		if i<=size(AGORAModelAbundMat,1)/2
		    diffAbunds(strcmp(uniqSubsystems,AGORAModel.subSystems{k}),1) = diffAbunds(strcmp(uniqSubsystems,AGORAModel.subSystems{k}),1) + AGORAModelAbundMat(i,j);
		else
		    diffAbunds(strcmp(uniqSubsystems,AGORAModel.subSystems{k}),2) = diffAbunds(strcmp(uniqSubsystems,AGORAModel.subSystems{k}),2) + AGORAModelAbundMat(i,j);
		end
	    end
	end
    end
elseif writeDiffAbund==0
    diffAbunds2 = diffAbunds(:,1)-diffAbunds(:,2);
    [~, sortIdxs] = sort(diffAbunds2,'descend');
    sortedSubsystems = uniqSubsystems(sortIdxs);
    writeData({addIdxStrings(sortedSubsystems(1:20)),diffAbunds2(sortIdxs(1:20))},[transferDir filesep 'IBDDiffSubsystems.txt'],'\t',{'sub','diffabund'});
    writeData({addIdxStrings(sortedSubsystems(end-20:end)),diffAbunds2(sortIdxs(end-20:end))},[transferDir filesep 'IBDDiffSubsystemsBottom.txt'],'\t',{'sub','diffabund'});
end
