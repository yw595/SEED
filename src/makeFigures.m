if 0
    AGORAModelArr = {};
    for i=1:length(AGORAMat)
	load(AGORAMat{i})
	i
	AGORAModelArr{i} = AGORAModel;
    end
end

if 0
sampFiles = {'/mnt/vdb/home/ubuntu2/MHnormal.txt','/mnt/vdb/home/ubuntu2/MHobese.txt'};
sampData = {};
for z1=1:2
    %system(['cat ' outputDir filesep 'simulateSmallModelsSeparatePicrustAGORA' filesep '*_1_* > ' outputDir filesep 'simulateSmallModelsSeparatePicrustAGORA' filesep 'AllFirst.txt']);
    a = importdata([outputDir filesep 'simulateSmallModelsSeparatePicrustAGORA' filesep 'AllFirst.txt']);
    samps = importdata(sampFiles{z1});
    allFiles = dir([outputDir filesep 'simulateSmallModelsSeparatePicrustAGORA']);
    currIdx = 0;
    for i=1:length(allFiles)
	i
	if ~isempty(regexp(allFiles(i).name,'_1_'))
	    words = strsplit(allFiles(i).name,'_');
            if strcmp(words{3},'1')
		ithField = str2num(words{1});
		jthField = str2num(words{2});
		sampField = words{4};
		sampField = strsplit(sampField,'.f');
		sampField = sampField{1};
		sampField = sampField(7:end);
		sampIdx = find(strcmp(samps,sampField));
		if ~isempty(sampIdx)
		    sampData{z1,sampIdx,ithField,jthField} = a(currIdx+1:currIdx+length(AGORAModelArr{ithField}.rxns));
		    currIdx = currIdx+length(AGORAModelArr{ithField}.rxns);
		end
	    end
        end
    end
end
end


if 0
    for z1=1:2
	ithTotalTermArrSort = [];
	compMets = containers.Map;
	coopMets = containers.Map;
	matrixSort = [];
	allMatrices = {};
	for z=1:size(sampData,2)
	    z
	    allMatrixZth = [];
	    for i=1:min(size(sampData,3),size(sampData,4))
		%i
		ithTotalTerm = 0;
		for j=1:min(size(sampData,3),size(sampData,4))
		    modelIth = AGORAModelArr{i};
		    modelJth = AGORAModelArr{j};
		    fluxesIth = sampData{z1,z,i,j};
		    fluxesJth = sampData{z1,z,j,i};
		    if ~isempty(fluxesIth) && ~isempty(fluxesJth)
			compTerm = 0;
			coopTerm = 0;
			for k=1:length(modelIth.rxns)
			    if ~isempty(regexp(modelIth.subSystems{k},'Exchange'))
				if ~isempty(strcmp(modelJth.rxns,modelIth.rxns{k}))
				    correspondK = find(strcmp(modelJth.rxns,modelIth.rxns{k}));
				    if fluxesIth(k) < 0 && mean(fluxesJth(correspondK)) < 0
					compTerm = compTerm + fluxesIth(k)*mean(fluxesJth(correspondK));
					if ~isKey(compMets,modelIth.rxns{k})
					    compMets(modelIth.rxns{k}) = fluxesIth(k)*mean(fluxesJth(correspondK));
					else
					    compMets(modelIth.rxns{k}) = compMets(modelIth.rxns{k}) + fluxesIth(k)*mean(fluxesJth(correspondK));
					end
				    end
				    if (fluxesIth(k) < 0 && mean(fluxesJth(correspondK)) > 0) || (fluxesIth(k) > 0 && mean(fluxesJth(correspondK)) < 0)
					coopTerm = coopTerm - fluxesIth(k)*mean(fluxesJth(correspondK));
					if ~isKey(coopMets,modelIth.rxns{k})
					    coopMets(modelIth.rxns{k}) = -fluxesIth(k)*mean(fluxesJth(correspondK));
					else
					    coopMets(modelIth.rxns{k}) = coopMets(modelIth.rxns{k}) - fluxesIth(k)*mean(fluxesJth(correspondK));
					end
				    end
				end
			    end
			end
			ithTotalTerm = ithTotalTerm - compTerm + coopTerm;
			matrixSort(i,j) = -compTerm + coopTerm;
			allMatrixZth(i,j) = -compTerm + coopTerm;
		    end
		end
		ithTotalTermArrSort(i) = ithTotalTerm;
	    end
	    allMatrices{z} = allMatrixZth;
	end
	if z1==1
	    allMatricesNormal = allMatrices;
        else
	    allMatricesObese = allMatrices;
	end
    end
end

groups = {};
obeseMatrixSums = [];
for i=1:length(allMatricesObese)
    obeseMatrixSums(i) = sum(sum(allMatricesObese{i}));
    groups{end+1} = 'obese';
end
normalMatrixSums = [];
for i=1:length(allMatricesNormal)
    normalMatrixSums(i) = sum(sum(allMatricesNormal{i}));
    groups{end+1} = 'normal';
end
ttest2(obeseMatrixSums,normalMatrixSums)
writeData({groups,[obeseMatrixSums normalMatrixSums]},[transferDir filesep 'samplecoops.txt'],'\t',{'group','coop'})
