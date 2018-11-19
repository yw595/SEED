if 1
useDiabetes = 0;
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
end
end

sampleVec = {};
speciesVec = {};
abundVec = [];
for i=1:length(sampleNames)
    disp(i)
    for j=1:length(AGORAMat)
	sampleVec{end+1} = sampleNames{i};
	AGORAName = strsplit(AGORAMat{j},'/');
        AGORAName = strrep(AGORAName{end},'_',' ');
	%AGORAName = strsplit(AGORAName{end},'_');
	%AGORAName = [AGORAName{1} ' ' AGORAName{2}];
        speciesVec{end+1} = AGORAName;
        abundVec(end+1) = AGORAModelAbundMat(i,j);
    end
end
if useDiabetes
    writeData({sampleVec,speciesVec,abundVec},'/mnt/vdb/home/ubuntu2/diabetesVec.txt','\t',{'sample','species','abund'});
else
    writeData({sampleVec,speciesVec,abundVec},'/mnt/vdb/home/ubuntu2/IBDVec.txt','\t',{'sample','species','abund'});
end

yokflux = runFluxMethod(ones(length(yok.rxns),1),yok.rxns,'',yok,'FALCON',ones(length(yok.rxns),1));
compare1 = textscan(fopen('/mnt/vdb/home/ubuntu2/compareFALCONs.txt','%f'));
corr(yokflux,compare1,'type','Spearman')
compare1nonzero = find(compare1~=0);
yokfluxnonzero = find(yokflux~=0);
intersectnonzero = intersect(compare1nonzero,yokfluxnonzero);
corr(yokflux(intersectnonzero),compare1(intersectnonzero),'type','Spearman')
