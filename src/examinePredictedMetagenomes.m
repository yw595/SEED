if 1
KOsArr = {};
for i=1:length(ucrFolders)
    ucrFolder = ucrFolders{i};
    system(['python /mnt/xvdf/home/ubuntu2/pickPresent.py ' inputDir filesep 'MGMData' filesep ucrFolder filesep 'normalized_otus.tsv']);
    presentSpecies = textscan(fopen('/mnt/xvdf/home/ubuntu2/pickPresent.txt'),'%s','Delimiter','\n','HeaderLines',0);
    presentSpecies = presentSpecies{1};
    KOsToAbunds = containers.Map;

    inFI = fopen([inputDir filesep 'MGMData' filesep ucrFolders{i} filesep 'normalized_otus_predicted_metagenome.tsv']);
    KOsToAbunds2 = containers.Map;
    line = fgetl(inFI);
    line = fgetl(inFI);
    line = fgetl(inFI);
    while line~=-1
        words = strsplit(line);
        if ~isKey(KOsToAbunds2,words{1})
            KOsToAbunds2(words{1}) = str2num(words{2});
        else
	    KOsToAbunds2(words{1}) = KOsToAbunds2(words{1})+str2num(words{2});
        end
	line = fgetl(inFI);
    end
    fclose(inFI);
    for j=1:length(AGORAMat)
	j
	i
        for k=1:length(presentSpecies)
	    if ~isempty(regexp(AGORAMat{j},presentSpecies{k}))
	        %inFI = fopen([inputDir filesep 'MGMData' filesep ucrFolders{i} filesep 'normalized_otus_predicted_metagenome.tsv']);
                inFI = fopen([inputDir filesep 'MGMData' filesep ucrFolders{i} filesep 'two_otus_test' num2str(j) '_predicted_metagenome.tsv']);
	        line = fgetl(inFI);
                line = fgetl(inFI);
                line = fgetl(inFI);
                while line~=-1
                    words = strsplit(line);
                    if ~isKey(KOsToAbunds,words{1})
                        KOsToAbunds(words{1}) = str2num(words{2});
		    else
		        KOsToAbunds(words{1}) = KOsToAbunds(words{1})+str2num(words{2});
                    end
		    line = fgetl(inFI);
                end
		fclose(inFI);
            end
        end
    end
    KOsArr{i,1} = KOsToAbunds;
    KOsArr{i,2} = KOsToAbunds2;
end
end

if 0
	    KOKeys = keys(KOsToAbunds);
KO1 = [];
KO2 = [];
for i=1:length(KOKeys)
	KO1(end+1) = KOsToAbunds(KOKeys{i});
KO2(end+1) = KOsToAbunds2(KOKeys{i});
end
end

if 0
KOSub1 = [];
KOSub2 = [];
for i=1:12
	KOSub1(i) = sum(KO1((i-1)*500+1:i*500+1));
        KOSub2(i) = sum(KO2((i-1)*500+1:i*500+1));
end
writeData({KOSub1,KOSub2},'/mnt/xvdf/home/ubuntu2/totalOxExprCompare.txt','\t',{'expr1','expr2'});
end
