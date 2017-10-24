inputFI1 = fopen([inputDir filesep 'ShoaieRefInfo' filesep 'GOmappings' filesep 'kegg2go']);
dataFields = textscan(inputFI1,'%s%s','Delimiter',{' > '}, 'HeaderLines',2);
fclose(inputFI1);
dataFields = [dataFields{:}];
keggIDs = cellfun(@(x) strsplitYiping(x,':'), dataFields(:,1),'UniformOutput',0); keggIDs = cellfun(@(x) x{2}, keggIDs, 'UniformOutput',0);
goTerms = cellfun(@(x) strsplitYiping(x,' ; '),dataFields(:,2),'UniformOutput',0); goTerms = cellfun(@(x) x{2}, goTerms, 'UniformOutput',0);
