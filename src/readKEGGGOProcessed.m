if 0

keggFI = fopen('keggmodulecopyprocessed.txt');
dataFields = textscan(keggFI,'%s%s%s','Delimiter','\t','HeaderLines',0);
fclose(keggFI);
dataFields = [dataFields{:}];
keggmodules = dataFields(:,1);
keggreactions = dataFields(:,3);
keggrxnsToModules = containers.Map;

for i=1:length(keggmodules)
    keggrxns = strsplit(keggreactions{i});
    for j=1:length(keggrxns)
        if ~isKey(keggrxnsToModules,keggrxns{j})
            keggrxnsToModules(keggrxns{j}) = {};
        end
        currentArr = keggrxnsToModules(keggrxns{j});
        currentArr{end+1} = keggmodules{i};
        keggrxnsToModules(keggrxns{j}) = currentArr;
    end
end

goFI = fopen('go-basic.obo.processed');
dataFields = textscan(goFI,'%s%s%s','Delimiter','\t','HeaderLines',0);
fclose(goFI);
dataFields = [dataFields{:}];
golayers = dataFields(:,1);
gonames = dataFields(:,2);
goids = dataFields(:,3);
golayersToIDs = containers.Map;

for i=1:length(golayers)
    ithgoids = strsplit(goids{i});
    for j=1:length(ithgoids)
	if ~isKey(golayersToIDs,ithgoids{j})
            golayersToIDs(ithgoids{j}) = {};
        end
        currentArr = golayersToIDs(ithgoids{j});
        currentArr{end+1} = golayers{i};
        golayersToIDs(ithgoids{j}) = currentArr;
    end
end

end

normalVector = normFluxesNormalArr{1};
normalVector = normalVector(:,3);
obeseVector = normFluxesObeseArr{11};
obeseVector = obeseVector(:,3);
    
if 0
goFI = fopen('go-basic.obo.processed');
dataFields = textscan(goFI,'%s%s%s','Delimiter','\t','HeaderLines',0);
fclose(goFI);
dataFields = [dataFields{:}];
golayers = dataFields(:,1);
gonames = dataFields(:,2);
goids = dataFields(:,3);
tempGOs = bigModelTableFlux.rxnGOs;
for i=1:length(tempGOs)
    if isempty(tempGOs{i})
	tempGOs{i} = '';
    end
end
testgo = unique(tempGOs);
golayersToOccurs = [];%containers.Map;
golayersToLengths = [];
%golayers = keys(golayersToIDs);
for i=1:length(golayers)
    if mod(i,1000)==0
	disp(i);
    end
    golayersToLengths(i) = length(strsplit(goids{i},' '));
    if golayersToLengths(i) >= 3
        [~,intersectIdxs] = intersect(testgo,strsplit(goids{i},' '));
        golayersToOccurs(i) = sum(abs(normalVector(intersectIdxs)-obeseVector(intersectIdxs)));%length(intersect(testgo,strsplit(goids{i},' ')));% / length(strsplit(goids{i},' '));
        %disp(golayersToOccurs(i))
    else
        golayersToOccurs(i) = 0;
    end
end
end

keggFI = fopen('keggmodulecopyprocessed.txt');
dataFields = textscan(keggFI,'%s%s%s','Delimiter','\t','HeaderLines',0);
fclose(keggFI);
dataFields = [dataFields{:}];
keggmodules = dataFields(:,1);
keggnames = dataFields(:,2);
keggreactions = dataFields(:,3);
testkegg = unique(bigModelTableFlux.rxnKEGGs);
for i=1:length(testkegg)
    ithKegg = testkegg{i};
    if length(ithKegg) > 3 && strcmp(ithKegg(1:3),'rxn')
        ithKegg(1:3)='';
        ithKegg = ['K' ithKegg];
    end
    testkegg{i} = ithKegg;
end
keggmodulesToOccurs = [];%containers.Map;
for i=1:length(keggmodules)
    if mod(i,1000)==0
	disp(i);
    end
    [~,intersectIdxs] = intersect(testkegg,strsplit(keggreactions{i},' '));
    keggmodulesToOccurs(i) = sum(abs(normalVector(intersectIdxs)-obeseVector(intersectIdxs)));%length(intersect(testkegg,strsplit(keggreactions{i},' '))) / length(strsplit(keggreactions{i},' '));
end

[~,sortIdxs] = sort(keggmodulesToOccurs,'descend');
keggmodulesToOccurs = keggmodulesToOccurs(sortIdxs);
keggnames = keggnames(sortIdxs);
keggnames = addIdxStrings(keggnames);
writeData({keggnames,keggmodulesToOccurs},[transferDir filesep 'keggoccurs.txt'],'\t',{'keggmodule','keggoccurs'});
[~,sortIdxs] = sort(golayersToOccurs,'descend');
golayersToOccurs = golayersToOccurs(sortIdxs);
gonames = gonames(sortIdxs);
gonames = addIdxStrings(gonames);
writeData({gonames,golayersToOccurs},[transferDir filesep 'gooccurs.txt'],'\t',{'golayer','gooccurs'});
