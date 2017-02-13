exchangeRxns = {};
for i=1:length(bigModelTableFlux.rxns)
    if strcmp(bigModelTableFlux.subSystems{i},'Exchange')
        % if any(strcmp(bigModelTableFlux.rxnNames{i}(4:end-2),keys(metsToIDs)))
        %     disp(bigModelTableFlux.rxnNames{i})
        %     disp(convertArr(1,i,4))
        % end
        metIDs = bigModelTableFlux.metKEGGs(bigModelTableFlux.S(:,i)~=0);
        savemets = metIDs;
        if ~any(strcmp(metIDs,''))
        metIDs = cellfun(@(x) [ x(2:6)], metIDs, 'UniformOutput',0);
        if length(intersect(metIDs,values(metsToIDs))) == length(metIDs)
            %disp(savemets)
            %disp(metIDs)
            disp(bigModelTableFlux.rxnNames{i})
            exchangeRxns{end+1} = bigModelTableFlux.rxnNames{i};
        end
        end
    end
end

hadzamets = {'Ile','Met','Glu','Orn','Val','Lys','Ala','Leu','Pro','Gly','Arg','Spermidine','Putrescine','Gln','Phe','Asp','Cit','His','Asn','Tyr','Thr','Trp','Histamine','Taurine','Dopamine','Spermine'};
exchangeDiffs = [];
for i=1:4
    for j=1:length(bigModelTableFlux.rxns)
        matchIdx = find(strcmp(bigModelTableFlux.rxnNames{j},exchangeRxns));
        if ~isempty(matchIdx)
            exchangeDiffs(matchIdx,i) = convertArr2(1,j,i) - convertArr3(1,j,i);
        end
    end
end
hadzadiffs = [];
for i=1:length(hadzamets)
    metDiffs = metsToAbunds(hadzamets{i});
    hadzadiffs(i) = metDiffs(1)-metDiffs(2);
end

%xvals = 1:5*length(hadzadiffs);
metAndMethod = {}; fluxChange = []; method = {};
for i=1:length(exchangeRxns)
    for j=1:4
        fluxChange(end+1) = exchangeDiffs(i,j);
        metAndMethod{end+1} = [exchangeRxns{i} ' ' methodsList{j}];
        method{end+1} = methodsList{j};
    end
    fluxChange(end+1) = hadzadiffs(i);
    metAndMethod{end+1} = [exchangeRxns{i} ' Experimental'];
    method{end+1} = 'Experimental';
end
writeData({metAndMethod,fluxChange,method},'/home/fs01/yw595/displayMetaCV.txt','\t',{'metAndMethod','fluxChange','method'});
            






        