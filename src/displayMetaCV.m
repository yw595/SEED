configSEED;
exchangeRxns = {};
allKeys = keys(metsToIDs);
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
            disp(i)
            for j=1:length(metIDs)
                for k=1:length(allKeys)
                    if strcmp(metsToIDs(allKeys{k}),metIDs{j})
                        disp(allKeys{k});
                    end
                end
            end
            disp(bigModelTableFlux.rxnNames{i})
            exchangeRxns{end+1} = bigModelTableFlux.rxnNames{i};
        end
        end
    end
end

if 1
useMatsumoto = 1;
hadzamets = {'Ile','C0','Met','Glu','Orn','Val','Lys','Ala','Leu','Pro','Gly','Arg','Spermidine','Putrescine','Gln','Phe','Asp','Cit','His','Asn','Tyr','Thr','Trp','Histamine','Taurine','Dopamine','Spermine'};
if useMatsumoto
    hadzamets = {'Methionine sulfoxide','Thiamine','Ile','Hypoxanthine','Choline','Carnitine','Lactic acid','Met','Glu','Uracil','Ornithine','Val','2-Oxoglutaric acid','Xanthine','Lys','Ala','Citric acid','Leu','Pro','Gly','Arg', 'Spermidine','Betaine','Gluconic acid','N-Acetylglucosamine','Putrescine','Gln','Malonic acid','N-Acetylmethionine','Quinic acid','Phe','Propionic acid','dTMP','Tyramine','Mannosamine','Cytidine','Asp','Xanthosine','Citrulline','Guanosine','His','3''-CMP; Cytidine-2;-monophosphate','Isethionic acid','Asn','Hexanoic acid','Guanine','Adenosine','Uridine','Tyr','CMP','Glucuronic acid','Fumaric acid','Homoserine','Succinic acid','N-Acetylglutamic acid','Taurocholic acid','Thr','Trp','Cytosine','Adenine','Taurine','Hypotaurine','Cadaverine','Uric acid','5-Oxoproline','5-Aminovaleric acid','Inosine','Glucosamine','2''Deoxyguanosine','Thymidine','2''-Deoxycytidine','N-Acetylneuraminic acid','3-Phenylpropionic acid','Trimethylamine N-oxide','Panthothenic acid','Nicotinamide','Pyridoxal','Nicotinic acid','Glycerophosphocholine','Glyceric acid','dCMP','Ala-Ala','dAMP','B-Ala','Cytsine','4-Methyl-2-oxopentanoic acid; 3-Methyl-2-oxovaleric acid','2-Oxoisovaleric acid','Spermine','Stachydrine','Pyridoxine','S-Adenosylmethionine'};
end
exchangeDiffs = [];
for i=1:4
    for j=1:length(bigModelTableFlux.rxns)
        matchIdx = find(strcmp(bigModelTableFlux.rxnNames{j},exchangeRxns));
        if ~isempty(matchIdx)
            if useMatsumoto
                exchangeDiffs(matchIdx,i) = convertArr2(1,j,i);
            else
                exchangeDiffs(matchIdx,i) = convertArr2(1,j,i) - convertArr3(1,j,i);
            end
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
if useMatsumoto
    writeData({metAndMethod,fluxChange,method},'/home/fs01/yw595/displayMetaCVMat.txt','\t',{'metAndMethod','fluxChange','method'});
else
    writeData({metAndMethod,fluxChange,method},'/home/fs01/yw595/displayMetaCV.txt','\t',{'metAndMethod','fluxChange','method'});
end
end










        