rxnsDiff = keys(rxnsToDiffExp);
exRxns = {};
for i=1:length(rxnsDiff)
	rxnDiff = rxnsDiff{i};
matchedSub = bigModelTableFlux.subSystems(strcmp(bigModelTableFlux.rxns,rxnDiff));
if ~isempty(matchedSub)
matchedSub = matchedSub{1};
end
matchedSub
if strcmp(matchedSub,'')
exRxns{end+1} = rxnsDiff{i};
end
end

exSubs = {};
nonEmptySubs = sort(unique(bigModelTableFlux.subSystems));
nonEmptySubs = nonEmptySubs(2:end);
for i=1:length(exRxns)
	
