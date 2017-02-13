function rxnMatrix = makeRxnConnMatrix(model,excludeTopMets,percentThresh)

if ~exist('percentThresh','var')
    percentThresh=.001;
end

[~, topMetIdxsOrig] = sort(sum(abs(model.S)~=0,2),'descend');
topMetIdxs = topMetIdxsOrig(1:floor(length(topMetIdxsOrig)*percentThresh));
excludeMets = {'ATP','ADP','Phosphate','NAD','NADH','CoA','PPi','NADP','NADPH','O2','AMP','NH3','CO2','UDP','S-Adenosyl-L-methionine','cpd00001','cpd00002','cpd00003','cpd00004','cpd00005','cpd00006','cpd00007','cpd00008','cpd00009','cpd00010','cpd00011','cpd00012','cpd00013','cpd00014','cpd00015','cpd00016','cpd00017','cpd00018','cpd00019','cpd00020','cpd00021'};
[~,intersectIdxs,~] = intersect(model.metNames,excludeMets);
topMetIdxs = union(topMetIdxs,intersectIdxs);

rxnMatrix = zeros(length(model.rxns),length(model.rxns));
for i=1:length(model.rxns)
    metIdxs = find(model.S(:,i)~=0);
    for j=1:length(metIdxs)
        if excludeTopMets==0 || ~any(metIdxs(j)==topMetIdxs)
            connRxnIdxs = find(model.S(metIdxs(j),:));
            for k=1:length(connRxnIdxs)
                rxnMatrix(i,connRxnIdxs(k))=1;
            end
        end
    end
end