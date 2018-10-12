function [returnrxns, returnfluxes] = traceFlux(model,fluxdist,specificrxn,excluderxns)

if ~exist('excluderxns','var')
    excluderxns = {};
end

posflux = 1;
if fluxdist(strcmp(model.rxns,specificrxn))<0
    posflux = 0;
end
metsToOccurs = containers.Map;
for i=1:length(model.mets)
    specificsum = sum(model.S(i,:)~=0);
    if posflux
	if model.S(i,strcmp(model.rxns,specificrxn))>0
	    metsToOccurs(model.mets{i}) = specificsum;
        end
    else
	if model.S(i,strcmp(model.rxns,specificrxn))<0
            metsToOccurs(model.mets{i}) = specificsum;
        end
    end	  
end
specificmets = keys(metsToOccurs);%model.mets(model.S(:,strcmp(model.rxns,specificrxn))~=0);
minmet = '';
minmetnum = length(model.rxns);
for i=1:length(specificmets)
    if metsToOccurs(specificmets{i})<=minmetnum && (~strcmp(specificmets{i},'akg[c]') || (sum(strcmp(specificmets,'ile_L[c]'))==0 && sum(strcmp(specificmets,'leu_L[c]'))==0))
	minmetnum = metsToOccurs(specificmets{i});
        minmet = specificmets{i};
    end
end
lessthanidxs = find(model.S(strcmp(model.mets,minmet),:)<0);
greaterthanidxs = find(model.S(strcmp(model.mets,minmet),:)>0);
posfluxidxs = find(fluxdist>0);
negfluxidxs = find(fluxdist<0);
goodidxs = union(intersect(lessthanidxs,posfluxidxs),intersect(greaterthanidxs,negfluxidxs));
specificfluxes = fluxdist(goodidxs);
specificrxns = model.rxns(goodidxs);
specificfluxes = specificfluxes(~strcmp(specificrxns,specificrxn));
specificrxns = specificrxns(~strcmp(specificrxns,specificrxn));
[~,diffIdxs] = setdiff(specificrxns,excluderxns);
specificfluxes = specificfluxes(diffIdxs);
specificrxns = specificrxns(diffIdxs);
[~,sortIdxs] = sort(abs(specificfluxes),'descend');
returnrxns = specificrxns(sortIdxs);
returnfluxes = specificfluxes(sortIdxs);
%[~,maxIdx] = max(abs(specificfluxes));
%returnrxn = specificrxns{maxIdx};
