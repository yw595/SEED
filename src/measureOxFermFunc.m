function [totalOx totalFerment oxRxnIdxs fermentRxnIdxs] = measureOxFermFunc(model,pseudoFlux)
%nadhIdxs = cellfun(@(x) ~isempty(regexp(x,'Nicotinamide|NADH')), mergedModel.metNames);
%fadhIdxs = cellfun(@(x) ~isempty(regexp(x,'Nicotinamide|FADH')), mergedModel.metNames);

if ~isfield(model,'metKEGGs') && isfield(model,'metKEGGID')
    model.metKEGGs = model.metKEGGID;
end

fermentKEGGs = {'C00186','C00469','C00246','C00033','C00163','C00084'};
oxKEGGs = {'C00007'};

fermentRxnIdxs = [];
oxRxnIdxs = [];
for i=1:length(model.rxns)
    consumeMetKEGGs = model.metKEGGs(model.S(:,i)<0);
    produceMetKEGGs = model.metKEGGs(model.S(:,i)>0);
    if length(intersect(produceMetKEGGs,fermentKEGGs))>0 && length(intersect(consumeMetKEGGs,fermentKEGGs))==0
	fermentRxnIdxs(end+1) = i;
    end
    if length(intersect(produceMetKEGGs,oxKEGGs))==0 && length(intersect(consumeMetKEGGs,oxKEGGs))>0
	oxRxnIdxs(end+1) = i;
    end
end
    
if ~exist('pseudoFlux','var')
    totalOx = length(oxRxnIdxs);
    totalFerment = length(fermentRxnIdxs);
else
    totalOx = sum(abs(pseudoFlux(oxRxnIdxs)));
    totalFerment = sum(abs(pseudoFlux(fermentRxnIdxs)));
end
