function [bootstrapVals bootstrapLabels] = bootstrap(origVals,numIters,numToDrop,origLabels)

bootstrapVals = [];
bootstrapLabels = {};    
for i=1:numIters
    ithPermIdxs = randperm(length(origVals));
    ithPerm = origVals(ithPermIdxs);
    bootstrapVals(end+1) = sum(ithPerm(1:length(origVals)-numToDrop));
    if exist('origLabels','var')
        ithPerm = origLabels(ithPermIdxs);
        bootstrapLabels{end+1} = ithPerm(1:length(origLabels)-numToDrop);
    end
end