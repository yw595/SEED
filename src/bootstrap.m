function bootstrapVals = bootstrap(origVals,numIters,numToDrop)

bootstrapVals = [];
for i=1:numIters
    ithPerm = origVals(randperm(length(origVals)));
    bootstrapVals(end+1) = sum(ithPerm(1:length(origVals)-numToDrop));
end

end