bigModelAddedTemp = bigModelAdded;
uniqueSubs = unique(bigModelAddedTemp.subSystems);
subsystemConnMatrix = zeros(length(uniqueSubs),length(uniqueSubs));
for i=1:length(bigModelAddedTemp.rxns)
    i
    ithSub = bigModelAddedTemp.subSystems{i};
    ithMetsIdxs = find(bigModelAddedTemp.S(:,i)~=0);
    for j=1:length(ithMetsIdxs)
        jthSubs = bigModelAddedTemp.subSystems(bigModelAddedTemp.S(ithMetsIdxs(j),:)~=0);
        for k=1:length(jthSubs)
            subsystemConnMatrix(strcmp(uniqueSubs,ithSub), strcmp(uniqueSubs,jthSubs{k})) = subsystemConnMatrix(strcmp(uniqueSubs,ithSub), strcmp(uniqueSubs,jthSubs{k}))+1;
        end
    end
end

subsystemConnMatrixNorm = subsystemConnMatrix;
for i=1:length(uniqueSubs)
    subsystemConnMatrixNorm(i,:) = subsystemConnMatrixNorm(i,:)/sum(subsystemConnMatrixNorm(i,:));
end