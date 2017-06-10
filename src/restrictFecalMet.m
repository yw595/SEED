function returnModel = restrictFecalMet(origModel,metabolomeData)

fecalMet = metabolomeData{2};
fecalMet = fecalMet(~strcmp(fecalMet,''));
fecalMet = cellfun(@(x) ['C' x], fecalMet, 'UniformOutput', 0);
returnModel = origModel;
for i=1:length(returnModel.rxns)
    if strcmp(returnModel.subSystems{i},'Exchange')
        metIDs = returnModel.metKEGGs(returnModel.S(:,i)~=0);
        if length(intersect(metIDs,fecalMet)) ~= length(metIDs)
            %disp(i)
            returnModel.lb(i) = 0;
            returnModel.ub(i) = 0;
        end
    end
end