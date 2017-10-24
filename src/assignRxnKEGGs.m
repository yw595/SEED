function assignRxnArr = assignRxnKEGGs(origModel,equations)

splitEqs = cellfun(@(x) strsplit(x,' '),equations,'UniformOutput',0);
  
assignRxnArr = {};
parfor i=1:length(origModel.rxnNames)
    i
    rxnName = origModel.rxnNames{i};
    if length(rxnName)>=3 && strcmp(rxnName(1:3),'rxn')
	%bigModelTableFlux.rxnKEGGs{i} = rxnName;
        assignRxnArr{i} = rxnName;
    else
	consumedMets = origModel.metKEGGs(origModel.S(:,i)<0);
	producedMets = origModel.metKEGGs(origModel.S(:,i)>0);
        consumedMets = cellfun(@(x) strrep(x,'C','cpd'), consumedMets,'UniformOutput',0);
        producedMets = cellfun(@(x) strrep(x,'C','cpd'), producedMets,'UniformOutput',0);
        allMets = union(consumedMets,producedMets);

        matchNumArr = cellfun(@(x) length(intersect(x,allMets))/length(x), splitEqs);
        [~,maxIdx]=max(matchNumArr);

        %disp(maxMatch)
        %if maxMatch>0
            assignRxnArr{i} = rxnData{maxIdx,1};
            if strcmp(assignRxnArr{i},'DATABASE')
                assignRxnArr{i} = '';
            end
            assignRxnArr{i}
        %end

	if 0
            plausibleEq = '';
            for j=1:length(consumedMets)-1
	        plausibleEq = [plausibleEq consumedMets{j} ' + '];
            end
            if length(consumedMets)>0
                plausibleEq = [plausibleEq consumedMets{end} ' <=> '];
            end
            for j=1:length(producedMets)-1
	        plausibleEq = [plausibleEq producedMets{j} ' + '];
            end
            if length(producedMets)>0
                plausibleEq = [plausibleEq producedMets{end}];
            end
            %plausibleEq
            matchPlausNum = sum(strcmp(equations,plausibleEq));
            if matchPlausNum>0
                plausibleEq
                matchPlausNum
            end
        end
    end
end
