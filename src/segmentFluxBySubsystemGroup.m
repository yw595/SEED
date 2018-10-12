function [subDiffArr subDiffArrNum subDiffMapIndividual rxnDiffArr rxnDiffArrNum rxnDiffMapIndividual] = segmentFluxBySubsystemGroup(NormalIdxs,IBDIdxs,picrustFluxesArr,tobemerged)

subDiffMapNormal = containers.Map;
subDiffMapIBD = containers.Map;
subDiffMapIndividual = containers.Map;
rxnDiffArr = {};
rxnDiffArrNum = [];
rxnDiffMapNormal = containers.Map;
rxnDiffMapIBD = containers.Map;
rxnDiffMapIndividual = containers.Map;
for z=1:length(picrustFluxesArr)
    for k=1:length(tobemerged.rxns)
	if ~isKey(rxnDiffMapIndividual,tobemerged.rxns{k})
	    rxnDiffMapIndividual(tobemerged.rxns{k}) = zeros(length(picrustFluxesArr),1);
        end
	temp = rxnDiffMapIndividual(tobemerged.rxns{k});
        temp(z) = picrustFluxesArr{z}(k);
        rxnDiffMapIndividual(tobemerged.rxns{k}) = temp;
	if sum(z==IBDIdxs)~=0
	    if ~isKey(rxnDiffMapIBD,tobemerged.rxns{k})
	        rxnDiffMapIBD(tobemerged.rxns{k}) = picrustFluxesArr{z}(k);
	    else
	        rxnDiffMapIBD(tobemerged.rxns{k}) = rxnDiffMapIBD(tobemerged.rxns{k}) + picrustFluxesArr{z}(k);
	    end
	else
	    if ~isKey(rxnDiffMapNormal,tobemerged.rxns{k})
	        rxnDiffMapNormal(tobemerged.rxns{k}) = picrustFluxesArr{z}(k);
	    else
	        rxnDiffMapNormal(tobemerged.rxns{k}) = rxnDiffMapNormal(tobemerged.rxns{k}) + picrustFluxesArr{z}(k);
	    end	  
	end
    end
	
    tobemerged2 = tobemerged;
    for k=1:length(tobemerged2.subSystems)
	if iscell(tobemerged2.subSystems{k})
	    tobemerged2.subSystems{k} = tobemerged2.subSystems{k}{1};
        end
    end
    [subLabelsArr,subFluxSums,~,~,~] = segmentFluxBySubsystem(tobemerged2,picrustFluxesArr{z},0,[],1);
    for k=1:length(subLabelsArr)
	if ~isKey(subDiffMapIndividual,subLabelsArr{k})
	    subDiffMapIndividual(subLabelsArr{k}) = zeros(length(picrustFluxesArr),1);
        end
	temp = subDiffMapIndividual(subLabelsArr{k});
        temp(z) = subFluxSums(k);
        subDiffMapIndividual(subLabelsArr{k}) = temp;
	if sum(z==IBDIdxs)~=0
	    if ~isKey(subDiffMapIBD,subLabelsArr{k})
		subDiffMapIBD(subLabelsArr{k}) = subFluxSums(k);
	    else
		subDiffMapIBD(subLabelsArr{k}) = subDiffMapIBD(subLabelsArr{k}) + subFluxSums(k);
	    end
	else
	    if ~isKey(subDiffMapNormal,subLabelsArr{k})
		subDiffMapNormal(subLabelsArr{k}) = subFluxSums(k);
	    else
		subDiffMapNormal(subLabelsArr{k}) = subDiffMapNormal(subLabelsArr{k}) + subFluxSums(k);
	    end	  
	end
    end
end
rxnKeys = keys(rxnDiffMapNormal);
rxnDiffArr = {};
rxnDiffArrNum = [];
for i=1:length(rxnKeys)
    rxnDiffArr{i} = rxnKeys{i};
    rxnDiffArrNum(i) = rxnDiffMapNormal(rxnKeys{i}) - rxnDiffMapIBD(rxnKeys{i});
end
subKeys = keys(subDiffMapNormal);
subDiffArr = {};
subDiffArrNum = [];
for i=1:length(subKeys)
    subDiffArr{i} = subKeys{i};
    subDiffArrNum(i) = subDiffMapNormal(subKeys{i}) - subDiffMapIBD(subKeys{i});
end
