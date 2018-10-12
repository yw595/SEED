function [subLabels,subFluxSums,subFluxDiffSums,subFluxStdSums,subFluxDiffScaledSums] = segmentFluxBySubsystem(model,fluxDist1,diffCond,fluxDist2,innerAbs,passedStds,fluxStds)

if ~exist('diffCond','var')
    diffCond=0;
end
if ~exist('passedStds','var')
    passedStds=0;
end
uniqSubs = unique(model.subSystems);

subLabels = {}; subFluxSums = []; subFluxDiffSums = []; subFluxStdSums = []; subFluxDiffScaledSums = [];
for j=1:length(uniqSubs)
    subLabels{end+1} = uniqSubs{j};
    if innerAbs==1
        subFluxSums(end+1) = sum(abs(fluxDist1(strcmp(model.subSystems,uniqSubs{j}))));
    else
        subFluxSums(end+1) = sum(fluxDist1(strcmp(model.subSystems,uniqSubs{j})));
    end
    if diffCond==1
        subFluxDiffSums(end+1) = abs(sum(fluxDist1(strcmp(model.subSystems,uniqSubs{j}))-fluxDist2(strcmp(model.subSystems,uniqSubs{j}))));
        if innerAbs
	    subFluxDiffSums(end) = sum(abs(fluxDist1(strcmp(model.subSystems,uniqSubs{j}))))-sum(abs(fluxDist2(strcmp(model.subSystems,uniqSubs{j}))));
        end
    end
    if passedStds==1
        subFluxStdSums(end+1) = abs(sum(fluxStds(strcmp(model.subSystems,uniqSubs{j}))));
    end
    if diffCond==1 && passedStds==1
        subFluxDiffScaledSums(end+1) = subFluxDiffSums(end)/subFluxStdSums(end);
    end
end
