
if 0
minModel = makeMinimalFBAModel(bigModelTableFlux);
minModelTrue = subselectModel(minModel,'minModelTrue',minModel.rxns(minModel.lb~=0 | minModel.ub~=0));
[reducedModelArrBig subsystemsAddedArrBig] = makeReducedModelArr(minModelTrue,bigModelTableFlux,1);
end

if 0
falconfluxesArr = {};
fbafluxesArr = {};
efluxesArr = {};
for i=5:length(reducedModelArrBig)
    efluxesArr{i} = runFluxMethod(expressionData,expressionIDs,'testeflux',reducedModelArrBig{i},'EFlux',ones(length(expressionData),1),'BIOMASS');
    fbasol = optimizeCbModel(reducedModelArrBig{i});
    fbafluxesArr{i} = fbasol.x;
    falconfluxesArr{i} = runFluxMethod(expressionData,expressionIDs,'testfalcon',reducedModelArrBig{i},'FALCON',expressionSDs);
end
end

minModelSubsystems = unique(minModelTrue.subSystems);
minModelSubsystemsCounts = cellfun(@(x) sum(strcmp(minModelTrue.subSystems,x)), minModelSubsystems);
[minModelSubsystemsCounts sortIdxs] = sort(minModelSubsystemsCounts);
minModelSubsystems = addIdxStrings(minModelSubsystems(sortIdxs));
writeData({minModelSubsystems,minModelSubsystemsCounts},'/home/fs01/yw595/minModelSubsystems.txt','\t',{'subsystem','count'});

subsystemsAddedArrBigLabel = addIdxStrings(subsystemsAddedArrBig);
corrMatrix = [];
corrLabels = {'FALCON vs. FBA','FALCON vs. E-Flux','FBA vs. E-Flux'};
for i=5:length(falconfluxesArr)
    corrMatrix(i,1) = corr(falconfluxesArr{i},fbafluxesArr{i},'type','Spearman');
    corrMatrix(i,2) = corr(falconfluxesArr{i},efluxesArr{i},'type','Spearman');
    corrMatrix(i,3) = corr(fbafluxesArr{i},efluxesArr{i},'type','Spearman');
end
xvals = []; yvals = []; xlabels = {}; group = {};
for i=1:length(corrLabels)
    for j=5:length(falconfluxesArr)
        if ~isnan(corrMatrix(j,i))
            xvals(end+1) = j;
            yvals(end+1) = corrMatrix(j,i);
            xlabel = subsystemsAddedArrBigLabel{j};
            xlabel = xlabel(1:min(30,length(xlabel)));
            xlabels{end+1} = xlabel;
            group{end+1} = corrLabels{i};
        end
    end
end

writeData({xvals,yvals,xlabels,group},'/home/fs01/yw595/reducedModelCorrs.txt','\t',{'xvals','yvals','xlabels','group'});





