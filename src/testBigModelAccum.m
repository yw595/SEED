
if 0
bigModelNoTE = makeEmptyModel();
for i=1:length(bigModelAccum.subSystems)
    if ~strcmp(bigModelAccum.subSystems{i},'Transport') && ~strcmp(bigModelAccum.subSystems{i},'Exchange')
        i
        bigModelNoTE = mergeModels(bigModelNoTE,bigModelAccum,bigModelAccum.rxns{i});
        bigModelNoTE.rxnECNums{strcmp(bigModelNoTE.rxns,bigModelAccum.rxns{i})} = rxnsToECsAccum(bigModelAccum.rxns{i});
        bigModelNoTE = checkModelDims(bigModelNoTE);
    end
end

end

if 0
bigModelNoTEUnique = makeEmptyModel();
seenECs = {};
for i=1:length(bigModelNoTE.subSystems)
    currentECs = bigModelNoTE.rxnECNums{i};
    notInSeenECs = 0;
    for j=1:length(currentECs)
        if ~any(strcmp(currentECs{j},seenECs))
            notInSeenECs = 1;
            seenECs{end+1} = currentECs{j};
        end
    end
    if notInSeenECs
        i
        bigModelNoTEUnique = mergeModels(bigModelNoTEUnique,bigModelNoTE,bigModelNoTE.rxns{i});
        bigModelNoTEUnique.rxnECNums{strcmp(bigModelNoTEUnique.rxns,bigModelNoTE.rxns{i})} = rxnsToECsAccum(bigModelNoTE.rxns{i});
        bigModelNoTEUnique = checkModelDims(bigModelNoTEUnique);
    end
end
end

if 0
[ecnumsAccum compsAccum centAccum] = compsAndCent(bigModelAccum);
[ecnumsNoTE compsNoTE centNoTE] = compsAndCent(bigModelNoTE);
[ecnumsNoTEUnique compsNoTEUnique centNoTEUnique] = compsAndCent(bigModelNoTEUnique);
writeData({ecnumsAccum,centAccum},[baseDir filesep 'ECsToCentsAccum.txt'],'\t');
writeData({ecnumsNoTE,centNoTE},[baseDir filesep 'ECsToCentsNoTE.txt'],'\t');
writeData({ecnumsNoTEUnique,centNoTEUnique},[baseDir filesep 'ECsToCentsNoTEUnique.txt'],'\t');
end

if 1
[ecnumsTable compsTable centTable] = compsAndCent(bigModelTable);
[ecnumsTableNoTE compsTableNoTE centTableNoTE] = compsAndCent(bigModelTableNoTE);
[ecnumsTableNoTEUnique compsTableNoTEUnique centTableNoTEUnique] = compsAndCent(bigModelTableNoTEUnique);
writeData({ecnumsAccum,centAccum},[baseDir filesep 'ECsToCentsTable.txt'],'\t');
writeData({ecnumsTableNoTE,centTableNoTE},[baseDir filesep 'ECsToCentsTableNoTE.txt'],'\t');
writeData({ecnumsTableNoTEUnique,centTableNoTEUnique},[baseDir filesep 'ECsToCentsTableNoTEUnique.txt'],'\t');
end

if 0
ECsToRxnsTableKeys = keys(ECsToRxnsTable);
xvals = 1:length(ECsToRxnsTableKeys);
yvals = [];
for i=1:length(ECsToRxnsTableKeys)
    yvals(i) = length(ECsToRxnsTable(ECsToRxnsTableKeys{i}));
end
outputFile = [baseDir filesep 'ECsToRxnsTableDist.txt'];
writeForGGPlot(xvals,yvals,outputFile,ECsToRxnsTableKeys);

rxnsToECsTableKeys = keys(rxnsToECsTable);
xvals = 1:length(rxnsToECsTableKeys);
yvals = [];
for i=1:length(rxnsToECsTableKeys)
    yvals(i) = length(rxnsToECsTable(rxnsToECsTableKeys{i}));
end
outputFile = [baseDir filesep 'rxnsToECsTableDist.txt'];
writeForGGPlot(xvals,yvals,outputFile,rxnsToECsTableKeys);
end