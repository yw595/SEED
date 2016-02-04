configSEED;
bigModelTableTest = addMustEx(bigModelTable);
percentageProd = 0;
minRxnsAdded = 0;
while percentageProd < 1
    sols = testBiomass(bigModelTableTest,biomassMets,biomassCoeffs,0);
    producesAll = cellfun(@(x) x.f,sols);
    percentageProd = sum(producesAll > 0)/length(producesAll);
    disp(percentageProd)
    disp(minRxnsAdded)
    if percentageProd < 1
        sortNum = 1;
        [bigModelTableTestCand, minRxnsAddedCand, addedModel] = mergeSmallest(bigModelTableTest,modelNamesToModelsValues,sortNum);
        
        metToFind = biomassMets{find(producesAll==0,1)}; hasMet = 0;
        while ~hasMet
            biomassIdxs = find(cellfun(@(x) ~isempty(regexpi(x,'(biomass|grow)')),addedModel.rxnNames));
            if ~isempty(biomassIdxs)
                for j=1:length(biomassIdxs)
                    biomassMetIdxs = find(addedModel.S(:,biomassIdxs(j))~=0);
                    for k=1:length(biomassMetIdxs)
                        if any(strcmp(addedModel.mets{biomassMetIdxs(k)},metToFind))
                            hasMet=1;
                        end
                    end
                end
            end
            if ~hasMet
                sortNum=sortNum+1;
                [bigModelTableTestCand, minRxnsAddedCand, addedModel] = mergeSmallest(bigModelTableTest,modelNamesToModelsValues,sortNum);
            end        
        end
        
        percentageProdCand = percentageProd;
        while percentageProdCand==percentageProd
            solsCand = testBiomass(bigModelTableTestCand,biomassMets,biomassCoeffs,0);
            producesAllCand = cellfun(@(x) x.f,solsCand);
            percentageProdCand = sum(producesAllCand > 0)/length(producesAllCand);
            if percentageProdCand==percentageProd
                sortNum = sortNum+1;
                [bigModelTableTestCand, minRxnsAddedCand] = mergeSmallest(bigModelTableTest,modelNamesToModelsValues,sortNum);
            end
        end
        bigModelTableTest = bigModelTableTestCand;
        minRxnsAdded = minRxnsAddedCand;
    end
end