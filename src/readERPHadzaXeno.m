
expressionDataArr = {};
expressionSDsArr = {};
for z = 1:zLim
    expressionIDs = {};
    expressionData = [];
    expressionSDs = [];

    if useERP
        for i=1:length(ERPData{1})
            i
            ECNums = ERPData{1}{i};
            for k =1:length(bigModelTableFlux.rxnECNums)
                bigModelECNums = bigModelTableFlux.rxnECNums{k};
                if any(strcmp(ECNums,bigModelECNums))
                    %if ~any(strcmp(expressionIDs,bigModelTableFlux.rxns{k}))
                        expressionIDs{end+1} = bigModelTableFlux.rxns{k};
                        if z==1
                            expressionData(end+1) = str2num(ERPData{2}{i});
                        else
                            expressionData(end+1) = str2num(ERPData{3}{i});
                        end
                        expressionSDs(end+1) = 1;
                        %end
                end
            end
        end
    elseif useXeno
        for i=1:length(xenoData{1})
            i
            ECNums = xenoData{1}{i};
            for k =1:length(bigModelTableFlux.rxnECNums)
                bigModelECNums = bigModelTableFlux.rxnECNums{k};
                if any(strcmp(ECNums,bigModelECNums))
                    if ~any(strcmp(expressionIDs,bigModelTableFlux.rxns{k}))
                        expressionIDs{end+1} = bigModelTableFlux.rxns{k};
                        expressionData(end+1) = str2num(xenoData{2}{i});
                        expressionSDs(end+1) = str2num(xenoData{3}{i});
                    end
                end
            end
        end
    elseif useHadza
        for i=1:length(hadzaData{1})
            i
            ECNums = hadzaData{1}{i};
            for k =1:length(bigModelTableFlux.rxnECNums)
                bigModelECNums = bigModelTableFlux.rxnECNums{k};
                if any(strcmp(ECNums,bigModelECNums))
                    if ~any(strcmp(expressionIDs,bigModelTableFlux.rxns{k}))
                        expressionIDs{end+1} = bigModelTableFlux.rxns{k};
                        expressionData(end+1) = str2num(hadzaData{2}{i});
                        expressionSDs(end+1) = 1;
                    end
                end
            end
        end
    else
        if z==1
            expressionIDs = rxnsToExpressNormKeys;
        else
            expressionIDs = rxnsToExpressObeseKeys;
        end
        expressionData = [];
        for i=1:length(expressionIDs)
            if z==1
                expressionData(i) = rxnsToExpressNorm(expressionIDs{i});
            else
                expressionData(i) = rxnsToExpressObese(expressionIDs{i});
            end
        end
        expressionSDs = ones(size(expressionData));
    end
    z
    expressionDataArr{z} = expressionData;
    expressionSDsArr{z} = expressionSDs;
end
