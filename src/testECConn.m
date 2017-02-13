for i=1:length(bigModelTableNoTE.rxns)
    if mod(i,100)==0
        disp(i)
    end
    for j=1:length(bigModelTableNoTE.rxns)
        sourceECNums = bigModelTableNoTE.rxnECNums{i};
        targetECNums = bigModelTableNoTE.rxnECNums{j};
        for k=1:length(sourceECNums)
            for l=1:length(targetECNums)
                if any(strcmp(sourceECNums,'6.4.1.1')) && any(strcmp(targetECNums,'2.7.1.69'))
                    sourceMetIdxs = find(bigModelTableNoTE.S(:,i)~=0);
                    targetMetIdxs = find(bigModelTableNoTE.S(:,j)~=0);
                    if ~isempty(intersect(sourceMetIdxs,targetMetIdxs))
                        disp(bigModelTableNoTE.metNames(sourceMetIdxs))
                        disp(bigModelTableNoTE.metNames(targetMetIdxs))
                        disp(bigModelTableNoTE.rxnNames{i});
                        disp(bigModelTableNoTE.rxnNames{j});
                    end
                end
            end
        end
    end
end