function [ECMatrix ECList] = makeECConnMatrix(model,excludeTopMets,ignoreDashes)

[~, topMetIdxsOrig] = sort(sum(abs(model.S)~=0,2),'descend');
topMetIdxs = topMetIdxsOrig(1:floor(length(topMetIdxsOrig)*.001));
excludeMets = {'H2O','H+','ATP','ADP','H2','H2CO3','Phosphate','NAD','NADH','H2O2','H2S','CoA','PPi','NADP','NADPH','O2','AMP','NH3','CO2','UDP','S-Adenosyl-L-methionine','Nitrate','Nitrite','NO','Sulfate','Sulfite','S-Adenosyl-homocysteine','cpd00001','cpd00067','cpd00002','cpd00008','cpd11640','cpd00242','cpd00009','cpd00003','cpd00004','cpd00025','cpd00239','cpd00010','cpd00012','cpd00006','cpd00005','cpd00007','cpd00018','cpd00013','cpd00011','cpd00014','cpd00017','cpd00209','cpd00075','cpd00418','cpd00048','cpd00081','cpd00019'};
[~,intersectIdxs,~] = intersect(model.metNames,excludeMets);
intersectIdxs
topMetIdxs = union(topMetIdxs,intersectIdxs);
model.metNames(topMetIdxs)

ECMatrix = [];
ECList = {};
count=0;
for i=1:length(model.rxns)
    if mod(i,100)==0
        disp(i)
        disp(length(ECList))
    end
    sourceECNums = model.rxnECNums{i};
    metIdxs = find(model.S(:,i)~=0);
    if i==757
        disp('HERE')
    end
    for j=1:length(metIdxs)
        if excludeTopMets==0 || ~any(metIdxs(j)==topMetIdxs)
            connRxnIdxs = find(model.S(metIdxs(j),:));
            for k=1:length(connRxnIdxs)
                rxnMatrix(i,connRxnIdxs(k))=1;
                targetECNums = model.rxnECNums{connRxnIdxs(k)};
                for l=1:length(sourceECNums)
                    for m=1:length(targetECNums)
                        if (isempty(regexp(sourceECNums{l},'-')) && isempty(regexp(targetECNums{m},'-'))) || ~ignoreDashes
                            lthIdx = strcmp(ECList,sourceECNums{l});
                            if ~any(lthIdx)
                                ECList{end+1} = sourceECNums{l};
                            end
                            lthIdx = strcmp(ECList,sourceECNums{l});
                            mthIdx = strcmp(ECList,targetECNums{m});
                            if ~any(mthIdx)
                                ECList{end+1} = targetECNums{m};
                            end
                            mthIdx = strcmp(ECList,targetECNums{m});
                            if count<100
                                count=count+1;
                                %disp(sourceECNums{l})
                                %disp(targetECNums{m})
                                %disp(lthIdx)
                                %disp(mthIdx)
                            end
                            ECMatrix(lthIdx,mthIdx)=1;
                        end
                    end
                end
            end
        end
    end
end