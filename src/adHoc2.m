metsToFluxDiff = zeros(length(bigModelTableFlux.mets),1);
for i=1:length(diffIdxs)
	diffSub = bigModelTableFlux.subSystems{diffIdxs(i)};
if strcmp(diffSub,'') || strcmp(diffSub,'Transport') || strcmp(diffSub,'Exchange') || strcmp(diffSub,'MUST_EX')
%disp(bigModelTableFlux.rxnNames{diffIdxs(i)});
metIdxs = find(bigModelTableFlux.S(:,diffIdxs(i)));
for j=1:length(metIdxs)
	metsToFluxDiff(metIdxs(j)) = metsToFluxDiff(metIdxs(j)) + abs(normFluxesNormalArr2(diffIdxs(i),1)-normFluxesNormalArr2(diffIdxs(i),3));
end
     end
     end

metsToFluxDiffSort = sort(metsToFluxDiff,'descend');
[~,sortIdxs] = sort(metsToFluxDiff);
metsToFluxDiffSort = metsToFluxDiffSort(10:end);
sortIdxs = sortIdxs(10:end);
writeData({addIdxStrings(bigModelTableFlux.metNames(sortIdxs)),metsToFluxDiffSort},'/mnt/vdb/home/ubuntu2/hostFluxDiff2.txt','\t',{'metname','fluxdiff'});

[subLabels,subFluxSubs,subFluxDiffSubs] = segmentFluxBySubsystem(bigModelTableFlux,normFluxesNormalArr2(:,1),1,normFluxesNormalArr2(:,3),1);
subFluxDiffSubs = abs(subFluxDiffSubs);
[subFluxDiffSubs,sortIdxs] = sort(subFluxDiffSubs,'descend');
subLabels = subLabels(sortIdxs);
subFluxDiffSubs = subFluxDiffSubs(4:end); subLabels = subLabels(4:end);
writeData({addIdxStrings(subLabels),subFluxDiffSubs},'/mnt/vdb/home/ubuntu2/hostFluxDiffSubs2.txt','\t',{'subname','fluxdiff'});
