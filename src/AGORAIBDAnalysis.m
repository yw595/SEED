if 0
IBDDataMat = {};%zeros(length(AGORAMat),length(AGORAMat),1);
sampleArr = {};
subDiffMapArr = {};
for i=1:length(AGORAModelArr)
    subDiffMapArr{i,1} = containers.Map;
    subDiffMapArr{i,2} = containers.Map;
end
speciesSubDiffArr = {};
outputDir1 = [outputDir filesep 'simulateSmallModelsSeparatePicrustAGORAIBD'];
files = dir(outputDir1);
for i=1:2000%length(files)
    file = files(i).name;
    i
    if ~strcmp(file,'.') && ~strcmp(file,'..') && strcmp(file(end-3:end),'flux') && ~isempty(regexp(file,'_2_'))
        fileStripped = file(1:end-5);
        words = strsplit(fileStripped,'_');
        ithTerm = str2num(words{1});
        jthTerm = str2num(words{2});
        sample = words{4};
        if ~any(strcmp(sample,sampleArr))
            sampleArr{end+1} = sample;
        end
	ithFileData = importdata([outputDir1 filesep file]);
        if strcmp(words{3},'1')
            AGORAModel = AGORAModelArr{ithTerm};
        else
            AGORAModel = AGORAModelArr{jthTerm};
        end
	for k=1:length(AGORAModel.subSystems)
	    AGORAModel.subSystems{k} = AGORAModel.subSystems{k}{1};
        end
        [subLabelsArr,subFluxSums,~,~,~] = segmentFluxBySubsystem(AGORAModel,ithFileData,0,[],1,0,[]);%,1,picrustFluxesArr{z+15},1,0);

    if exist('sampleNames','var')
        subDiffMapNormal = subDiffMapArr{ithTerm,1};
        subDiffMapIBD = subDiffMapArr{ithTerm,2};
        for k=1:length(subLabelsArr)
	    if sum(strcmp(sample,sampleNames(1:2:end)))~=0
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
	subDiffMapArr{ithTerm,1} = subDiffMapNormal;
        subDiffMapArr{ithTerm,2} = subDiffMapIBD;

        speciesSubDiffArr{ithTerm,1} = subLabelsArr;
        speciesSubDiffArr{ithTerm,2} = subFluxSums;
        IBDDataMat{ithTerm,jthTerm,strcmp(sampleArr,sample)} = ithFileData;
        IBDDataMat{jthTerm,ithTerm,strcmp(sampleArr,sample)} = ithFileData;
    end
end
end
end

    %nonsense = nonsense+1;
    
if 1
allInteractMat = [];
for k=1:10%length(sampleArr)
    k
    for i=1:size(IBDDataMat,1)
	firstIdx = 0;
        disp(i)
	for j=1:size(IBDDataMat,2)
	    if ~isempty(IBDDataMat{i,j,k})% && ~isempty(IBDDataMat{j,i,k})
		disp(j)
		if firstIdx==0
		    firstIdx = j;
                end
	        [compTerm coopTerm coopBasicTerm compBasicTerm] = metInteractFuncTemp2(AGORAModelArr{firstIdx},AGORAModelArr{j},IBDDataMat{i,firstIdx,k},IBDDataMat{i,j,k});
 	        allInteractMat(k,i,j,1) = compTerm;
 	        allInteractMat(k,i,j,2) = coopTerm;
 	        allInteractMat(k,i,j,3) = coopBasicTerm;
 	        allInteractMat(k,i,j,4) = compBasicTerm;
            end
        end
    end
end
end

	    nonsense = nonsense+1;
	    
%IBDIdxs = [1,2,3];
%NormalIdxs = [4,5];

useIBD = 0;
useEFlux = 0;
if useIBD
    [~,NormalIdxs,~] = intersect(sampleNames,HMP2Normal);
    [~,IBDIdxs,~] = intersect(sampleNames,HMP2IBD);
    NormalIdxs = NormalIdxs(randperm(length(NormalIdxs),24));
    %IBDIdxs = IBDIdxs(randperm(length(IBDIdxs),93))
    %NormalIdxs = 1:2:length(picrustFluxesArr);
    %IBDIdxs = 2:2:length(picrustFluxesArr);
    picrustFluxesArr = picrustFluxesArrIBD;
    %subDiffKthMap = subDiffKthMapIBD;
    suffix = 'IBD';
else
    [~,IBDIdxs,~] = intersect(ucrFoldersShort,MHobese);
    [~,NormalIdxs,~] = intersect(ucrFoldersShort,MHnormal);
    picrustFluxesArr = picrustFluxesArrNormObese;
    if useEFlux
        picrustFluxesArr = picrustFluxesArrNormObeseEFlux;
    end
    %subDiffKthMap = subDiffKthMapNormObese;
    suffix = 'NormObese';
    if useEFlux
        suffix = 'NormObeseEFlux';
    end
    acArr = {};
    for i=1:length(ucrFolders)
        acArr{i} = importdata([inputDir '/MGMData/' ucrFolders{i} '/normalized_otus.modelexpr']);
    end
end
if 1
if 1
[subDiffArr subDiffArrNum subDiffMapIndividual rxnDiffArr rxnDiffArrNum rxnDiffMapIndividual] = segmentFluxBySubsystemGroup(NormalIdxs,IBDIdxs,picrustFluxesArr,tobemerged);
end
if useIBD
    for z=1:length(subDiffArr)
        subvals = subDiffMapIndividual(subDiffArr{z});
        subvals = subvals*3;
        subDiffMapIndividual(subDiffArr{z}) = subvals;
    end
end
if 1
[expSubDiffArr expSubDiffArrNum expSubDiffMapIndividual expRxnDiffArr expRxnDiffArrNum expRxnDiffMapIndividual] = segmentFluxBySubsystemGroup(NormalIdxs,IBDIdxs,acArr,tobemerged);
end
randIters = 200;
if 0
    for k=1:randIters
	k
	randidxs1 = randperm(length(picrustFluxesArr),length(NormalIdxs));
	randidxs2 = [];
	for k1=1:length(picrustFluxesArr)
	    if sum(randidxs1==k1)==0
		randidxs2(end+1) = k1;
	    end
	end
	[expSubDiffArrKth expSubDiffArrNumKth expSubDiffMapIndividualKth expRxnDiffArrKth expRxnDiffArrNumKth expRxnDiffMapIndividualKth] = segmentFluxBySubsystemGroup(randidxs1,randidxs2,acArr,tobemerged);
	for k1=1:length(expRxnDiffArrKth)
	    if ~isKey(expRxnDiffKthMap,expRxnDiffArrKth{k1})
		 expRxnDiffKthMap(expRxnDiffArrKth{k1}) = [expRxnDiffArrNumKth(k1)];
	    else
		temp = expRxnDiffKthMap(expRxnDiffArrKth{k1});
		temp(end+1) = expRxnDiffArrNumKth(k1);
                expRxnDiffKthMap(expRxnDiffArrKth{k1}) = temp;
	    end
	end
    end
end
if 0
    subDiffKthMap = containers.Map;
    rxnDiffKthMap = containers.Map;
    for k=1:randIters
	k
	randidxs1 = randperm(length(picrustFluxesArr),length(NormalIdxs));
	randidxs2 = [];
	for k1=1:length(picrustFluxesArr)
	    if sum(randidxs1==k1)==0
		randidxs2(end+1) = k1;
	    end
	end
	[subDiffArrKth subDiffArrNumKth subDiffMapIndividualKth rxnDiffArrKth rxnDiffArrNumKth rxnDiffMapIndividualKth] = segmentFluxBySubsystemGroup(randidxs1,randidxs2,picrustFluxesArr,tobemerged);
        for k1=1:length(rxnDiffArrKth)
	    k1thfluxes = rxnDiffMapIndividualKth(rxnDiffArrKth{k1});
            metric = varianceMetric(k1thfluxes(randidxs1),k1thfluxes(randidxs2));
	    if ~isKey(rxnMetricKthMap,rxnDiffArrKth{k1})
		 rxnMetricKthMap(rxnDiffArrKth{k1}) = [metric];
	    else
		temp = rxnMetricKthMap(rxnDiffArrKth{k1});
		temp(end+1) = metric;
                rxnMetricKthMap(rxnDiffArrKth{k1}) = temp;
	    end
	end
	for k1=1:length(rxnDiffArrKth)
	    if ~isKey(rxnDiffKthMap,rxnDiffArrKth{k1})
		 rxnDiffKthMap(rxnDiffArrKth{k1}) = [rxnDiffArrNumKth(k1)];
	    else
		temp = rxnDiffKthMap(rxnDiffArrKth{k1});
		temp(end+1) = rxnDiffArrNumKth(k1);
                rxnDiffKthMap(rxnDiffArrKth{k1}) = temp;
	    end
	end
	for k1=1:length(subDiffArrKth)
	    if ~isKey(subDiffKthMap,subDiffArrKth{k1})
		 subDiffKthMap(subDiffArrKth{k1}) = [subDiffArrNumKth(k1)];
	    else
		temp = subDiffKthMap(subDiffArrKth{k1});
		temp(end+1) = subDiffArrNumKth(k1);
                subDiffKthMap(subDiffArrKth{k1}) = temp;
	    end
	end
    end
end
for k1=1:length(rxnDiffArr)
    rxnDiffKthMap(rxnDiffArr{k1}) = sort(rxnDiffKthMap(rxnDiffArr{k1}),'descend');
end
for k1=1:length(rxnDiffArr)
    rxnMetricKthMap(rxnDiffArr{k1}) = sort(rxnMetricKthMap(rxnDiffArr{k1}),'ascend');
end
for k1=1:length(subDiffArr)
    subDiffKthMap(subDiffArr{k1}) = sort(subDiffKthMap(subDiffArr{k1}),'descend');
end

topSubs = {};
bottomSubs = {};
topPsSub = [];
bottomPsSub = [];
topRxns = {};
bottomRxns = {};
topPsRxn = [];
bottomPsRxn = [];
topRxnsT = {};
topPsRxnT = [];
topRxnsMetric = {};
topNumsRxnMetric = [];
topPsRxnMetric = [];
topPsRxnMetricNames = {};
topExpRxns = {};
bottomExpRxns = {};
topPsExpRxn = [];
bottomPsExpRxn = [];
topExpRxnsT = {};
topPsExpRxnT = [];
topExpRxnsMetric = {};
topNumsExpRxnMetric = [];
for k1=1:length(expRxnDiffArr)
    individualMap = expRxnDiffMapIndividual(expRxnDiffArr{k1});
    [h p] = ttest2(individualMap(NormalIdxs),individualMap(IBDIdxs));
    topExpRxnsT{k1} = expRxnDiffArr{k1};
    topPsExpRxnT(k1) = p;
    normalexps1 = individualMap(NormalIdxs);
    obeseexps1 = individualMap(IBDIdxs);
    metric = varianceMetric(normalexps1,obeseexps1);
    topExpRxnsMetric{k1} = expRxnDiffArr{k1};
    topNumsExpRxnMetric(k1) = metric;
    expRxnDiffArrNumKth = expRxnDiffArrNum(k1);
    expRxnDiffKthMapVals = expRxnDiffKthMap(expRxnDiffArr{k1});
    top = max(find(expRxnDiffKthMapVals >= expRxnDiffArrNumKth));
    if isempty(top)
        top = 0;
    end
    bottom = min(find(expRxnDiffKthMapVals <= expRxnDiffArrNumKth));
    if isempty(bottom)
        bottom = length(expRxnDiffKthMapVals)+1;
    end
    if 1%bottom-top==1
        topExpRxns{end+1} = expRxnDiffArr{k1};
        topPsExpRxn(end+1) = (top+1)/randIters;
        bottomExpRxns{end+1} = expRxnDiffArr{k1};
        bottomPsExpRxn(end+1) = (randIters-(bottom-1))/randIters;
    else
	disp(expRxnDiffKthMapVals)
	disp(bottom)
	disp(top)
    end
end
for k1=1:length(rxnDiffArr)
    individualMap = rxnDiffMapIndividual(rxnDiffArr{k1});
    [h p] = ttest2(individualMap(NormalIdxs),individualMap(IBDIdxs));
    topRxnsT{k1} = rxnDiffArr{k1};
    topPsRxnT(k1) = p;
    normalfluxes1 = individualMap(NormalIdxs);
    obesefluxes1 = individualMap(IBDIdxs);
    metric = varianceMetric(normalfluxes1,obesefluxes1);
    topRxnsMetric{k1} = rxnDiffArr{k1};
    topNumsRxnMetric(k1) = metric;
    rxnDiffArrNumKth = rxnDiffArrNum(k1);
    rxnDiffKthMapVals = rxnDiffKthMap(rxnDiffArr{k1});
    top = max(find(rxnDiffKthMapVals >= rxnDiffArrNumKth));
    if isempty(top)
        top = 0;
    end
    bottom = min(find(rxnDiffKthMapVals <= rxnDiffArrNumKth));
    if isempty(bottom)
        bottom = length(rxnDiffKthMapVals)+1;
    end
    if 1%bottom-top==1
        topRxns{end+1} = rxnDiffArr{k1};
        topPsRxn(end+1) = (top+1)/randIters;
        bottomRxns{end+1} = rxnDiffArr{k1};
        bottomPsRxn(end+1) = (randIters-(bottom-1))/randIters;
    else
	disp(rxnDiffKthMapVals)
	disp(bottom)
	disp(top)
    end

    rxnMetricKthMapVals = rxnMetricKthMap(rxnDiffArr{k1});
    rxnMetricKthMapVals = rxnMetricKthMapVals(~isnan(rxnMetricKthMapVals));
    if ~isempty(rxnMetricKthMapVals)
        top = max(find(rxnMetricKthMapVals <= metric));
        if isempty(top)
            top = 0;
        end
        topPsRxnMetric(end+1) = (top+1)/randIters;
        topPsRxnMetricNames{end+1} = rxnDiffArr{k1};
    else
	topPsRxnMetric(end+1) = 1;
        topPsRxnMetricNames{end+1} = rxnDiffArr{k1};
    end
end
for k1=1:length(subDiffArr)
    subDiffArrNumKth = subDiffArrNum(k1);
    subDiffKthMapVals = subDiffKthMap(subDiffArr{k1});
    top = max(find(subDiffKthMapVals >= subDiffArrNumKth));
    if isempty(top)
        top = 0;
    end
    bottom = min(find(subDiffKthMapVals <= subDiffArrNumKth));
    if isempty(bottom)
        bottom = length(subDiffKthMapVals)+1;
    end
    if bottom-top==1
        topSubs{end+1} = subDiffArr{k1};
        topPsSub(end+1) = (top+1)/randIters;
        bottomSubs{end+1} = subDiffArr{k1};
        bottomPsSub(end+1) = (randIters-(bottom-1))/randIters;
    end
end
writeData({topPsSub},['/mnt/vdb/home/ubuntu2/allSubDiffPs' suffix 'Top.txt'],'\t',{'pval'});
writeData({bottomPsSub},['/mnt/vdb/home/ubuntu2/allSubDiffPs' suffix 'Bottom.txt'],'\t',{'pval'});
for k=1:length(subDiffArr)
    strrepString = strrep(strrep(subDiffArr{k},'/','_'),' ','_');
    subDiffKthMapVals = subDiffKthMap(subDiffArr{k});
    writeData({subDiffKthMapVals},['/mnt/vdb/home/ubuntu2/mergedModelDists' suffix '/' strrepString 'KthMap.txt'],'\t',{'fluxdiff'});
    subDiffMapIndividualVals = subDiffMapIndividual(subDiffArr{k});
    highlightArr = {};
    for k1=1:length(subDiffMapIndividualVals)
	if any(k1==NormalIdxs)
	    highlightArr{k1} = 'normal';
        else
            %subDiffMapIndividualVals(k1) = subDivvMapIndividualVals(k1)*3;
            if useIBD
	        highlightArr{k1} = 'IBD';
            else
                highlightArr{k1} = 'obese';
            end
	end
    end
    writeData({subDiffMapIndividualVals,highlightArr},['/mnt/vdb/home/ubuntu2/mergedModelDists' suffix '/' strrepString 'DiffMap.txt'],'\t',{'flux','highlight'});
end
[~, sortIdxs] = sort(abs(subDiffArrNum),'descend');
subDiffArrNum = subDiffArrNum(sortIdxs);
subDiffArr = subDiffArr(sortIdxs);
highlightArr = {};
for k=1:length(subDiffArrNum)
    if topPsSub(strcmp(subDiffArr{k},topSubs)) < .05
	highlightArr{k} = 'yes';
    elseif bottomPsSub(strcmp(subDiffArr{k},bottomSubs)) < .05
        highlightArr{k} = 'yes';
    else
        highlightArr{k} = 'no';
    end
end
writeData({addIdxStrings(subDiffArr),subDiffArrNum,highlightArr},['/mnt/vdb/home/ubuntu2/mergedModel' suffix 'SubDiffs.txt'],'\t',{'sub','diffflux','highlight'});
[~, sortIdxs] = sort(abs(rxnDiffArrNum),'descend');
rxnDiffArrNum = rxnDiffArrNum(sortIdxs);
rxnDiffArr = rxnDiffArr(sortIdxs);
rxnNameDiffArr = {};
rxnTsArr = [];
rxnNumsMetricArr = [];
rxnPsMetricArr = [];
rxnPermPValsArr = [];
rxnExpDiffArr = [];
rxnExpPValsArr = [];
rxnExpNumsMetricArr = [];
rxnSubArr = {};
highlightArr = {};
for k=1:length(rxnDiffArrNum)
    rxnSubArr{k} = tobemerged2.subSystems{strcmp(tobemerged.rxns,rxnDiffArr{k})};
    rxnTsArr(k) = topPsRxnT(strcmp(topRxnsT,rxnDiffArr{k}));
    rxnNumsMetricArr(k) = topNumsRxnMetric(strcmp(topRxnsMetric,rxnDiffArr{k}));
    rxnPsMetricArr(k) = topPsRxnMetric(strcmp(topPsRxnMetricNames,rxnDiffArr{k}));
    rxnExpNumsMetricArr(k) = topNumsExpRxnMetric(strcmp(topExpRxnsMetric,rxnDiffArr{k}));
    rxnNameDiffArr{k} = tobemerged.rxnNames{strcmp(tobemerged.rxns,rxnDiffArr{k})};
    rxnPermPValsArr(k) = min(topPsRxn(strcmp(rxnDiffArr{k},topRxns)),bottomPsRxn(strcmp(rxnDiffArr{k},bottomRxns)));
    rxnExpDiffArr(k) = expRxnDiffArrNum(strcmp(rxnDiffArr{k},expRxnDiffArr));
    rxnExpPValsArr(k) = min(topPsExpRxn(strcmp(rxnDiffArr{k},topExpRxns)),bottomPsExpRxn(strcmp(rxnDiffArr{k},bottomExpRxns)));
    if topPsRxn(strcmp(rxnDiffArr{k},topRxns)) < .05
	highlightArr{k} = 'yes';
    elseif bottomPsRxn(strcmp(rxnDiffArr{k},bottomRxns)) < .05
        highlightArr{k} = 'yes';
    else
        highlightArr{k} = 'no';
    end
end
manovaarr = [];
manovastatusarr = {};
group = {};
for i=1:length(NormalIdxs)
    group{NormalIdxs(i)} = 'normal';
end
for i=1:length(IBDIdxs)
    group{IBDIdxs(i)} = 'IBD';
end
for i=1:length(rxnDiffArr)
    X = [expRxnDiffMapIndividual(rxnDiffArr{i}) rxnDiffMapIndividual(rxnDiffArr{i})];
    if sum(sum(abs(X)))==0
        manovastatusarr{i} = 'both zero';
        manovaarr(i)=NaN;
    elseif sum(sum(abs(X(:,1))))==0
        manovastatusarr{i} = 'expression zero';
        manovaarr(i)=NaN;
    elseif sum(sum(abs(X(:,2))))==0
        manovastatusarr{i} = 'flux zero';
        manovaarr(i)=NaN;
    else
        [d p] = manova1(X,group);
        manovastatusarr{i} = 'both nonzero';
        manovaarr(i) = p;
    end
end
[r p] = corr(rxnPermPValsArr',rxnExpPValsArr','type','Spearman');
r
p
[r p] = corr(rxnPermPValsArr(rxnPermPValsArr~=1 & rxnExpPValsArr~=1)',rxnExpPValsArr(rxnPermPValsArr~=1 & rxnExpPValsArr~=1)','type','Spearman');
r
p
[r p] = corr(rxnTsArr(~isnan(rxnTsArr))',rxnPermPValsArr(~isnan(rxnTsArr))','type','Spearman');
r
p
writeData({rxnPermPValsArr,rxnExpPValsArr},'/mnt/vdb/home/ubuntu2/fluxExpNormObeseCorr.txt','\t',{'fluxpval','exppval'});
[r p] = corr(rxnPermPValsArr(~isnan(rxnNumsMetricArr))',rxnNumsMetricArr(~isnan(rxnNumsMetricArr))','type','Spearman');
r
p
writeData({addIdxStrings(rxnDiffArr),manovaarr,manovastatusarr,rxnExpDiffArr},['/mnt/vdb/home/ubuntu2/manova' suffix '.txt'],'\t',{'rxn','manova p-val','manova status','expression diff'});
writeData({rxnPermPValsArr(~isnan(rxnNumsMetricArr))',rxnNumsMetricArr(~isnan(rxnNumsMetricArr))'},['/mnt/vdb/home/ubuntu2/NormObesePermMetricCorr.txt'],'\t',{'permpval','metricnum'});
writeData({rxnPermPValsArr(~isnan(rxnPsMetricArr))',rxnPsMetricArr(~isnan(rxnPsMetricArr))'},['/mnt/vdb/home/ubuntu2/NormObesePermMetricCorr2.txt'],'\t',{'permpval','metricpval'});
writeData({addIdxStrings(rxnDiffArr),rxnDiffArrNum,rxnSubArr,rxnNameDiffArr,rxnPermPValsArr,rxnTsArr,rxnExpPValsArr,rxnNumsMetricArr,rxnExpNumsMetricArr,highlightArr},['/mnt/vdb/home/ubuntu2/mergedModel' suffix 'RxnDiffs.txt'],'\t',{'rxn','diffflux','subsystem','rxnname','permutation p-val','t-test p-val','expression t-test p-val','flux metric','expression metric','highlight'});
meanPValSubDiffArr = [];
rxnNumSubArr = [];
rxnNumSigSubArr = [];
for k=1:length(subDiffArr)
    subRxns = tobemerged2.rxns(strcmp(tobemerged2.subSystems,subDiffArr{k}));
    M = length(subRxns);
    N = length(rxnDiffArr);
    K = sum(rxnPermPValsArr < .05);
    [~,~,intersectIdxs] = intersect(subRxns,rxnDiffArr);
    x = sum(rxnPermPValsArr(intersectIdxs) < .05);
    rxnNumSubArr(k) = M;
    rxnNumSigSubArr(k) = x;
    meanPValSubDiffArr(k) = mean(rxnPermPValsArr(intersectIdxs));
end
writeData({addIdxStrings(subDiffArr),subDiffArrNum,rxnNumSubArr,rxnNumSigSubArr,meanPValSubDiffArr},['/mnt/vdb/home/ubuntu2/mergedModel' suffix 'SubDiffs2.txt'],'\t',{'sub','diffflux','total number of rxns','number of significant rxns','mean p-val of reactions in subsystem'});

for k=1:length(rxnDiffArr)
    if rxnPermPValsArr(k) < .05
    strrepString = strrep(strrep(rxnDiffArr{k},'/','_'),' ','_');
    rxnDiffKthMapVals = rxnDiffKthMap(rxnDiffArr{k});
    writeData({rxnDiffKthMapVals},['/mnt/vdb/home/ubuntu2/mergedModelRxnDists' suffix '/' strrepString 'KthMap.txt'],'\t',{'fluxdiff'});
    rxnDiffMapIndividualVals = rxnDiffMapIndividual(rxnDiffArr{k});
    highlightArr = {};
    for k1=1:length(rxnDiffMapIndividualVals)
	if any(k1==NormalIdxs)
	    highlightArr{k1} = 'normal';
        else
            %rxnDiffMapIndividualVals(k1) = rxnDivvMapIndividualVals(k1)*3;
            if useIBD
	        highlightArr{k1} = 'IBD';
            else
                highlightArr{k1} = 'obese';
            end
	end
    end
    writeData({rxnDiffMapIndividualVals,highlightArr},['/mnt/vdb/home/ubuntu2/mergedModelRxnDists' suffix '/' strrepString 'DiffMap.txt'],'\t',{'flux','highlight'});
    end
end
end

if 0
disp('HERE')
subDiffMapNormal = containers.Map;
subDiffMapIBD = containers.Map;
for z=1:5%length(sampleArr)
    for i=1:size(IBDDataMat,1)
	for j=1:size(IBDDataMat,2)
	    AGORAModel = AGORAModelArr{j}
	    for k=1:length(AGORAModel.subSystems)
		AGORAModel.subSystems{k} = AGORAModel.subSystems{k}{1};
            end
	    [subLabelsArr,subFluxSums] = segmentFluxBySubsystem(AGORAModel,IBDDataMat{i,j,z});
            for k=1:length(subLabelsArr)
	        if sum(strcmp(z,IBDIdxs))~=0
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
    end
end
subKeys = keys(subDiffMapNormal);
subDiffArr = {};
subDiffArrNum = [];
for i=1:length(subKeys)
	subDiffArr{i} = subKeys{i};
subDiffArrNum(i) = subDiffMapNormal(subKeys{i}) - subDiffMapIBD(subKeys{i});
end
[~, sortIdxs] = sort(abs(subDiffArrNum),'descend');
subDiffArrNum = subDiffArrNum(sortIdxs);
subDiffArr = subDiffArr(sortIdxs);
writeData({subDiffArr,subDiffArrNum},['/mnt/vdb/home/ubuntu2/IBDSubDiffs.txt'],'\t',{'sub','diffflux'});
end

IBDCentsMat = [];
NormalCentsMat = [];
for z=1:size(allInteractMat,1)
    dependencyZth = abs(allInteractMat(z,:,:,2)-allInteractMat(z,:,:,1));
    dependencyZth = squeeze(dependencyZth);
    dependencyZth = dependencyZth(1:min(size(dependencyZth)),1:min(size(dependencyZth)));
    dependencyCentsZth = betweenness_centrality(sparse(squeeze(dependencyZth)*1.0));
    if sum(IBDIdxs==z)~=0
        IBDCentsMat(end+1,:) = dependencyCentsZth;
    else
        NormalCentsMat(end+1,:) = dependencyCentsZth;
    end
end
[NormalCentsArr, sortIdxs] = sort(mean(NormalCentsMat,1),'descend');
IBDCentsArr = mean(IBDCentsMat,1); IBDCentsArr = IBDCentsArr(sortIdxs);

for i=length(NormalCentsArr):773
    NormalCentsArr(i) = 0;
    IBDCentsArr(i) = 0;
end

NormalCentsArr = NormalCentsArr./sum(NormalCentsArr);
IBDCentsArr = IBDCentsArr./sum(IBDCentsArr);

doubledXCoords = ones(length(IBDCentsArr)*2,1);
doubledXCoords(1:2:length(doubledXCoords)) = 1:length(IBDCentsArr);
doubledXCoords(2:2:length(doubledXCoords)) = 1:length(IBDCentsArr);
groupArr = {};
for i=1:length(doubledXCoords)
    if mod(i,2)==1
	groupArr{end+1} = 'Normal';
    else
        groupArr{end+1} = 'IBD';
    end
end
writeData({doubledXCoords,[NormalCentsArr IBDCentsArr],groupArr},['/mnt/vdb/home/ubuntu2/IBDCents.txt'],'\t',{'modelid','cent','group'});

eigs = [];
for z=1:5%length(sampleArr)
    eigMatrix = squeeze(allInteractMat(z,:,:,2)-allInteractMat(z,:,:,1));
    eigMatrix = eigMatrix(1:min(size(eigMatrix)),1:min(size(eigMatrix)));
    eigZth = eig(eigMatrix);
    eigZth = real(eigZth); [~, sortIdxs] = sort(abs(eigZth),'descend'); eigZth = eigZth(sortIdxs);
    eigs(:,z) = eigZth*-1;
end
eiglabels = {'eig01','eig02','eig03','eig04','eig05','eig06','eig07','eig08','eig09','eig10'};
writeData({eiglabels,mean(eigs(1:10,IBDIdxs),2),max(eigs(1:10,IBDIdxs),[],2),min(eigs(1:10,IBDIdxs),[],2)},'/mnt/vdb/home/ubuntu2/AGORAIBDEigsIBD.txt','\t',{'eigidx','meaneig','upper','lower'});
writeData({eiglabels,mean(eigs(1:10,NormalIdxs),2),max(eigs(1:10,NormalIdxs),[],2),min(eigs(1:10,NormalIdxs),[],2)},'/mnt/vdb/home/ubuntu2/AGORAIBDEigsNormal.txt','\t',{'eigidx','meaneig','upper','lower'});
