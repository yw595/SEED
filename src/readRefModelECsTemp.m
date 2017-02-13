configSEED;
inputFI = fopen([inputDir filesep 'ShoaieRefInfo' filesep 'iAF1260.txt']);
dataFields = textscan(inputFI,'%s');
dataFields = [dataFields{:}];
iAF1260ECs = unique(dataFields);
fclose(inputFI);

inputFI = fopen([inputDir filesep 'ShoaieRefInfo' filesep 'iMH551.txt']);
dataFields = textscan(inputFI,'%s');
dataFields = [dataFields{:}];
iMH551ECs = unique(dataFields);
iMH551ECs = cellfun(@(x) strsplitYiping(x,', '), iMH551ECs, 'UniformOutput',0); iMH551ECs = [iMH551ECs{:}]; iMH551ECs = unique(iMH551ECs);
iMH551ECs = cellfun(@(x) strsplitYiping(x,'/ '), iMH551ECs,'UniformOutput',0); iMH551ECs = [iMH551ECs{:}]; iMH551ECs = unique(iMH551ECs);
iMH551ECs = cellfun(@(x) strsplitYiping(x,'/'), iMH551ECs,'UniformOutput',0); iMH551ECs = [iMH551ECs{:}]; iMH551ECs = unique(iMH551ECs);
iMH551ECs = cellfun(@(x) strsplitYiping(x,' and '), iMH551ECs,'UniformOutput',0); iMH551ECs = [iMH551ECs{:}]; iMH551ECs = unique(iMH551ECs);
iMH551ECs = cellfun(@(x) strsplitYiping(x,' or '), iMH551ECs,'UniformOutput',0); iMH551ECs = [iMH551ECs{:}]; iMH551ECs = unique(iMH551ECs);
iMH551ECs = cellfun(@(x) strsplitYiping(x,','), iMH551ECs,'UniformOutput',0); iMH551ECs = [iMH551ECs{:}]; iMH551ECs = unique(iMH551ECs);
iMH551ECs = cellfun(@(x) strsplitYiping(x,';'), iMH551ECs,'UniformOutput',0); iMH551ECs = [iMH551ECs{:}]; iMH551ECs = unique(iMH551ECs);
iMH551ECs = iMH551ECs(cellfun(@(x) isempty(regexp(x,'A')), iMH551ECs));

inputFI = fopen([inputDir filesep 'ShoaieRefInfo' filesep 'iBif452.txt']);
dataFields = textscan(inputFI,'%s');
dataFields = [dataFields{:}];
iBif452ECs = unique(dataFields);
iBif452ECs = cellfun(@(x) strsplitYiping(x,';'), iBif452ECs, 'UniformOutput',0); iBif452ECs = [iBif452ECs{:}]; iBif452ECs = unique(iBif452ECs);
iBif452ECs = cellfun(@(x) strsplitYiping(x,':'), iBif452ECs, 'UniformOutput',0); iBif452ECs = [iBif452ECs{:}]; iBif452ECs = unique(iBif452ECs);
iBif452ECs = iBif452ECs(cellfun(@(x) isempty(regexp(x,'K')), iBif452ECs));
iBif452ECs = iBif452ECs(cellfun(@(x) isempty(regexp(x,'A')),iBif452ECs));

modelRAVENKEGG = load([inputDir filesep 'RAVEN' filesep 'kegg' filesep 'keggRxns.mat'],'model'); modelRAVENKEGG = modelRAVENKEGG.model;
RAVENKEGGECs = modelRAVENKEGG.eccodes;
RAVENKEGGECs = cellfun(@(x) strsplitYiping(x,' '), RAVENKEGGECs,'UniformOutput',0); RAVENKEGGECs = [RAVENKEGGECs{:}]; RAVENKEGGECs = unique(RAVENKEGGECs);

flatTableECNums = unique([bigModelTable.rxnECNums{:}]);

iAF1260intersect = length(intersect(iAF1260ECs,flatTableECNums));
iMH551intersect = length(intersect(iMH551ECs,flatTableECNums));
iBif452intersect = length(intersect(iBif452ECs,flatTableECNums));
RAVENKEGGintersect = length(intersect(RAVENKEGGECs,flatTableECNums));

xvals = 1:3; yvals = [iAF1260intersect length(iAF1260ECs)-iAF1260intersect;iMH551intersect length(iMH551ECs)-iMH551intersect;iBif452intersect length(iBif452ECs)-iBif452intersect];%;RAVENKEGGintersect length(RAVENKEGGECs)-RAVENKEGGintersect];
titleString = 'Coverage of Previous Models';
xlabels = {'iAF1260', 'iMH551','iBif452'}; grouplabels = {'shared with large model','not shared with large model'};
makeBar(xvals,yvals,titleString,outputDir,'isStackBar',1, 'ylabelString','Num ECs covered','xlabels',xlabels,'legendLabels',grouplabels,'titleFontSize',5,'ylabelFontSize',5);
xvalsTemp = []; yvalsTemp = []; xlabelsTemp = {}; grouplabelsTemp = {};
for i=1:size(yvals,1)
    for j=1:size(yvals,2)
        xvalsTemp(end+1) = xvals(i);
        yvalsTemp(end+1) = yvals(i,j);
        xlabelsTemp{end+1} = xlabels{i};
        grouplabelsTemp{end+1} = grouplabels{j};
    end
end
xvals = xvalsTemp; yvals = yvalsTemp; xlabels = xlabelsTemp; grouplabels = grouplabelsTemp;
dataFields = {};
dataFields{1} = arrayfun(@(x) num2str(x), xvals, 'UniformOutput', 0);
dataFields{2} = arrayfun(@(x) num2str(x), yvals, 'UniformOutput', 0);
dataFields{3} = xlabels;
dataFields{4} = grouplabels
colHeaders = {'xvals','yvals','xlabels','grouplabels'};
for i=1:length(dataFields)
    ithData = dataFields{i};
    ithData(2:end+1) = ithData(1:end);
    ithData{1} = colHeaders{i};
    dataFields{i} = ithData;
end
writeData(dataFields,[outputDir filesep 'Coverage_of_Previous_Models'],'\t');

xvals = 1:length(speciesNames);
yvals = percentageECCovs(:,2);
titleString = 'EC Coverage Among HMP Species';
makeBar(xvals,yvals,titleString,outputDir,'ylabelString','% ECs covered','xlabels',speciesNames);

xvals = 1:length(speciesNames);
yvals = percentageGeneCovs(:,2);
titleString = 'Gene Coverage Among HMP Species';
makeBar(xvals,yvals,titleString,outputDir,'ylabelString','% Genes covered','xlabels',speciesNames);










