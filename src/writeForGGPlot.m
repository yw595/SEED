function writeForGGPlot(xvals,yvals,outputFile,xlabels,grouplabels)

if length(size(yvals))==2 && ~any(size(yvals)==1)
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
end
dataFields = {};
dataFields{1} = arrayfun(@(x) num2str(x), xvals, 'UniformOutput', 0);
dataFields{2} = arrayfun(@(x) num2str(x), yvals, 'UniformOutput', 0);
if exist('xlabels','var')
    dataFields{3} = xlabels;
end
if exist('grouplabels','var')
    dataFields{4} = grouplabels;
end
colHeaders = {'xvals','yvals','xlabels','grouplabels'};
for i=1:length(dataFields)
    ithData = dataFields{i};
    ithData(2:end+1) = ithData(1:end);
    ithData{1} = colHeaders{i};
    dataFields{i} = ithData;
end
writeData(dataFields,outputFile,'\t');

end