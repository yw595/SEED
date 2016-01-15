inputDir = '/home/ubuntu/sas';
outputDir = '/mnt/extra/SEED';
FI = fopen([inputDir filesep 'allModels.txt']);
dataFields = textscan(FI,'%s','Delimiter','\n');
dataFields = dataFields{1};
for i=1:length(dataFields)
    status=system(sprintf('curl --data "model=%s&file=XLS" seed-viewer.theseed.org/ModelSEEDdownload.cgi > %s',dataFields{i},[outputDir filesep dataFields{i} '.xls']));
    status=system(sprintf('curl --data "model=%s&file=SBML" seed-viewer.theseed.org/ModelSEEDdownload.cgi > %s',dataFields{i},[outputDir filesep dataFields{i} '.xml']));
end
