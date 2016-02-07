outputDir1 = [inputDir filesep 'SEEDModels'];
FI = fopen([inputDir filesep 'allModels.txt']);
dataFields = textscan(FI,'%s','Delimiter','\n');
dataFields = dataFields{1};
diagnostic = 1;
for i=1:length(dataFields)
    if diagnostic
        dispString = dataFields{i};
        [status result] = system(sprintf('cat %s | wc -l',[outputDir1 filesep dataFields{i} '.tsv']));
        if ~isempty(regexp(result,'No such file')) || str2num(result)<=10
            dispString = [dispString ' xls missing'];
            %disp(dispString)
        end
        [status result] = system(sprintf('tail -n 1 %s',[outputDir1 filesep dataFields{i} '.xml']));
        %disp(result)
        if ~isempty(regexp(result,'No such file')) || ~isempty(regexp(result,'html'))
            dispString = [dispString ' xml missing'];
        end
        disp(dispString);
    else        
        status=system(sprintf('curl --data "model=%s&file=XLS" seed-viewer.theseed.org/ModelSEEDdownload.cgi > %s',dataFields{i},[outputDir1 filesep dataFields{i} '.xls']));
        status=system(sprintf('curl --data "model=%s&file=SBML" seed-viewer.theseed.org/ModelSEEDdownload.cgi > %s',dataFields{i},[outputDir1 filesep dataFields{i} '.xml']));
    end
end
