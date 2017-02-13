baseDir = '/home/fs01/yw595/MATLAB/SEED';
inputDir = [baseDir filesep 'input'];
outputDir = [baseDir filesep 'output'];
modelsDir = [inputDir filesep 'SEEDModels'];
GreenblumDir = [inputDir filesep 'GreenblumData'];
if ~exist('GreenblumEC','var')
    FI = fopen([GreenblumDir filesep 'GreenblumECs.txt']);
    GreenblumEC = textscan(FI,'%s\n');
    fclose(FI);
    GreenblumEC = GreenblumEC{1};
end
if ~exist('CobraLPSolver','var')
    %initCobraToolbox;
end
if exist([outputDir filesep 'makeBigModelAccum' filesep 'makeBigModelAccum.mat'],'file') && ~exist('bigModelAccum','var')
    load([outputDir filesep 'makeBigModelAccum' filesep 'makeBigModelAccum.mat'],'bigModelAccum','modelNamesToModels','rxnsToECsAccum','ECsToRxnsAccum');
end
if exist([outputDir filesep 'makeBigModelTable' filesep 'makeBigModelTable.mat'],'file') && ~exist('bigModelTable','var')
    load([outputDir filesep 'makeBigModelTable' filesep 'makeBigModelTable.mat'],'bigModelTable','rxnsToECsTable','ECsToRxnsTable','bigModel','bigModelAdded','testModel');
end
load([outputDir filesep 'makeBigModelAccum' filesep 'extra.mat']);
if exist([outputDir filesep 'readMD5Express' filesep 'readMD5Express.mat'],'file') && ~exist('rxnsToExpressObese','var')
    load([outputDir filesep 'readMD5Express' filesep 'readMD5Express.mat'],'rxnsToExpressObese','rxnsToExpressNorm','rxnsToExpressNormKeys','rxnsToDiffExp');
end