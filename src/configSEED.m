baseDir = '/home/ubuntu/MATLAB/SEED';
inputDir = '/home/ubuntu/MATLAB/SEED/input';
outputDir = '/home/ubuntu/MATLAB/SEED/output';
modelsDir = '/home/ubuntu/MATLAB/SEED/input/SEEDModels';
GreenblumDir = '/home/ubuntu/MATLAB/SEED/input/GreenblumData';
FI = fopen([GreenblumDir filesep 'GreenblumECs.txt']);
GreenblumEC = textscan(FI,'%s\n');
fclose(FI);
GreenblumEC = GreenblumEC{1};
if ~exist('CobraLPSolver','var')
    %initCobraToolbox;
end
if exist([outputDir filesep 'makeBigModelAccum' filesep 'makeBigModelAccum.mat'],'file')
    load([outputDir filesep 'makeBigModelAccum' filesep 'makeBigModelAccum.mat'],'bigModelAccum','bigModelsAccum','modelNamesToModels','rxnsToECsAccum','ECsToRxnsAccum');
end
if exist([outputDir filesep 'makeBigModelTable' filesep 'makeBigModelTable.mat'],'file')
    load([outputDir filesep 'makeBigModelTable' filesep 'makeBigModelTable.mat'],'bigModelTable','rxnsToECsTable','ECsToRxnsTable');
end
if exist([outputDir filesep 'readMD5Express' filesep 'readMD5Express.mat'],'file')
    load([outputDir filesep 'readMD5Express' filesep 'readMD5Express.mat'],'rxnsToExpressObese','rxnsToExpressNorm','rxnsToExpressNormKeys','rxnsToDiffExp');
end
