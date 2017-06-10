baseDir = '/home/fs01/yw595/MATLAB/SEED';
inputDir = [baseDir filesep 'input'];
outputDir = [baseDir filesep 'output'];
transferDir = [baseDir filesep 'transfer'];
modelsDir = [inputDir filesep 'SEEDModels'];
GreenblumDir = [inputDir filesep 'GreenblumData'];
if ~exist('GreenblumEC','var')
    FI = fopen([GreenblumDir filesep 'GreenblumECs.txt']);
    GreenblumEC = textscan(FI,'%s\n');
    fclose(FI);
    GreenblumEC = GreenblumEC{1};
end
if ~exist('cpdData','var')
    FI = fopen([inputDir filesep 'compoundsCleaned.csv']);
    dataFields = textscan(FI,repmat('%s',1,10),'Delimiter',',');
    cpdData = [dataFields{:}];
    cpdIDs = cpdData(:,1); cpdNames = cpdData(:,2); cpdAbbrvs = cpdData(:,3); cpdKEGGs = cpdData(:,5);
    fclose(FI);
    FI = fopen([inputDir filesep 'reactionsCleaned.csv']);
    dataFields = textscan(FI,repmat('%s',1,9),'Delimiter',',');
    rxnData = [dataFields{:}]; fclose(FI);
    rxnIDs = rxnData(:,1); rxnNames = rxnData(:,2); equations = rxnData(:,7);
    equations = cellfun(@(x) strrep(strrep(x,'(',''),')',''), equations, 'UniformOutput',0);
    fclose(FI);
end
if ~exist('metabolomeData','var')
    ERPData = textscan(fopen('/home/fs01/yw595/MATLAB/SEED/input/MGMData/ERPBowtieComplete.txt'),'%s%s%s%s','Delimiter','\t','HeaderLines',0);
    xenoData = textscan(fopen('/home/fs01/yw595/MATLAB/SEED/output/writeXenoData/GSM935962_A1_EtOH_CDS.EC.txt'),'%s%s%s','Delimiter','\t','HeaderLines',0);
    metabolomeData = textscan(fopen('/home/fs01/yw595/MATLAB/SEED/input/fecalmicrobiomeauto.txt'),'%s%s','Delimiter','\t','HeaderLines',0);
    hadzaData = textscan(fopen('/home/fs01/yw595/output3.EC.txt'),'%s%s','Delimiter','\t','HeaderLines',0);
end
if ~exist('closestFamiliesData','var')
    closestFamiliesData = textscan(fopen('/home/fs01/yw595/closestFamilies.txt'),'%s%s','Delimiter','|','HeaderLines',0);
    abundsData = textscan(fopen('/home/fs01/yw595/MATLAB/SEED/src/HMPFamiliesAbundsAndOccurs.txt'),'%s%s%s','Delimiter','\t','HeaderLines',0);
    allAbundsData = textscan(fopen('/home/fs01/yw595/MATLAB/SEED/src/HMPAllAbunds.txt'),repmat('%s',1,3840),'Delimiter','\t','HeaderLines',0);
    closestFamiliesRevData = textscan(fopen('/home/fs01/yw595/closestFamiliesRev.txt'),'%s%s','Delimiter','|','HeaderLines',0);
    allAbundsDataMatrix = cell2mat(cellfun(@(x) cellfun(@(y) str2num(y), x), allAbundsData(2:end),'UniformOutput',0));
    allAbundsDataNames = allAbundsData{1};
end
if ~exist('initedCobra','var') || initedCobra==0
    initCobraToolbox;
    initedCobra=1;
end
if exist([outputDir filesep 'readTIGRFams' filesep 'readTIGRFams.mat'],'file') && ~exist('TIGRIDsToECNums','var')
    load([outputDir filesep 'readTIGRFams' filesep 'readTIGRFams.mat'],'TIGRIDsToECNums','TIGRIDsToDEs','DEsToTIGRIDs','flatDEs');
end
if exist([outputDir filesep 'readGOMaps' filesep 'readGOMaps.mat'],'file') && ~exist('pfamToEC','var')
    load([outputDir filesep 'readGOMaps' filesep 'readGOMaps.mat'],'pfamToEC','pfamToPfamLong','pfamLongToPfam','flatPfamLong');
end
if exist([outputDir filesep 'mapHMPEC' filesep 'mapHMPEC.mat'],'file') && ~exist('allECs','var')
    load([outputDir filesep 'mapHMPEC' filesep 'mapHMPEC.mat'],'allECs','percentageGeneCovs','percentageECCovs','speciesNames','allECLengths','flatTableECNums','geneNameToEC','HMPIDsToGeneNames');
end
if exist([outputDir filesep 'writeXenoECs' filesep 'writeXenoECs.mat'],'file') && ~exist('xenoExpECs','var')
    load([outputDir filesep 'writeXenoECs' filesep 'writeXenoECs.mat'],'xenoExpECs','xenoExps','xenoExpStds');
end
if exist([outputDir filesep 'makeBigModelAccum' filesep 'makeBigModelAccum.mat'],'file') && ~exist('bigModelAccum','var')
    load([outputDir filesep 'makeBigModelAccum' filesep 'makeBigModelAccum.mat'],'bigModelAccum','bigModelsAccumSubsystems','modelNamesToModels','modelNames','modelNamesShort','rxnsToECsAccum','ECsToRxnsAccum');
end
if exist([outputDir filesep 'makeBigModelTable' filesep 'makeBigModelTable.mat'],'file') && ~exist('bigModelTable','var')
    load([outputDir filesep 'makeBigModelTable' filesep 'makeBigModelTable.mat'],'bigModelTable','rxnsToECsTable','ECsToRxnsTable','bigModel','bigModelAdded','testModel');
end
if exist([outputDir filesep 'readMD5Express' filesep 'readMD5Express.mat'],'file') && ~exist('rxnsToExpressObese','var')
    load([outputDir filesep 'readMD5Express' filesep 'readMD5Express.mat'],'rxnsToExpressObese','rxnsToExpressNorm','rxnsToExpressNormKeys','rxnsToDiffExp');
end
if exist([outputDir filesep 'reconcBigModelRecon2' filesep 'reconcBigModelRecon2.mat'],'file') && ~exist('bigModelReconc','var')
    load([outputDir filesep 'reconcBigModelRecon2' filesep 'reconcBigModelRecon2.mat'],'bigModelReconc');
end
if exist([outputDir filesep 'examineBiomass' filesep 'examineBiomass.mat'],'file') && ~exist('allBiomassRates','var')
    load([outputDir filesep 'examineBiomass' filesep 'examineBiomass.mat'],'allBiomassRates','allShadMets','allBiomassDists');
end
if exist([outputDir filesep 'mergeSmallModels' filesep 'mergeSmallModels.mat'],'file') && ~exist('allTwoBiomasses','var')
    load([outputDir filesep 'mergeSmallModels' filesep 'mergeSmallModels.mat'],'allTwoBiomasses');
    allTwoBiomassesTemp = [];
    while length(allTwoBiomasses)<length(modelNames)*length(modelNames)
        allTwoBiomasses(end+1) = 0;%rand(1)*1000;
    end
    zCount = 1;
    for i=1:length(modelNames)
        for j=1:length(modelNames)
            allTwoBiomassesTemp(i,j) = allTwoBiomasses(zCount);
            zCount = zCount+1;
        end
    end
    allTwoBiomasses = allTwoBiomassesTemp;
end
if exist([outputDir filesep 'readHMPTaxaData' filesep 'readHMPTaxaData.mat'],'file') && ~exist('closestFamiliesData','var')
    load([outputDir filesep 'readHMPTaxaData' filesep 'readHMPTaxaData.mat']);
end






