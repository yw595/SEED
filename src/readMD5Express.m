configSEED;
outputDir1 = [outputDir filesep 'readMD5Express'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1);
end
rxnsToExpressObese = mapExpToRxns(ECsToRxns,[baseDir filesep 'testObese.txt']);
rxnsToExpressNorm = mapExpToRxns(ECsToRxns,[baseDir filesep 'testNorm.txt']);
rxnsToDiffExp = containers.Map;
rxnsToExpressNormKeys = keys(rxnsToExpressNorm);
for i=1:length(rxnsToExpressNormKeys)
    expKey = rxnsToExpressNormKeys{i};
    if isKey(rxnsToExpressObese,expKey)
        rxnsToDiffExp(expKey) = abs(log2( (rxnsToExpressObese(expKey)/37)/(rxnsToExpressNorm(expKey)/87) ));
    end
    i
end
save([outputDir1 filesep 'readMD5Express.mat'],'rxnsToExpressObese','rxnsToExpressNorm','rxnsToExpressNormKeys','rxnsToDiffExp');