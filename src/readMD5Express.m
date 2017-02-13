configSEED;
outputDir1 = [outputDir filesep 'readMD5Express'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1);
end
%rxnsToExpressObese = mapExpToRxns(ECsToRxnsTable,[baseDir filesep 'testObese.txt']);
rxnsToExpressObese = mapExpToRxns(ECsToRxnsTable,[baseDir filesep 'expressObese.txt']);
%rxnsToExpressNorm = mapExpToRxns(ECsToRxnsTable,[baseDir filesep 'testNorm.txt']);
rxnsToExpressNorm = mapExpToRxns(ECsToRxnsTable,[baseDir filesep 'expressNorm.txt']);
rxnsToDiffExp = containers.Map;
rxnsToExpressNormKeys = keys(rxnsToExpressNorm);
for i=1:length(rxnsToExpressNormKeys)
    expKey = rxnsToExpressNormKeys{i};
    if isKey(rxnsToExpressObese,expKey)
        %rxnsToDiffExp(expKey) = abs(log2( (rxnsToExpressObese(expKey)/37)/(rxnsToExpressNorm(expKey)/87) ));
        rxnsToDiffExp(expKey) = abs(log2( (rxnsToExpressObese(expKey))/(rxnsToExpressNorm(expKey)) ));
    end
    i
end
save([outputDir1 filesep 'readMD5Express.mat'],'rxnsToExpressObese','rxnsToExpressNorm','rxnsToExpressNormKeys','rxnsToDiffExp');