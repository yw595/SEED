configSEED;
outputDir1 = [outputDir filesep 'readHMPTaxaData'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1);
end

closestAbunds = [];
closestOccurs = [];

for i=1:length(modelNamesShort)
    speciesName = modelNamesShort{i};
    if sum(strcmp(closestFamiliesData{1},speciesName))~=0
        closestFamily = closestFamiliesData{2}{strcmp(closestFamiliesData{1},speciesName)};
        disp(closestFamily);
        closestAbund = abundsData{2}{strcmp(abundsData{1},closestFamily)};
        closestOccur = abundsData{3}{strcmp(abundsData{1},closestFamily)};
        closestAbunds(i) = str2num(closestAbund);
        closestOccurs(i) = str2num(closestOccur);
        disp(closestAbund);
        disp(closestOccur);
    end
end

save([outputDir1 filesep 'readHMPTaxaData.mat'],'closestFamiliesData','abundsData','allAbundsData','closestFamiliesRevData','closestAbunds','closestOccurs');