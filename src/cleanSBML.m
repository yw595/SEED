configSEED;
outputDir1 = modelsDir;
filenames = dir(outputDir1);
for i=1:length(filenames)
    if ~isempty(regexp(filenames(i).name,'^i.*.xml$'))
        disp(filenames(i).name)
        modelName = filenames(i).name(1:regexp(filenames(i).name,'.xml')-1);
        disp(modelName)
        trCmd = sprintf('remove250_1.sh %s %s',[outputDir1 filesep filenames(i).name],[outputDir1 filesep modelName '_2.xml']);
        status=system(trCmd);
        trCmd = sprintf('remove250_2.sh %s %s',[outputDir1 filesep modelName '_2.xml'],[outputDir1 filesep modelName '_3.xml']);
        status=system(trCmd);
        rawMatrix = readVaryFile([outputDir1 filesep modelName '_3.xml']);
        rawMatrix = regexprep(rawMatrix,'reaction id="rxn\d+', '$0_2');
        writeData({rawMatrix},[outputDir1 filesep modelName '_4.xml']);
    end
end