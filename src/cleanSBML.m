outputDir = '/mnt/extra/SEED';
filenames = dir(outputDir);
for i=1:length(filenames)
    if ~isempty(regexp(filenames(i).name,'^i.*.xml$'))
        disp(filenames(i).name)
        modelName = filenames(i).name(1:regexp(filenames(i).name,'.xml')-1);
        disp(modelName)
        trCmd = sprintf('remove250_1.sh %s %s %s',[outputDir filesep filenames(i).name],[outputDir filesep modelName '_2.xml']);
        status=system(trCmd);
        trCmd = sprintf('remove250_2.sh %s %s %s',[outputDir filesep modelName '_2.xml'],[outputDir filesep modelName '_3.xml']);
        status=system(trCmd);
        rawMatrix = readVaryFile([outputDir filesep modelName '_3.xml']);
        rawMatrix = regexprep(rawMatrix,'reaction id="rxn\d+', '$0_2');
        writeData({rawMatrix},[outputDir filesep modelName '_4.xml']);
    end
end