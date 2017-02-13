CytoFI = fopen([baseDir filesep 'SEEDMine.txt'],'w');
for i=1:length(ECListTableNoTE)
    for j=1:length(ECListTableNoTE)
        if ECMatrixTableNoTE(i,j)~=0 && i~=j
            fprintf(CytoFI,sprintf('%s pp %s\n',ECListTableNoTE{i},ECListTableNoTE{j}));
        end
    end
end
fclose(CytoFI);