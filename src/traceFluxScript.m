beginRxns = {'VALTA','LEUTA','ILETA','HISTD'};
returnsizes = [];
for j=1:length(beginRxns)
    beginRxn = beginRxns{j};
    excluderxns = {};
    outFI = fopen(['/mnt/vdb/home/ubuntu2/' beginRxn 'fluxes.txt'],'w');
    for i=1:10
	excluderxns{end+1} = beginRxn;
	excluderxns = unique(excluderxns);
	beginRxn
	[returnrxns, returnfluxes] = traceFlux(tobemerged3,fluxdist,beginRxn,excluderxns);
returnsizes(end+1) = length(returnrxns);
        if isempty(returnrxns)
            break;
        end
	fprintf(outFI,'%s\t',beginRxn);
	fprintf(outFI,'%s\n',printRxnEq(tobemerged3,beginRxn));
        fprintf(outFI,'%s\n',tobemerged3.subSystems{strcmp(tobemerged3.rxns,beginRxn)});
	fprintf(outFI,'%f\n',fluxdist(strcmp(tobemerged3.rxns,beginRxn)));
	beginRxn = returnrxns{1};
    end
    fclose(outFI);
end
