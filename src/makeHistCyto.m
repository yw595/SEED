function makeHistCyto(model,fluxdist,outputfile,subsystem)
    rxnarr = {};
    metarr = {};
    fluxarr = [];
    for i=1:length(model.rxns)
	if strcmp(model.subSystems{i},subsystem)
	    for j=1:length(model.mets)
		if model.S(j,i)~=0
		    rxnarr{end+1} = [model.rxns{i} model.rxnNames{i}];
                    metarr{end+1} = [model.mets{j} model.metNames{j}];
                    fluxarr(end+1) = fluxdist(i)*model.S(j,i);
                end
            end
        end
    end
    writeData({rxnarr,metarr,fluxarr},outputfile,'\t',{'rxn','met','flux'});
