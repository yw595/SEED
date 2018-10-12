function [compTerm coopTerm coopBasicTerm compBasicTerm] = metInteractFuncTemp2(modelIth,modelJth,fluxesIth,fluxesJth)
            compTerm = 0;
	    coopTerm = 0;
            coopBasicTerm = 0;
            compBasicTerm = 0;
	    for k=1:length(modelIth.rxns)
		if ~isempty(regexp(modelIth.subSystems{k},'Exchange'))
		    if ~isempty(strcmp(modelJth.rxns,modelIth.rxns{k}))
			correspondK = find(strcmp(modelJth.rxns,modelIth.rxns{k}));
			if fluxesIth(k) < 0 && mean(fluxesJth(correspondK)) < 0
			    compTerm = compTerm + fluxesIth(k)*mean(fluxesJth(correspondK));
                            compBasicTerm = compBasicTerm + 1;
			end
			if (fluxesIth(k) < 0 && mean(fluxesJth(correspondK)) > 0) || (fluxesIth(k) > 0 && mean(fluxesJth(correspondK)) < 0)
			    coopTerm = coopTerm - fluxesIth(k)*mean(fluxesJth(correspondK));
                            coopBasicTerm = coopBasicTerm + 1;
                        end
		    end
		end
	    end
