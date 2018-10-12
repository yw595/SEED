if 0

AGORAFiles = dir([inputDir filesep 'AGORAModels/Western-Diet-Paper']);
AGORAMat = {};
count = 0;
for i = 1:length(AGORAFiles)
    dotIdx = regexp(AGORAFiles(i).name,'\.');
    if ~strcmp(AGORAFiles(i).name,'.') && ~strcmp(AGORAFiles(i).name,'..') && strcmp(AGORAFiles(i).name(dotIdx+1:end),'mat')
        AGORAMat{end+1} = [inputDir filesep 'AGORAModels/Western-Diet-Paper' filesep AGORAFiles(i).name];
    end
end

outputDir1 = [outputDir filesep 'simulateSmallModelsSeparatePicrustAGORA'];

ucrFolderSort = 'ucrC97sortmerna';
presentSpeciesSort = {'Bacteroides_massiliensis','Bacteroides_clarus','Lactobacillus_paracasei','Parabacteroides_merdae','Alistipes_shahii','Faecalibacterium_prausnitzii','Clostridium_hathewayi','Parvimonas_micra','Prevotella_salivae','Dysgonomonas_gadei','Paraprevotella_clara','Prevotella_bryantii'};
ucrFolderMH0007 = 'ucrC97MH0007'
presentSpeciesMH0007 = {'Lactobacillus_fermentum','Bacteroides_clarus','Collinsella_tanakaei','Acidaminococcus_intestini','Anaerotruncus_colihominis','Faecalibacterium_prausnitzii','Trueperella_pyogenes','Pseudomonas_aeruginosa','Clostridium_hathewayi','Parvimonas_micra','Staphylococcus_epidermidis','Hafnia_alvei','Paraprevotella_clara','Prevotella_bryantii','Alistipes_shahii','Fusobacterium_nucleatum','Slackia_exigua'};
modelsSort = {};
fluxesSort = {};
modelsMH0007 = {};
fluxesMH0007 = {};
for i=1:length(AGORAMat)
    i
    for j=1:length(AGORAMat)
	if i~=j && ~isempty(AGORAMat{i}) && ~isempty(AGORAMat{j})
	    isPresent1 = 0;
            isPresent2 = 0;
	    for k=1:length(presentSpeciesSort)
		if ~isempty(regexp(AGORAMat{i},presentSpeciesSort{k}))
		    isPresent1 = 1;
		end
		if ~isempty(regexp(AGORAMat{j},presentSpeciesSort{k}))
		    isPresent2 = 1;
		end
	    end

	    if isPresent1==1 && isPresent2==1
	        %file1 = [outputDir1 filesep num2str(i) '_' num2str(j) '_1_' ucrFolderSort '.flux'];
                file2 = [outputDir1 filesep num2str(i) '_' num2str(j) '_2_' ucrFolderSort '.flux'];
                if exist(file2,'file')
	            %disp('HERE')
	            AGORA2Fluxes = importdata(file2);
                    if sum(AGORA2Fluxes~=0) > 50
		        load(AGORAMat{j});
                        modelsSort{end+1} = AGORAModel;
		        fluxesSort{end+1} = AGORA2Fluxes;
                    end
                end
            end


	    isPresent1 = 0;
            isPresent2 = 0;
	    for k=1:length(presentSpeciesMH0007)
		if ~isempty(regexp(AGORAMat{i},presentSpeciesMH0007{k}))
		    isPresent1 = 1;
		end
		if ~isempty(regexp(AGORAMat{j},presentSpeciesMH0007{k}))
		    isPresent2 = 1;
		end
	    end

	    if isPresent1==1 && isPresent2==1
	        %file1 = [outputDir1 filesep num2str(i) '_' num2str(j) '_1_' ucrFolderMH0007 '.flux'];
                file2 = [outputDir1 filesep num2str(i) '_' num2str(j) '_2_' ucrFolderMH0007 '.flux'];
                if exist(file2,'file')
	            %disp('HERE')
	            AGORA2Fluxes = importdata(file2);
                    if sum(AGORA2Fluxes~=0) > 50
		        load(AGORAMat{j});
                        modelsMH0007{end+1} = AGORAModel;
		        fluxesMH0007{end+1} = AGORA2Fluxes;
                    end
                end
            end
	end
    end
end
  
end


ithTotalTermArrSort = [];
ithTotalTermArrMH0007 = [];
compMets = containers.Map;
coopMets = containers.Map;
matrixSort = [];
matrixMH0007 = [];
for z=1:2
    if z==1
	modelsPick = modelsSort;
        fluxesPick = fluxesSort;
    else
	modelsPick = modelsMH0007;
        fluxesPick = fluxesMH0007;
    end
    for i=1:length(modelsPick)
	    i
	ithTotalTerm = 0;
	for j=1:length(modelsPick)
	    modelIth = modelsPick{i};
	    modelJth = modelsPick{j};
	    fluxesIth = fluxesPick{i};
	    fluxesJth = fluxesPick{j};
	    %i
	    %j
	    compTerm = 0;
	    coopTerm = 0;
	    for k=1:length(modelIth.rxns)
		if ~isempty(regexp(modelIth.subSystems{k},'Exchange'))
		    if ~isempty(strcmp(modelJth.rxns,modelIth.rxns{k}))
			correspondK = find(strcmp(modelJth.rxns,modelIth.rxns{j}));
			if fluxesIth(k) < 0 && mean(fluxesJth(correspondK)) < 0
			    compTerm = compTerm + fluxesIth(k)*mean(fluxesJth(correspondK));
                            if ~isKey(compMets,modelIth.rxns{k})
                                compMets(modelIth.rxns{k}) = fluxesIth(k)*mean(fluxesJth(correspondK));
                            else
                                compMets(modelIth.rxnNames{k}) = compMets(modelIth.rxns{k}) + fluxesIth(k)*mean(fluxesJth(correspondK));
                            end
			end
			if (fluxesIth(k) < 0 && mean(fluxesJth(correspondK)) > 0) || (fluxesIth(k) > 0 && mean(fluxesJth(correspondK)) < 0)
			    coopTerm = coopTerm - fluxesIth(k)*mean(fluxesJth(correspondK));
                            if ~isKey(coopMets,modelIth.rxns{k})
                                coopMets(modelIth.rxnNames{k}) = -fluxesIth(k)*mean(fluxesJth(correspondK));
                            else
                                coopMets(modelIth.rxns{k}) = coopMets(modelIth.rxns{k}) - fluxesIth(k)*mean(fluxesJth(correspondK));
                            end
                        end
		    end
		end
	    end
	    %compTerm
	    %coopTerm
	    ithTotalTerm = ithTotalTerm - compTerm + coopTerm;
            if z==1
	        matrixSort(i,j) = -compTerm + coopTerm;
	    else
	      matrixMH0007(i,j) = -compTerm + coopTerm;
	    end
	end
	if z==1
	    ithTotalTermArrSort(i) = ithTotalTerm;
	else
	    ithTotalTermArrMH0007(i) = ithTotalTerm;
	end
    end
end

diffsFake = [ithTotalTermArrSort ithTotalTermArrMH0007];
diffsFake = diffsFake(1:10);
speciesFake = presentSpeciesMH0007(1:length(diffsFake));
[diffsFake sortIdxs] = sort(diffsFake,'descend');
speciesFake = speciesFake(sortIdxs);
writeData({addIdxStrings(speciesFake),diffsFake},[transferDir filesep 'speciesInfluenceDiffs.txt'],'\t',{'species','diff'});

keysComp = keys(compMets);
valuesComp = [];
for i=1:length(keysComp);
valuesComp(i) = compMets(keysComp{i});
end
[valuesComp sortIdxs] = sort(valuesComp,'descend');
keysComp = keysComp(sortIdxs);
writeData({addIdxStrings(keysComp),valuesComp},[transferDir filesep 'metsInfluenceComp.txt'],'\t',{'metname','influence'})

keysCoop = keys(coopMets);
valuesCoop = [];
for i=1:length(keysCoop);
valuesCoop(i) = coopMets(keysCoop{i});
end
[valuesCoop sortIdxs] = sort(valuesCoop,'descend');
keysCoop = keysCoop(sortIdxs);
writeData({addIdxStrings(keysCoop),valuesCoop},[transferDir filesep 'metsInfluenceCoop.txt'],'\t',{'metname','influence'})

rxnsToFluxesSort = containers.Map;
rxnsToSubsSort = containers.Map;
		    for i=1:length(modelsSort)
			    AGORAModel = modelsSort{i};
AGORAFluxes = fluxesSort{i};
for j=1:length(AGORAModel.rxns)
if isKey(rxnsToFluxesSort,AGORAModel.rxnNames{j})
	rxnsToFluxesSort(AGORAModel.rxnNames{j}) = rxnsToFluxesSort(AGORAModel.rxnNames{j}) + abs(AGORAFluxes(j));
 else
   rxnsToFluxesSort(AGORAModel.rxnNames{j}) = abs(AGORAFluxes(j));
temp = AGORAModel.subSystems{j};
rxnsToSubsSort(AGORAModel.rxnNames{j}) = temp{1};
end
end
end

keysSort = keys(rxnsToFluxesSort);
subsSort = {};
valuesSort = [];
for i=1:length(keysSort);
valuesSort(i) = rxnsToFluxesSort(keysSort{i});
subsSort{i} = rxnsToSubsSort(keysSort{i});
end
[valuesSort sortIdxs] = sort(valuesSort);
keysSort = keysSort(sortIdxs);
subsSort = subsSort(sortIdxs);

		    rxnsToFluxesMH0007 = containers.Map;
rxnsToSubsMH0007 = containers.Map;
		    for i=1:length(modelsMH0007)
			    AGORAModel = modelsMH0007{i};
AGORAFluxes = fluxesMH0007{i};
for j=1:length(AGORAModel.rxns)
if isKey(rxnsToFluxesMH0007,AGORAModel.rxnNames{j})
	rxnsToFluxesMH0007(AGORAModel.rxnNames{j}) = rxnsToFluxesMH0007(AGORAModel.rxnNames{j}) + abs(AGORAFluxes(j));
 else
   rxnsToFluxesMH0007(AGORAModel.rxnNames{j}) = abs(AGORAFluxes(j));
temp = AGORAModel.subSystems{j};
rxnsToSubsMH0007(AGORAModel.rxnNames{j}) = temp{1};
end
end
end

keysMH0007 = keys(rxnsToFluxesMH0007);
valuesMH0007 = [];
subsMH0007 = {};
for i=1:length(keysMH0007);
valuesMH0007(i) = rxnsToFluxesMH0007(keysMH0007{i});
subsMH0007{i} = rxnsToSubsMH0007(keysMH0007{i});
end
[valuesMH0007 sortIdxs] = sort(valuesMH0007);
keysMH0007 = keysMH0007(sortIdxs);
subsMH0007 = subsMH0007(sortIdxs);

keysDiff = {};
valuesDiff = [];
subsDiff = {};
for i=1:length(keysSort)
	if isKey(rxnsToFluxesMH0007,keysSort{i})
	keysDiff{end+1} = keysSort{i};
valuesDiff(end+1) = valuesSort(i)-rxnsToFluxesMH0007(keysSort{i});
subsDiff{end+1} = subsSort{i};
end
end

[valuesDiff sortIdxs] = sort(valuesDiff,'descend');
keysDiff = keysDiff(sortIdxs);
subsDiff = subsDiff(sortIdxs);

keysDiffExch = {};
valuesDiffExch = [];
for i=1:length(keysDiff)
	if ~isempty(regexp(subsDiff{i},'Exchange'))
	keysDiffExch{end+1} = keysDiff{i};
valuesDiffExch(end+1) = valuesDiff(i);
  end
  end

  for i=1:length(keysDiffExch)
	  if ~isempty(regexp(keysDiffExch{i},'NAD|icotinamide|lutathione|scorb|ocopherol|arnosine|ilirubin|itamin'))
	  disp(keysDiffExch{i})
	  disp(valuesDiffExch(i))
	  end
	  end

	  writeData({addIdxStrings(keysDiffExch),valuesDiffExch},[transferDir filesep 'diffExchRxnsAGORA.txt'],'\t',{'rxnname','difflux'})
	  	  writeData({addIdxStrings(keysDiff),valuesDiff},[transferDir filesep 'diffRxnsAGORA.txt'],'\t',{'rxnname','difflux'})
