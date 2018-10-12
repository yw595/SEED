outputDir1 = [outputDir filesep 'simulateSmallModelsSeparatePicrustAGORA'];

countOxFermExpr = 0;

if 0
moreThanTwenty = [];
exactCount = [];
for i=1:length(AGORAMat)
    load(AGORAMat{i})
    i
    [~,~,oxRxnIdxs,fermentRxnIdxs] = measureOxFermFunc(AGORAModel);
    exactCount(i,1) = length(oxRxnIdxs);
    exactCount(i,2) = length(fermentRxnIdxs);
    if length(oxRxnIdxs)>20
         moreThanTwenty(i) = 1;
     end
 end
 end


 MHnormal = importdata('/mnt/xvdf/home/ubuntu2/MHnormal.txt');
 MHobese = importdata('/mnt/xvdf/home/ubuntu2/MHobese.txt');
 if 0
 MHnormox = zeros(115,1);
 MHobeseox = zeros(115,1);
 MHnormoxcount = zeros(115,1);
 MHobeseoxcount = zeros(115,1);
 for i=1:115%length(ucrFolders)
 	i
     for j=1:size(oxFermMat,2)
 	%if moreThanTwenty(j)==1
 	    MHID = ucrFolders{i};
             MHID = MHID(7:end);
             if sum(strcmp(MHnormal,MHID))~=0
 	        MHnormox(i) = MHnormox(i)+oxFermMat(i,j,1);
 	        if exist([inputDir filesep 'MGMData' filesep ucrFolders{i} filesep 'speciesSep' num2str(j) '.modelexpr'],'file')
                     MHnormoxcount(i) = MHnormoxcount(i)+exactCount(j,1);
                 end
             end
             if sum(strcmp(MHobese,MHID))~=0
 	        MHobeseox(i) = MHobeseox(i)+oxFermMat(i,j,1);
 	        if exist([inputDir filesep 'MGMData' filesep ucrFolders{i} filesep 'speciesSep' num2str(j) '.modelexpr'],'file')
                     MHobeseoxcount(i) = MHobeseoxcount(i)+exactCount(j,1);
                 end
             end
         %end
     end
 end
 	    grouplabels = {};
 for i=1:length(MHnormox)
 	grouplabels{end+1} = 'normal';
 end
 for i=1:length(MHobeseox)
 	grouplabels{end+1} = 'obese';
 end
 %writeData({[MHnormox;MHobeseox],grouplabels},'/mnt/xvdf/home/ubuntu2/moreThanTwentyOxExpr.txt','\t',{'expr','grouplabel'});
 end

 if 0
 if ~countOxFermExpr
     fluxesMat = {};
 end
 if countOxFermExpr
     oxFermMat = [];
 end
 for z=1:length(ucrFolders)
     system(['python /mnt/xvdf/home/ubuntu2/pickPresent.py ' inputDir filesep 'MGMData' filesep ucrFolders{z} filesep 'normalized_otus.tsv 800']);
     presentSpecies = textscan(fopen('/mnt/xvdf/home/ubuntu2/pickPresent800.txt'),'%s','Delimiter','\n','HeaderLines',0);
     presentSpecies = presentSpecies{1};
     presentIdxs = zeros(length(AGORAMat),1);
     for i=1:length(AGORAMat)
 	    if ~isempty(AGORAMat{i})
 	        isPresent = 0;
 		for k=1:length(presentSpecies)
 		    if ~isempty(regexp(AGORAMat{i},presentSpecies{k}))
 			presentIdxs(i) = 1;
                    end
                end
            end
     end
     if countOxFermExpr
 	for i=1:length(AGORAMat)
 	    i
 	    z
 	    if exist([inputDir filesep 'MGMData' filesep ucrFolders{z} filesep 'speciesSep' num2str(i) '.modelexpr'],'file')
 		load(AGORAMat{z})
 		AGORAFluxes = importdata([inputDir filesep 'MGMData' filesep ucrFolders{z} filesep 'speciesSep' num2str(i) '.modelexpr']);
 		if length(AGORAModel.rxns)<=length(AGORAFluxes)
 		    [totalOx,totalFerm] = measureOxFermFunc(AGORAModel,AGORAFluxes);
 		    oxFermMat(z,i,1) = totalOx;
 		    oxFermMat(z,i,2) = totalFerm;
 		end
 	    end
 	end
     else
     for i=1:length(AGORAMat)
         i
 	z
         for j=1:length(AGORAMat)
 	    if i~=j && ~isempty(AGORAMat{i}) && ~isempty(AGORAMat{j})
 	        isPresent1 = presentIdxs(i);
                 isPresent2 = presentIdxs(j);
 if 0
 		for k=1:length(presentSpecies)
 		    if ~isempty(regexp(AGORAMat{i},presentSpecies{k}))
 			isPresent1 = 1;
 		    end
 		    if ~isempty(regexp(AGORAMat{j},presentSpecies{k}))
 		        isPresent2 = 1;
 		    end
 		end
 end
 	        if isPresent1==1 && isPresent2==1
 		      file1 = [outputDir1 filesep num2str(i) '_' num2str(j) '_1_' ucrFolders{z} '.flux'];
                     %file2 = [outputDir1 filesep num2str(i) '_' num2str(j) '_2_' ucrFolderSort '.flux'];
                     if exist(file1,'file')
                         AGORAFluxes = importdata(file1);
                         fluxesMat{z,i} = AGORAFluxes;
                     end
 		    %if exist(file2,'file')
                     %    AGORAFluxes = importdata(file2);
                     %    fluxesMat{z,i,j,2} = AGORAFluxes;
                     %end
 		end
 	    end
 	end
     end
     end
 end
 end

 		    if 0
 		    AGORAModelArr = {};
 for i=1:length(AGORAMat)
 	i
 	load(AGORAMat{i});
 AGORAModelArr{i} = AGORAModel;
 end
 end

 if 0
 allInteractMat = [];
 for z=1:10
     for i=1:size(fluxesMat,2)
 	    z
 	    i
 	for j=1:size(fluxesMat,2)
 		if ~isempty(fluxesMat{z,i}) && ~isempty(fluxesMat{z,j})
 		[compTerm coopTerm coopBasicTerm compBasicTerm] = metInteractFuncTemp2(AGORAModelArr{i},AGORAModelArr{j},fluxesMat{z,i},fluxesMat{z,j});
 	    allInteractMat(z,i,j,1) = compTerm;
 	    allInteractMat(z,i,j,2) = coopTerm;
 	    allInteractMat(z,i,j,3) = coopBasicTerm;
 	    allInteractMat(z,i,j,4) = compBasicTerm;
 end
 	end
     end
 end
 end

 if 0
 allOxRxnIdxs = {};
 for i=1:length(AGORAMat)
     if exist(AGORAMat{i},'file')
 	i
         load(AGORAMat{i});
         AGORAModelArr{i} = AGORAModel;
         [~,~,oxRxnIdxs,~] = measureOxFermFunc(AGORAModel);
         allOxRxnIdxs{i} = oxRxnIdxs;
     end
 end
    
 end

 if 0
 oxExprMat = [];
for i=1:115
    %i=ceil(z/length(AGORAModelArr));
    %j=mod(z,length(AGORAModelArr));
    for j=1:length(AGORAModelArr)
	i
	j
	exprFile = [inputDir filesep 'MGMData' filesep ucrFolders{i} filesep 'speciesSep' num2str(j) '.modelexpr'];
        if exist(exprFile,'file')
            %[~,~,oxRxnIdxs,~] = measureOxFermFunc(AGORAModelArr{i});
oxRxnIdxs = allOxRxnIdxs{j};
            exprData = importdata(exprFile);
            if max(oxRxnIdxs)<=length(exprData)
                oxExprMat(i,j) = sum(exprData(oxRxnIdxs));
            end
        end
    end
end
end

allOxECs = {};
for i=1:length(AGORAModelArr)
	i
	allOxECs = union(allOxECs,AGORAModelArr{i}.rxnECNumbers(allOxRxnIdxs{i}));
end

allOxRxnCounts = [];
 for z=1:length(ucrFolders)
     system(['python /mnt/xvdf/home/ubuntu2/pickPresent.py ' inputDir filesep 'MGMData' filesep ucrFolders{z} filesep 'normalized_otus.tsv 800']);
     presentSpecies = textscan(fopen('/mnt/xvdf/home/ubuntu2/pickPresent800.txt'),'%s','Delimiter','\n','HeaderLines',0);
     presentSpecies = presentSpecies{1};
     presentIdxs = zeros(length(AGORAMat),1);
     for i=1:length(AGORAMat)
 	    if ~isempty(AGORAMat{i})
 	        isPresent = 0;
 		for k=1:length(presentSpecies)
 		    if ~isempty(regexp(AGORAMat{i},presentSpecies{k}))
 			presentIdxs(i) = 1;
                    end
                end
            end
     end
z
		    allOxRxnCounts(z) = 0;
		    for i=1:length(allOxRxnIdxs)
			    if presentIdxs(i)==1
			    allOxRxnCounts(z) = allOxRxnCounts(z) + length(allOxRxnIdxs{i});
end
end
 end
		    
interactMat = [];
for z=1:10
interactMatZth = allInteractMat(z,:,:,2)-allInteractMat(z,:,:,1);
interactMat(z,:) = interactMatZth(:);
end
interactMat = interactMat/1000;
eigs = [];
for z=1:10
    eigZth = eig(squeeze(allInteractMat(z,:,:,2)-allInteractMat(z,:,:,1)));
    eigZth = real(eigZth); [~, sortIdxs] = sort(abs(eigZth),'descend'); eigZth = eigZth(sortIdxs);
    eigs(:,z) = eigZth*-1;
end
eiglabels = {'eig01','eig02','eig03','eig04','eig05','eig06','eig07','eig08','eig09','eig10'};
writeData({eiglabels,mean(eigs(1:10,[5 6 8]),2),max(eigs(1:10,[1 2 3 4 7 9 10]),[],2),min(eigs(1:10,[1 2 3 4 7 9 10]),[],2)},'/mnt/xvdf/home/ubuntu2/AGORAEigsRenewObese.txt','\t',{'eigidx','meaneig','upper','lower'});
writeData({eiglabels,mean(eigs(1:10,[5 6 8]),2),max(eigs(1:10,[5 6 8]),[],2),min(eigs(1:10,[5 6 8]),[],2)},'/mnt/xvdf/home/ubuntu2/AGORAEigsRenewNormal.txt','\t',{'eigidx','meaneig','upper','lower'});
histsteps = min(min(interactMat)):(max(max(interactMat))-min(min(interactMat)))/100:max(max(interactMat));
