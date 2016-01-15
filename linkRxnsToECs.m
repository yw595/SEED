if 1
if 0
outputDir = '/mnt/extra/SEED';
filenames = dir(outputDir);
rxnsToECs = containers.Map;
ECsToRxns = containers.Map;
allRxnNames = containers.Map;
rxnsToMets = containers.Map;
rxnsToCoeffs = containers.Map;
end
for i=810:length(filenames)
    if ~isempty(regexp(filenames(i).name,'.tsv')) %&& count < 1
        %count=count+1;
        modelName = filenames(i).name; modelName = modelName(1:regexp(modelName,'.tsv')-1);
        if 1
        status=system(sprintf(['/home/ubuntu/MATLAB/SEED/makeTemp.sh %s %s'],[outputDir filesep modelName '.tsv'],[outputDir filesep modelName '.temp']));
        disp(modelName);
        FI = fopen([outputDir filesep modelName '.temp']);
        line = fgetl(FI);
        count = 0;
        while line ~= -1
            if mod(count,25)==0
                disp(modelName)
                disp(count)
                disp(line)
            end
            count = count+1;
            [regex1 regex2] = regexp(line,'(\d|-)+\.(\d|-)+\.(\d|-)+\.(\d|-)+');
            if ~isempty(regex1) 
                ECNums = arrayfun(@(x,y) line(x:y), regex1,regex2, 'UniformOutput',0);
                words = strsplit(line,',');
                rxnName = words{1};
                if strcmp(rxnName,'rxn00792')
                    %disp(line)
                    %disp(regex1)
                    %disp(ECNums)
                    if isKey(rxnsToECs,rxnName)
                        %disp(rxnsToECs(rxnName))
                    end
                end
                if ~isKey(rxnsToECs,rxnName)
                    rxnsToECs(rxnName) = ECNums;
                else
                    temp = rxnsToECs(rxnName);
                    if ~iscell(temp)
                        temp = {temp};
                    end
                    for j=1:length(ECNums)
                        temp{end+1} = ECNums{j};
                    end
                    temp = unique(temp);
                    rxnsToECs(rxnName) = temp;
                    % if ~all(ismember(ECNums,rxnsToECs(rxnName)))
                    %     disp('WARNING: INCONSISTENT EC NUMBERS');
                    %     temp = unique([ECNums rxnsToECs(rxnName)]);
                    %     disp(temp);
                    %     rxnsToECs(rxnName) = temp;
                    % end 
                end
                for j=1:length(ECNums)
                    if ~isKey(ECsToRxns,ECNums{j})
                        ECsToRxns(ECNums{j}) = rxnName;
                    else
                        temp = ECsToRxns(ECNums{j});
                        if ~iscell(temp) 
                            temp = {temp};
                        end
                        %for k=1:length(ECNums)
                            temp{end+1} = rxnName;
                            %end
                        temp = unique(temp);
                        ECsToRxns(ECNums{j}) = temp;
                        % if ~all(ismember(ECNums,rxnsToECs(rxnName)))
                        %     disp('WARNING: INCONSISTENT EC NUMBERS');
                        %     temp = unique([ECNums rxnsToECs(rxnName)]);
                        %     disp(temp);
                        %     rxnsToECs(rxnName) = temp;
                        % end 
                    end
                end
            end
            line = fgetl(FI);
        end
        %disp(count)
        fclose(FI);
        end
        
        if count > 10
            %if ~exist(strrep(modelName,'.','_'),'var')
                modelTemp = readCbModel([outputDir filesep modelName '.xml']);
                %eval([strrep(modelName,'.','_') ' = modelTemp;']);
                %end
                %modelTemp = eval(strrep(modelName,'.','_'));
            %notInAll = find(~isKey(allRxnNames,modelTemp.rxns));
            %notInAll = find(~ismember(modelTemp.rxns,));
            for j=1:length(modelTemp.rxns)
                %notMets = modelTemp.mets(modelTemp.S(:,notInAll(j))~=0);
                %notCoeffs = modelTemp.S(:,modelTemp.S(:,notInAll(j))~=0);
                rxnsToMets(modelTemp.rxns{j}) = modelTemp.mets(modelTemp.S(:,j)~=0);
                rxnsToCoeffs(modelTemp.rxns{j}) = modelTemp.S(modelTemp.S(:,j)~=0,j);
            end           
        end
    end
end
end

if 0
FI = fopen([outputDir filesep 'GreenblumECs.txt']);
GreenblumEC = textscan(FI,'%s\n');
fclose(FI);
GreenblumEC = GreenblumEC{1};
bigModel = struct(); bigModel.S = []; bigModel.rxns = {}; bigModel.mets = {};
ECs = keys(ECsToRxns);
for i=1:length(ECs)
    if any(strcmp(ECs{i},GreenblumEC))
        rxns = ECsToRxns(ECs{i});
        if ~iscell(rxns)
            rxns = {rxns};
        end
        for j=1:length(rxns)
            if ~any(strcmp(bigModel.rxns,rxns{j}))
                bigModel.rxns{end+1} = rxns{j};
                mets = rxnsToMets(rxns{j});
                coeffs = rxnsToCoeffs(rxns{j});
                for k=1:length(mets)
                    if ~any(strcmp(bigModel.mets,mets{k}))
                        bigModel.mets{end+1} = mets{k};
                    end
                    bigModel.S(strcmp(bigModel.mets,mets{k}),strcmp(bigModel.rxns,rxns{j})) = coeffs(k);
                end
            end
        end
    end
end
end
            
% ECsToRxns = containers.Map;
% rxnsToECsKeys = keys(rxnsToECs);
% for i=1:length(rxnsToECsKeys)
%     temp1 = rxnsToECsKeys{i};
%     temp2 = rxnsToECs(rxnsToECsKeys{i}); %temp2 = temp2{1};
%     for j = 1:length(temp2)
%         temp3 = temp2{j};
%         %disp(i); disp(temp1); disp(temp3);
%         if isKey(ECsToRxns,temp3)
%             ECsToRxns(temp3) = temp1;
%         else
%             if ~all(ismember(temp1,ECsToRxns(temp3)))
%                 disp('WARNING: INCONSISTENT RXN NAMES');
%                 temp = unique([temp1 ECsToRxns(temp3)]);
%                 disp(temp);
%                 ECsToRxns(temp3) = temp;
%             end
%         end
%     end
% end