mfunction newModel = addMetSubsystems(origModel)

newModel = origModel;
newModel.metSubsystems = {};
for i=1:length(newModel.metKEGGs)
    correspondingSubsystems = newModel.subSystems(newModel.S(i,:)~=0);
subsToNums = containers.Map;
largestSub = '';
largestNum = 0;
for j=1:length(correspondingSubsystems)
    if largestNum==0
        largestNum = 1;
        largestSub = correspondingSubsystems{j};
    end
    if ~isKey(subsToNums,correspondingSubsystems{j})
	subsToNums(correspondingSubsystems{j}) = 1;
    else
	subsToNums(correspondingSubsystems{j}) = subsToNums(correspondingSubsystems{j})+1;
        if subsToNums(correspondingSubsystems{j}) > largestNum
            largestNum = subsToNums(correspondingSubsystems{j});
            largestSub = correspondingSubsystems{j};
        end
    end
end
newModel.metSubsystems{i} = largestSub;
end
