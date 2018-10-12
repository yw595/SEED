checkRandom = 0;
if checkRandom
    checkRandomArr = [];
    zEnd = 200;
else
    zEnd = 1;
end
for z=1:zEnd
    idxsToIterate = 1:length(AGORAModelArr);
    if checkRandom
        idxsToIterate = randperm(length(AGORAModelArr));
    end
    tobemerged = AGORAModelArr{idxsToIterate(1)};
    for i=1:length(idxsToIterate)
	checkRandomArr(z,i) = length(tobemerged.rxns);
        if ~isempty(AGORAModelArr{idxsToIterate(i)})
	    tobemerged = mergeModels(tobemerged,AGORAModelArr{idxsToIterate(i)});
	    disp(z)
	    disp(i)
	    disp(length(tobemerged.rxns))
	end
    end
    checkRandomArr(z,i+1) = length(tobemerged.rxns);
end
