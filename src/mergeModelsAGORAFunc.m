function returnName = mergeModelsAGORAFunc(model,idx,isRxn,addNum)

if isRxn
    name = model.rxns{idx};
else
    name = model.mets{idx};
end

if isRxn
    if ~isempty(regexp('xchange',model.subSystems{idx}))
        returnName = name;
    else
        returnName = [name '_' addNum];
    end
else
    if ~isempty(regexp('\[e\]',name))
        returnName = name;
    else
        returnName = [name '_' addNum];
    end
end

end
