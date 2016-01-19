function connMatrix = makeConnMatrix(model)

connMatrix = zeros(length(model.rxns)+length(model.mets),length(model.rxns)+length(model.mets));
for i=1:length(model.rxns)
    metIdxs = find(model.S(:,i)~=0);
    for j=1:length(metIdxs)
        connMatrix(i,length(model.rxns)+j)=1;
        connMatrix(length(model.rxns)+j,i)=1;
    end
end

end