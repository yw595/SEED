fbaOxMat = [];
fbaFermentMat = [];
for i=1:length(AGORAMat)
    load(AGORAMat{i})
    for j=1:length(AGORAModel.rxns)
        metIdxs = find(AGORAModel.S(:,j)~=0);
        metKEGGs = AGORAModel.metKEGGID(metIdxs);
if length(metKEGGs)==1 && (strcmp(metKEGGs,'C00007'))% || strcmp(metKEGGs,'C00186') || strcmp(metKEGGs,'C00469') || strcmp(metKEGGs,'C00246') || strcmp(metKEGGs,'C00033') || strcmp(metKEGGs,'C00163') || strcmp(metKEGGs,'C00084'))
            i
	    %disp(AGORAModel.rxns{j})
	    %disp(AGORAModel.lb(j))
	    %disp(AGORAModel.ub(j))
	    AGORAModel = assignSortedBiom(AGORAModel);
soln = optimizeCbModel(AGORAModel,1);
soln = soln.x;
[returnArr1 returnArr2] = measureOxFermFunc(AGORAModel,soln);
fbaOxMat(i,1) = returnArr1;
fbaFermentMat(i,1) = returnArr2;
AGORAModel.lb(j)=-1000;
soln = optimizeCbModel(AGORAModel,1);
soln = soln.x;
[returnArr1 returnArr2] = measureOxFermFunc(AGORAModel,soln);
fbaOxMat(i,2) = returnArr1;
fbaFermentMat(i,2) = returnArr2;
disp(fbaOxMat(i,1))
disp(fbaOxMat(i,2))
disp(fbaFermentMat(i,1))
disp(fbaFermentMat(i,2))
        end
    end
end
