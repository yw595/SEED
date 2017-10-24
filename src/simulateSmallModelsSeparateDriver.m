if 1
modelLimit = 212;
iterLimit = 1;
FBAGapMatrixAvg = zeros(modelLimit,modelLimit);
FBAGapMatrixLims = [];
for z=1:iterLimit
	FBAGapMatrixSingle = simulateSmallModelsSeparateFactored(modelLimit,modelNames,modelNamesToModels,allBiomassRates,allBiomassDists);
    FBAGapMatrixAvg = FBAGapMatrixAvg+FBAGapMatrixSingle;
    if z==1
        FBAGapMatrixLims(:,:,1) = FBAGapMatrixSingle;
        FBAGapMatrixLims(:,:,2) = FBAGapMatrixSingle;
    else
        FBAGapMatrixLims(:,:,1) = min(FBAGapMatrixLims(:,:,1),FBAGapMatrixSingle);
        FBAGapMatrixLims(:,:,2) = max(FBAGapMatrixLims(:,:,2),FBAGapMatrixSingle);
    end
end
FBAGapMatrixAvg = FBAGapMatrixAvg/iterLimit;
end

modelArr1 = {};
modelArr2 = {};
multipleArr = [];
disp(modelLimit)
denomsArr = [];
for i=1:modelLimit
    for j=1:modelLimit
	denomsArr(i,j) = max(abs(FBAGapMatrixAvg(i,j)-FBAGapMatrixLims(i,j,1)),abs(FBAGapMatrixAvg(i,j)-FBAGapMatrixLims(i,j,2)));
    end
end
denomsArr(denomsArr==0) = min(min(denomsArr(denomsArr~=0)));
for i=1:modelLimit
    for j=1:modelLimit
	i
	j
	disp(length(multipleArr))
	modelArr1{end+1} = modelNamesShort{i};
        modelArr2{end+1} = modelNamesShort{j};
        denom = denomsArr(i,j);
        multipleArr(end+1) = abs(FBAGapMatrixAvg(i,j)/denom);
    end
end
%writeData({modelArr1,modelArr2,multipleArr},[transferDir filesep 'FBAGapMultiple.txt'],'\t',{'species1','species2','fbagapmultiple'});

    
