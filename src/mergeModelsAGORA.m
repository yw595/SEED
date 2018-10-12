function returnModel = mergeModelsAGORA(AGORAModel1Temp,AGORAModel2Temp)

%AGORAModel1Temp = AGORAModel2;
%AGORAModel2Temp = AGORAModel2;
AGORAModel1Temp.rxns = arrayfun(@(x) mergeModelsAGORAFunc(AGORAModel1Temp,x,1,'1'), 1:length(AGORAModel1Temp.rxns), 'UniformOutput', 0)';
AGORAModel2Temp.rxns = arrayfun(@(x) mergeModelsAGORAFunc(AGORAModel2Temp,x,1,'2'), 1:length(AGORAModel2Temp.rxns), 'UniformOutput', 0)';
AGORAModel1Temp.mets = arrayfun(@(x) mergeModelsAGORAFunc(AGORAModel1Temp,x,0,'1'), 1:length(AGORAModel1Temp.mets), 'UniformOutput', 0)';
AGORAModel2Temp.mets = arrayfun(@(x) mergeModelsAGORAFunc(AGORAModel2Temp,x,0,'2'), 1:length(AGORAModel2Temp.mets), 'UniformOutput', 0)';
returnModel = mergeTwoModels(AGORAModel1Temp,AGORAModel2Temp);
