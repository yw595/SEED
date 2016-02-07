HMPECs = unique([splittest{:}]);
bigModelECs = unique([keys(ECsToRxnsAccum) keys(ECsToRxnsTable)]);

xvals = 1:2; yvals = [length(intersect(HMPECs,bigModelECs)) length(setdiff(HMPECs,bigModelECs));  length(intersect(HMPECs,bigModelECs)) length(setdiff(bigModelECs,HMPECs))];
titleString = 'compareBigModelHMPRef';
makeBar(xvals,yvals,titleString,outputDir,'ylabelString','Num ECs','xlabels',{'HMPRef','bigModel'},'legendLabels',{'HMPRef','bigModel'});