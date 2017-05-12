sum(sign(hadzadiffs')==sign(exchangeDiffs(:,1)))
sum(sign(hadzadiffs')==sign(exchangeDiffs(:,2)))
sum(sign(hadzadiffs')==sign(exchangeDiffs(:,3)))
sum(sign(hadzadiffs')==sign(exchangeDiffs(:,4)))


corr(hadzadiffs',exchangeDiffs(:,1),'type','Spearman')
corr(hadzadiffs',exchangeDiffs(:,2),'type','Spearman')
corr(hadzadiffs',exchangeDiffs(:,3),'type','Spearman')
corr(hadzadiffs',exchangeDiffs(:,4),'type','Spearman')

writeData({hadzadiffs,exchangeDiffs(:,1)},'/home/fs01/yw595/expVsSimFluxFBA.txt','\t',{'expchange','simchange'});
writeData({hadzadiffs,exchangeDiffs(:,2)},'/home/fs01/yw595/expVsSimFluxEFlux.txt','\t',{'expchange','simchange'});
writeData({hadzadiffs,exchangeDiffs(:,3)},'/home/fs01/yw595/expVsSimFluxFALCON.txt','\t',{'expchange','simchange'});
writeData({hadzadiffs,exchangeDiffs(:,4)},'/home/fs01/yw595/expVsSimFluxGXFBA.txt','\t',{'expchange','simchange'});

corr(bigModelReconc.express,abs(convertArr(1,:,1))','type','Spearman')
corr(bigModelReconc.express,abs(convertArr(1,:,2))','type','Spearman')
corr(bigModelReconc.express,abs(convertArr(1,:,3))','type','Spearman')
corr(bigModelReconc.express,abs(convertArr(1,:,4))','type','Spearman')

writeData({bigModelReconc.express,abs(convertArr(1,:,1))},'/home/fs01/yw595/expressionVsFluxFBA.txt','\t',{'expression','flux'});
writeData({bigModelReconc.express,abs(convertArr(1,:,2))},'/home/fs01/yw595/expressionVsFluxEFlux.txt','\t',{'expression','flux'});
writeData({bigModelReconc.express,abs(convertArr(1,:,3))},'/home/fs01/yw595/expressionVsFluxFALCON.txt','\t',{'expression','flux'});
writeData({bigModelReconc.express,abs(convertArr(1,:,4))},'/home/fs01/yw595/expressionVsFluxGXFBA.txt','\t',{'expression','flux'});