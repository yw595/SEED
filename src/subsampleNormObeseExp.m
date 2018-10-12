commonKeys = intersect(rxnsToExpressNormKeys,keys(rxnsToExpressObese));
randObese = {};
randNorm = {};
for i=1:10
	selectIdxs = randperm(length(commonKeys));
selectIdxs = selectIdxs(1:1000);
ithRandNorm = [];
ithRandObese = [];
for j=1:1000
	ithRandNorm(j) = rxnsToExpressNorm(commonKeys{selectIdxs(j)});
	ithRandObese(j) = rxnsToExpressObese(commonKeys{selectIdxs(j)});
end
ithRandNorm = sort(ithRandNorm); ithRandNorm = ithRandNorm(101:900);
ithRandObese = sort(ithRandObese); ithRandObese = ithRandObese(101:900);
randNorm{i} = ithRandNorm;
randObese{i} = ithRandObese;
end
writeFI = fopen('/mnt/vdb/home/ubuntu2/otu_table_rand.txt','w');
writeFI2 = fopen('/mnt/vdb/home/ubuntu2/metadata_rand.txt','w');
for i=1:length(randNorm)
	fprintf(writeFI,'%f\t',randNorm{i});
fprintf(writeFI,'\n');
fprintf(writeFI2,'1\n');
end
for i=1:length(randObese)
	fprintf(writeFI,'%f\t',randObese{i});
fprintf(writeFI,'\n');
fprintf(writeFI2,'0\n');
end
fclose(writeFI);
fclose(writeFI2);
