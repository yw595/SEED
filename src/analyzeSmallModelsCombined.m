diabetesctrl = {'BGI-06A','N044A','SZEY-75A'};
diabetestype2nomet = {'NG-5636_551','DOM024','DLM001'};
HMP2Normal = {'206703','206704','206700'};
HMP2IBD = {'206701','206708','206709'};
MHnormal = {'MH0005','MH0006','MH0008'};
MHobese = {'MH0001','MH0002','MH0003'};
grouparr = {diabetesctrl,diabetestype2nomet,HMP2Normal,HMP2IBD,MHnormal,MHobese};
grouplabelarr = {'diabetesctrl','diabetestype2nomet','HMP2Normal','HMP2IBD','MHnormal','MHobese'};
centsarr = {};
for z1=1:length(grouparr)
	 centsarr{z1} = {};
    for z=1:length(grouparr{z1})
	    interactMat2 = importdata(['/mnt/vdb/home/ubuntu2/interactMatTemp' grouparr{z1}{z} '.txt']);
cents = betweenness_centrality(abs(sparse(interactMat2)));
centsarr{z1}{z} = cents;
end
end

sumcentsarr = {};
for z1=1:length(centsarr)
    sumcents = centsarr{z1}{1};
    for i=2:length(centsarr{z1})
	sumcents = sumcents + centsarr{z1}{i};
    end
    sumcents = sumcents/length(centsarr{z1});
    sumcentsarr{z1} = sumcents;
    if mod(z1,2)==0
        sumcentsdiff = sumcentsarr{z1}-sumcentsarr{z1-1};
        [~,sortIdxs] = sort(abs(sumcentsdiff),'descend');
        diffarr = {};
        namearr = {};
        for k=1:10
	    diffarr{k} = sumcentsdiff(sortIdxs(k));
            name = AGORAMat{sortIdxs(k)};
            name = strsplit(name,'/');
            namearr{k} = name{end};
	end
	diffarr
	namearr
    end
end
