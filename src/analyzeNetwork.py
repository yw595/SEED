import networkx as nx
import os
import numpy as np
from bootstrapAGORAIBD import writeData

diabetesctrl = ['BGI-06A','N044A','SZEY-75A']
diabetestype2nomet = ['NG-5636_551','DOM024','DLM001']
HMP2Normal = ['206703','206704','206700']
HMP2IBD = ['206701','206708','206709']
MHnormal = ['MH0005','MH0006','MH0008']
MHobese = ['MH0001','MH0002','MH0003']
grouparr = [diabetesctrl,diabetestype2nomet,HMP2Normal,HMP2IBD,MHnormal,MHobese]
grouplabelarr = ['diabetesctrl','diabetestype2nomet','HMP2Normal','HMP2IBD','MHnormal','MHobese']
groupmap = {}
for i in range(len(grouplabelarr)):
    groupmap[grouplabelarr[i]] = []
files = os.listdir('/mnt/vdb/home/ubuntu2/')
AGORAMat = []
AGORAMatFI = open('/mnt/vdb/home/ubuntu2/AGORAMat.txt')
for line in AGORAMatFI:
    AGORAMat.append(line.strip().split('/')[10])
AGORAMatFI.close()
for z1 in range(len(grouparr)):
    for z in range(len(grouparr[z1])):
        prefix = grouparr[z1][z]
        afile = 'interactMatTemp'+prefix+'.txt'
        inFI = open('/mnt/vdb/home/ubuntu2/'+afile)
        ithIdx = 0
        G = nx.Graph()
        for line in inFI:
            words = line.split('\t')
            for j in range(len(words)):
                if float(words[j])!=0:
                    G.add_edge(str(ithIdx),str(j),weight=abs(float(words[j])))
            ithIdx += 1
        print('HERE')
        out = nx.betweenness_centrality(G)
        print(afile)
        print(sum(out.values()))
        print(sum(np.array(out.values())!=0))
        valuesarr = [0 for k in range(773)]
        for key in out:
            valuesarr[int(key)] = out[key]
        groupmap[grouplabelarr[z1]].append(valuesarr)
        inFI.close()

    meanarr = []
    threelists = groupmap[grouplabelarr[z1]]
    for k in range(len(threelists[0])):
        meanarr.append((threelists[0][k]+threelists[1][k]+threelists[2][k])/3)
    writeData([AGORAMat,meanarr],'/mnt/vdb/home/ubuntu2/cents'+grouplabelarr[z1]+'.txt',delimiter='\t')
