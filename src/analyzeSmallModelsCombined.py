import os
import numpy as np
import copy

inputDir = '/mnt/vdb/home/ubuntu2/MATLAB/SEED/output/simulateSmallModelsCombined/'

mergedFI = open('/mnt/vdb/home/ubuntu2/tobemerged2.txt')
mergedFI.readline()
rxnstosubs = {}
for line in mergedFI:
    words = line.strip().split('\t')
    if len(words)==3:
        rxnstosubs[words[0]] = words[2]
mergedFI.close()

diabetesctrl = ['BGI-06A','N044A','SZEY-75A']
diabetestype2nomet = ['NG-5636_551','DOM024','DLM001']
HMP2Normal = ['206703','206704','206700']
HMP2IBD = ['206701','206708','206709']
MHnormal = ['MH0005','MH0006','MH0008']
MHobese = ['MH0001','MH0002','MH0003']
grouparr = [diabetesctrl,diabetestype2nomet,HMP2Normal,HMP2IBD,MHnormal,MHobese]
grouplabelarr = ['diabetesctrl','diabetestype2nomet','HMP2Normal','HMP2IBD','MHnormal','MHobese']
for z1 in range(len(grouparr)):
    for z in range(len(grouparr[z1])):
        interactMat = np.zeros([773,773])
        bigmap = {}
        for afileorig in os.listdir(inputDir):
            if afileorig.endswith('.txt'):
                afile = afileorig[:len(afileorig)-4]
            words = afile.split('_')
            if words[1]==grouparr[z1][z]:
                ithIdx = int(words[2])
                jthIdx = int(words[3])
                inFI = open(inputDir+afileorig)
                inFI.readline()
                for line in inFI:
                    linewords = line.strip().split('\t')
                    speciesnum = -1
                    if linewords[2].split('_')[0]=='1' or linewords[2].split('_')[0]=='2':
                        speciesnum = int(linewords[2].split('_')[0])
                        rxn = linewords[2].split('_')[1]
                        #print(rxn)
                        if ithIdx not in bigmap:
                            bigmap[ithIdx] = {}
                        if jthIdx not in bigmap[ithIdx]:
                            bigmap[ithIdx][jthIdx] = {}
                        if rxn not in bigmap[ithIdx][jthIdx]:
                            bigmap[ithIdx][jthIdx][rxn] = [0,0]
                        #print('HERE')
                        bigmap[ithIdx][jthIdx][rxn][speciesnum-1] = float(linewords[0])

        #nonsense = nonsense+1
        for ithIdx in bigmap:
            for jthIdx in bigmap[ithIdx]:
                competecount = 0
                coopcount = 0
                for rxn in bigmap[ithIdx][jthIdx]:
                    if rxn in rxnstosubs and (rxnstosubs[rxn]=='Transport, extracellular' or rxnstosubs[rxn]=='Transport'):
                        flux1 = bigmap[ithIdx][jthIdx][rxn][0]
                        flux2 = bigmap[ithIdx][jthIdx][rxn][1]
                        if flux1!=0 and flux2!=0:
                            if (flux1>0 and flux2>0) or (flux1<0 and flux2<0):
                                competecount += 1
                                interactMat[ithIdx,jthIdx] += -abs(flux1)-abs(flux2)
                            elif (flux1>0 and flux2<0) or (flux1<0 and flux2>0):
                                coopcount += 1
                                interactMat[ithIdx,jthIdx] += abs(flux1)+abs(flux2)
                    if competecount>0 or coopcount>0:
                        pass
                        #print(ithIdx)
                        #print(jthIdx)
                        #print(competecount)
                        #print(coopcount)
        interactMat2 = copy.deepcopy(interactMat)
        for i in range(len(interactMat)):
            nonzerovals = interactMat[i,interactMat[i,:]!=0]
            if len(nonzerovals)!=0:
                nonzeromean = np.mean(nonzerovals)
                for j in range(len(interactMat[i])):
                    interactMat2[i,j] += nonzeromean
        for i in range(len(interactMat[0])):
            nonzerovals = interactMat[interactMat[:,i]!=0,i]
            if len(nonzerovals)!=0:
                nonzeromean = np.mean(nonzerovals)
                for j in range(len(interactMat)):
                    interactMat2[j,i] += nonzeromean
        [eigvals,eigvectors] = np.linalg.eig(interactMat2)
        print(grouplabelarr[z1])
        print(grouparr[z1][z])
        print(eigvals[:10])
        outFI = open('/mnt/vdb/home/ubuntu2/interactMatTemp'+grouparr[z1][z]+'.txt','w')
        for i in range(len(interactMat2)):
            for j in range(len(interactMat2)-1):
                outFI.write(str(interactMat2[i,j])+'\t')
            outFI.write(str(interactMat2[i,772])+'\n')
        outFI.close()
