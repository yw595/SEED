import os
import numpy as np
import copy
import scipy.stats
from scipy.stats import hypergeom
from bootstrapAGORAIBD import writeData

def calcEigs(bigmap,eigminarr,eigmaxarr,eigmeanarr):
    mergedFI = open('/mnt/vdb/home/ubuntu2/tobemerged2.txt')
    mergedFI.readline()
    rxnstosubs = {}
    for line in mergedFI:
        words = line.strip().split('\t')
        if len(words)==3:
            rxnstosubs[words[0]] = words[2]
    mergedFI.close()
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
    if len(eigminarr)==0:
        for k in range(len(eigvals)):
            eigminarr.append(eigvals[k].real)
    else:
        for k in range(len(eigminarr)):
            if eigvals[k]<eigminarr[k]:
                eigminarr[k] = eigvals[k].real
    if len(eigmaxarr)==0:
        for k in range(len(eigvals)):
            eigmaxarr.append(eigvals[k].real)
    else:
        for k in range(len(eigmaxarr)):
            if eigvals[k]>eigmaxarr[k]:
                eigmaxarr[k] = eigvals[k].real
    if len(eigmeanarr)==0:
        for k in range(len(eigvals)):
            eigmeanarr.append(eigvals[k].real)
    else:
        for k in range(len(eigmeanarr)):
            eigmeanarr[k] += eigvals[k].real
    print(grouplabelarr[z1])
    print(grouparr[z1][z])
    print(eigvals[:10])
    outFI = open('/mnt/vdb/home/ubuntu2/interactMatTemp'+grouparr[z1][z]+'.txt','w')
    for i in range(len(interactMat2)):
        for j in range(len(interactMat2)-1):
            outFI.write(str(interactMat2[i,j])+'\t')
        outFI.write(str(interactMat2[i,772])+'\n')
    outFI.close()
    return [eigminarr,eigmaxarr,eigmeanarr]
    
if __name__=='__main__':

    inFI = open('/mnt/vdb/home/ubuntu2/tobemerged2.txt')
    tobemerged2map = {}
    for line in inFI:
        words = line.strip().split('\t')
        if len(words)==2:
            words.append('')
        tobemerged2map[words[0]] = [words[1],words[2]]
    inFI.close()

    inputDir = '/mnt/vdb/home/ubuntu2/MATLAB/SEED/output/simulateSmallModelsCombined/'

    #diabetesctrl = ['BGI-06A','N044A','SZEY-75A']
    #diabetestype2nomet = ['NG-5636_551','DOM024','DLM001']
    diabetesctrl = ['BGI-06A','N044A','SZEY-75A','MH0189','MH0181','N029A','NLM031','MH0439','MH0431','MH0024']
    diabetestype2nomet = ['NG-5636_551','DOM024','DLM001','MH0334','DLM027','MH0370','DOM023','DLM019','DLM028','DLM012','MH0345']
    HMP2Normal = ['206703','206704','206700']
    HMP2IBD = ['206701','206708','206709']
    MHnormal = ['MH0005','MH0006','MH0008']
    MHobese = ['MH0001','MH0002','MH0003']
    grouparr = [diabetesctrl,diabetestype2nomet]#,HMP2Normal,HMP2IBD,MHnormal,MHobese]
    grouplabelarr = ['diabetesctrl','diabetestype2nomet']#,'HMP2Normal','HMP2IBD','MHnormal','MHobese']
    analyzeSingle = True
    for z1 in range(len(grouparr)):
        if not analyzeSingle:
            eigminarr = []
            eigmaxarr = []
            eigmeanarr = []
        if z1%2==0:
            rxndiffmap = {}
        for z in range(len(grouparr[z1])):
            count = 0
            interactMat = np.zeros([773,773])
            bigmap = {}
            for afileorig in os.listdir(inputDir):
                if afileorig.endswith('.txt'):
                    afile = afileorig[:len(afileorig)-4]
                words = afile.split('_')
                if words[1]==grouparr[z1][z] and ((analyzeSingle and len(words)==3) or (not analyzeSingle and len(words)==4)):
                    count = count+1
                    print(count)
                    ithIdx = int(words[2])
                    if not analyzeSingle:
                        jthIdx = int(words[3])
                    inFI = open(inputDir+afileorig)
                    inFI.readline()
                    for line in inFI:
                        linewords = line.strip().split('\t')
                        speciesnum = -1
                        if (linewords[2].split('_')[0]=='1' or linewords[2].split('_')[0]=='2') or analyzeSingle:
                            if (linewords[2].split('_')[0]=='1' or linewords[2].split('_')[0]=='2'):
                                speciesnum = int(linewords[2].split('_')[0])
                                rxn = linewords[2].split('_')
                                rxn = '_'.join(rxn[1:])
                            elif analyzeSingle:
                                rxn = linewords[2]
                            rxn = rxn.replace('_LSQBKT','')
                            rxn = rxn.replace('_RSQBKT','')
                            rxn = rxn.replace('_LPAREN','')
                            rxn = rxn.replace('_RPAREN','')
                            rxn = rxn.replace('DASH','-')
                            if rxn not in rxndiffmap:
                                rxndiffmap[rxn] = [{},{}]
                            if z1%2==0:
                                if ithIdx not in rxndiffmap[rxn][0]:
                                    rxndiffmap[rxn][0][ithIdx] = []
                                rxndiffmap[rxn][0][ithIdx].append(float(linewords[0]))
                            else:
                                if ithIdx not in rxndiffmap[rxn][1]:
                                    rxndiffmap[rxn][1][ithIdx] = []
                                rxndiffmap[rxn][1][ithIdx].append(float(linewords[0]))
                            if not analyzeSingle:
                                if ithIdx not in bigmap:
                                    bigmap[ithIdx] = {}
                                if jthIdx not in bigmap[ithIdx]:
                                    bigmap[ithIdx][jthIdx] = {}
                                if rxn not in bigmap[ithIdx][jthIdx]:
                                    bigmap[ithIdx][jthIdx][rxn] = [0,0]
                                #print('HERE')
                                bigmap[ithIdx][jthIdx][rxn][speciesnum-1] = float(linewords[0])

            if not analyzeSingle:
                [eigminarr,eigmaxarr,eigmeanarr] = calcEigs(bigmap,eigminarr,eigmaxarr,eigmeanarr)

        if not analyzeSingle:
            onetoten = ['_01','_02','_03','_04','_05','_06','_07','_08','_09','_10']
            for k in range(len(eigmeanarr)):
                eigmeanarr[k] = eigmeanarr[k]/3
            writeData([onetoten,eigmeanarr,eigminarr,eigmaxarr],'/mnt/vdb/home/ubuntu2/eig'+grouplabelarr[z1]+'.txt',delimiter='\t',headers=['eignum','eigval','lower','upper'])

        if z1%2==1:
            if analyzeSingle:
                rxndiffmaptemp = {}
                for rxn in rxndiffmap:
                    rxndiffmaptemp[rxn] = [[],[]]
                    for k in [0,1]:
                        for ithIdx in rxndiffmap[rxn][k]:
                            rxndiffmaptemp[rxn][k].append(np.mean(rxndiffmap[rxn][k][ithIdx]))
                            if k==0:
                                k1=1
                            else:
                                k1=0
                            if ithIdx in rxndiffmap[rxn][k1]:
                                rxndiffmaptemp[rxn][k1].append(np.mean(rxndiffmap[rxn][k1][ithIdx]))
                            else:
                                rxndiffmaptemp[rxn][k1].append(0)
                rxndiffmap = rxndiffmaptemp
            else:
                rxndiffmaptemp = {}
                for rxn in rxndiffmap:
                    rxndiffmaptemp[rxn] = [[],[]]
                    for k in [0,1]:
                        for ithIdx in rxndiffmap[rxn][k]:
                            for l in range(len(rxndiffmap[rxn][k][ithIdx])):
                                rxndiffmaptemp[rxn][k].append(rxndiffmap[rxn][k][ithIdx][l])
                rxndiffmap = rxndiffmaptemp
                
        if z1==1:
            diabetesrxndiffmap = rxndiffmap
        if z1==3:
            IBDrxndiffmap = rxndiffmap
        if z1==5:
            normobeserxndiffmap = rxndiffmap
        if z1%2==1:
            rxnarr = rxndiffmap.keys()
            rxnarr.sort()
            pvalarr = []
            rxnnamearr = []
            subarr = []
            meandiffarr = []
            for rxn in rxnarr:
                maxlen = max(len(rxndiffmap[rxn][0]),len(rxndiffmap[rxn][1]))
                if len(rxndiffmap[rxn][0])==0:
                    for k in range(len(rxndiffmap[rxn][0]),maxlen):
                        rxndiffmap[rxn][0].append(0)
                if len(rxndiffmap[rxn][1])==0:
                    for k in range(len(rxndiffmap[rxn][1]),maxlen):
                        rxndiffmap[rxn][1].append(0)  
                if max(rxndiffmap[rxn][0])!=max(rxndiffmap[rxn][1]) or min(rxndiffmap[rxn][0])!=min(rxndiffmap[rxn][1]):
                    if not analyzeSingle:
                        [stat, pval] = scipy.stats.mannwhitneyu(rxndiffmap[rxn][0],rxndiffmap[rxn][1])
                    else:
                        statdiffarr = []
                        for k in range(len(rxndiffmap[rxn][0])):
                            statdiffarr.append(rxndiffmap[rxn][0][k]-rxndiffmap[rxn][1][k])
                        [stat, pval] = scipy.stats.ttest_1samp(statdiffarr,0)
                else:
                    pval = 1
                pvalarr.append(pval)
                rxncand = rxn
                if rxncand not in tobemerged2map:
                    rxncand = rxn.replace('_c_','c')
                if rxncand not in tobemerged2map:
                    rxncand = rxn.replace('_m_','m')
                if rxncand not in tobemerged2map:
                    rxncand = rxn.replace('_e_','e')
                if rxncand not in tobemerged2map:
                    rxncand = rxn.replace('_p_','p')
                if rxncand not in tobemerged2map:
                    rxncand = rxn.replace('EX','EX_')
                if rxncand not in tobemerged2map and rxn.startswith('EX_'):
                    rxncand = rxn.replace('_','')
                    rxncand = rxncand.replace('EX','EX_')
                    
                if rxncand not in tobemerged2map:
                    rxncand = rxn.replace('_','')
                rxnnamearr.append(tobemerged2map[rxncand][0])
                subarr.append(tobemerged2map[rxncand][1])
                meandiffarr.append(np.mean(rxndiffmap[rxn][0])-np.mean(rxndiffmap[rxn][1]))
                diffout = '/mnt/vdb/home/ubuntu2/pairwiseSingleSpeciesRxnDistsDiabetes/'+rxncand+'DiffMap.txt'
                rxndiffarr = []
                highlightarr = []
                for k in range(len(rxndiffmap[rxn][0])):
                    rxndiffarr.append(rxndiffmap[rxn][0][k])
                    highlightarr.append('ctrl')
                for k in range(len(rxndiffmap[rxn][1])):
                    rxndiffarr.append(rxndiffmap[rxn][1][k])
                    highlightarr.append('t2d')
                writeData([rxndiffarr,highlightarr],diffout,delimiter='\t',headers=['flux','highlight'])

            uniqSubs = np.unique(subarr)
            numrxnsall = []
            numsigrxnsall = []
            hypergeomparrall = []
            for sub in uniqSubs:
                N = len(rxnarr)
                M = sum(np.array(subarr)==sub)
                K = sum(np.array(pvalarr)<.05)
                x = sum(np.logical_and(np.array(pvalarr)<.05,np.array(subarr)==sub))
                numrxnsall.append(M)
                numsigrxnsall.append(x)
                hypergeomparrall.append(1-hypergeom.cdf(x-1,N,M,K))
            numrxnspos = []
            numsigrxnspos = []
            hypergeomparrpos = []
            for sub in uniqSubs:
                N = len(rxnarr)
                M = sum(np.array(subarr)==sub)
                K = sum(np.logical_and(np.array(pvalarr)<.05,np.array(meandiffarr)>0))
                x = sum(np.logical_and(np.logical_and(np.array(pvalarr)<.05,np.array(meandiffarr)>0),np.array(subarr)==sub))
                numrxnspos.append(M)
                numsigrxnspos.append(x)
                hypergeomparrpos.append(1-hypergeom.cdf(x-1,N,M,K))
            numrxnsneg = []
            numsigrxnsneg = []
            hypergeomparrneg = []
            for sub in uniqSubs:
                N = len(rxnarr)
                M = sum(np.array(subarr)==sub)
                K = sum(np.logical_and(np.array(pvalarr)<.05,np.array(meandiffarr)<0))
                x = sum(np.logical_and(np.logical_and(np.array(pvalarr)<.05,np.array(meandiffarr)<0),np.array(subarr)==sub))
                numrxnsneg.append(M)
                numsigrxnsneg.append(x)
                hypergeomparrneg.append(1-hypergeom.cdf(x-1,N,M,K))

            if z1==1:
                disease='diabetes'
            if z1==3:
                disease='IBD'
            if z1==5:
                disease='normobese'
            if analyzeSingle:
                disease = 'SingleSpecies'+disease
            writeData([rxnarr,rxnnamearr,subarr,pvalarr,meandiffarr],'/mnt/vdb/home/ubuntu2/pairwiseRxnDiff'+disease+'.txt',delimiter='\t',headers=['rxn','rxnname','sub','wilcoxon p-val','mean flux diff'])
            writeData([uniqSubs,numrxnsall,numsigrxnsall,hypergeomparrall],'/mnt/vdb/home/ubuntu2/pairwiseSubDiffAll'+disease+'.txt',delimiter='\t',headers=['sub','num reactions','num sig reactions','hypergeometric p-val'])
            writeData([uniqSubs,numrxnspos,numsigrxnspos,hypergeomparrpos],'/mnt/vdb/home/ubuntu2/pairwiseSubDiffPos'+disease+'.txt',delimiter='\t',headers=['sub','num reactions','num sig reactions','hypergeometric p-val'])
            writeData([uniqSubs,numrxnsneg,numsigrxnsneg,hypergeomparrneg],'/mnt/vdb/home/ubuntu2/pairwiseSubDiffNeg'+disease+'.txt',delimiter='\t',headers=['sub','num reactions','num sig reactions','hypergeometric p-val'])
