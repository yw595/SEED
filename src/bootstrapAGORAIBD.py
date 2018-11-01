import numpy as np
import os
import math
from scipy.stats import hypergeom
import time
import scipy
import scipy.stats

def writeData(dataarrs,outfile,headers=[],delimiter='\t'):
    outFI = open(outfile,'w')
    if len(headers)!=0:
        for i in range(len(headers)-1):
            outFI.write(headers[i]+delimiter)
        outFI.write(headers[len(headers)-1]+'\n')
    for i in range(len(dataarrs[0])):
        for j in range(len(dataarrs)-1):
            dataitem = dataarrs[j][i]
            if not isinstance(dataitem,str):
                dataitem = str(dataitem)
            outFI.write(dataitem+delimiter)
        dataitem = dataarrs[len(dataarrs)-1][i]
        if not isinstance(dataitem,str):
            dataitem = str(dataitem)
        outFI.write(dataitem+'\n')
    outFI.close()

def addIdxStrings(array1):
    maxDigits = int(math.ceil(math.log(len(array1),10)))
    array2 = []
    for i in range(len(array1)):
        ithString = str(i)
        while len(ithString)<maxDigits:
            ithString = '0'+ithString
        array2.append(ithString+'_'+array1[i])
    return array2
    
def absArr(array1):
    array2 = []
    for i in range(len(array1)):
        array2.append(abs(array1[i]))
    return array2
    
def minArr(array1,array2):
    array3 = []
    for i in range(len(array1)):
        if array1[i] <= array2[i]:
            array3.append(array1[i])
        else:
            array3.append(array2[i])
    return array3

def sortByIdx(array,sortIdxs):
    array2 = []
    for i in range(len(sortIdxs)):
        array2.append(array[sortIdxs[i]])
    return array2

def varianceMetric(normalfluxes1,obesefluxes1):
    numer = np.var(normalfluxes1)+np.var(obesefluxes1)
    allfluxes = []
    for i in range(len(normalfluxes1)):
        allfluxes.append(normalfluxes1[i])
    for i in range(len(obesefluxes1)):
        allfluxes.append(obesefluxes1[i])
    denom = np.var(allfluxes)
    if denom==0:
        metric = -1
    else:
        metric = numer/denom
    return metric

def segmentFluxBySubsystem(modelrxns,modelsubs,fluxDist1,diffCond=False,innerAbs=True):
    uniqSubs = list(set(modelsubs))

    subLabels = []
    subFluxSums = []
    subFluxDiffSums = []
    subFluxStdSums = []
    subFluxDiffScaledSums = []
    for j in range(len(uniqSubs)):
        subLabels.append(uniqSubs[j])
        fluxsum = 0
        for i in range(len(fluxDist1)):
            if modelsubs[i]==uniqSubs[j]:
                if innerAbs:
                    fluxsum += abs(fluxDist1[i])
                else:
                    fluxsum += fluxDist1[i]
        subFluxSums.append(fluxsum)
    return [subLabels,subFluxSums]

def segmentFluxBySubsystemGroup(NormalIdxs,IBDIdxs,picrustFluxesArr,tobemergedrxns,tobemergedrxnnames,tobemergedsubsystems):
    subDiffMapNormal = {}
    subDiffMapIBD = {}
    subDiffMapIndividual = {}
    rxnDiffArr = []
    rxnDiffArrNum = []
    rxnDiffMapNormal = {}
    rxnDiffMapIBD = {}
    rxnDiffMapIndividual = {}
    for z in range(len(picrustFluxesArr)):
        for k in range(len(tobemergedrxns)):
            if z in NormalIdxs or z in IBDIdxs:
                #print('k'+tobemergedrxns[k])
                if tobemergedrxns[k] not in rxnDiffMapIndividual:
                    rxnDiffMapIndividual[tobemergedrxns[k]] = [0 for k1 in range(len(picrustFluxesArr))]
                rxnDiffMapIndividual[tobemergedrxns[k]][z] = picrustFluxesArr[z][k]
            if z in IBDIdxs:
                if tobemergedrxns[k] not in rxnDiffMapIBD:
                    #print(rxnDiffMapIBD)
                    #print(tobemergedrxns[k])
                    rxnDiffMapIBD[tobemergedrxns[k]] = picrustFluxesArr[z][k]
                else:
                    #print(rxnDiffMapIBD)
                    rxnDiffMapIBD[tobemergedrxns[k]] += picrustFluxesArr[z][k]
            elif z in NormalIdxs:
                if tobemergedrxns[k] not in rxnDiffMapNormal:
                    #print(rxnDiffMapNormal)
                    #print(tobemergedrxns[k])
                    rxnDiffMapNormal[tobemergedrxns[k]] = picrustFluxesArr[z][k]
                else:
                    rxnDiffMapNormal[tobemergedrxns[k]] += picrustFluxesArr[z][k]

        [subLabelsArr,subFluxSums] = segmentFluxBySubsystem(tobemergedrxns,tobemergedsubsystems,picrustFluxesArr[z],diffCond=False,innerAbs=True)
        for k in range(len(subLabelsArr)):
	    if subLabelsArr[k] not in subDiffMapIndividual:
	        subDiffMapIndividual[subLabelsArr[k]] = [0 for k1 in range(len(picrustFluxesArr))]
	    subDiffMapIndividual[subLabelsArr[k]][z] = subFluxSums[k]
	    if z in IBDIdxs:
	        if subLabelsArr[k] not in subDiffMapIBD:
		    subDiffMapIBD[subLabelsArr[k]] = subFluxSums[k]
	        else:
		    subDiffMapIBD[subLabelsArr[k]] += subFluxSums[k]
	    elif z in NormalIdxs:
	        if subLabelsArr[k] not in subDiffMapNormal:
		    subDiffMapNormal[subLabelsArr[k]] = subFluxSums[k]
	        else:
		    subDiffMapNormal[subLabelsArr[k]] += subFluxSums[k]

    #nonsense += 1
    rxnKeys = rxnDiffMapNormal.keys()
    rxnDiffArr = []
    rxnDiffArrNum = []
    for i in range(len(rxnKeys)):
        rxnDiffArr.append(rxnKeys[i])
        rxnDiffArrNum.append(rxnDiffMapNormal[rxnKeys[i]]/len(NormalIdxs) - rxnDiffMapIBD[rxnKeys[i]]/len(IBDIdxs))
    subKeys = subDiffMapNormal.keys()
    subDiffArr = []
    subDiffArrNum = []
    subDiffWilcPVals = []
    for i in range(len(subKeys)):
        subDiffArr.append(subKeys[i])
        subDiffArrNum.append(subDiffMapNormal[subKeys[i]]/len(NormalIdxs) - subDiffMapIBD[subKeys[i]]/len(IBDIdxs))

    return [subDiffArr,subDiffArrNum,subDiffMapIndividual,rxnDiffArr,rxnDiffArrNum,rxnDiffMapIndividual]

if __name__=='__main__':
    useDiabetes = True
    useNormObese = False
    useWithMet = False
    if useDiabetes:
        ctrlFile = '/mnt/vdb/home/ubuntu2/diabetesctrl.txt'
        if useWithMet:
            type2File = '/mnt/vdb/home/ubuntu2/diabetestype2withmet.txt'
        else:
            type2File = '/mnt/vdb/home/ubuntu2/diabetestype2nomet.txt'
        inputsuffix = 'diabetesnormmerged'
        if useWithMet:
            outputsuffix = 'DiabetesWithMetNorm'
        else:
            outputsuffix = 'DiabetesNorm'
        inputdir = '/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData/Diabetes/'
    elif useNormObese:
        ctrlFile = '/mnt/vdb/home/ubuntu2/MHnormal.txt'
        type2File = '/mnt/vdb/home/ubuntu2/MHobese.txt'
        inputsuffix = 'normobesenormmerged'
        outputsuffix = 'NormObeseNorm'
        inputdir = '/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData/'
    ctrlFI = open(ctrlFile)
    ctrlsamps = []
    for line in ctrlFI:
        ctrlsamps.append(line.strip())
    ctrlFI.close()
    type2FI = open(type2File)
    type2samps = []
    for line in type2FI:
        type2samps.append(line.strip())
    type2FI.close()

    allfiles = os.listdir('/mnt/vdb/home/ubuntu2')
    picrustFluxesArr = []
    acArr = []
    NormalIdxs = []
    IBDIdxs = []
    for i in range(len(allfiles)):
        if allfiles[i].find(inputsuffix)!=-1:
            sampname = allfiles[i][0:allfiles[i].find(inputsuffix)]
            inFI = open('/mnt/vdb/home/ubuntu2/'+allfiles[i])
            fluxes = []
            for line in inFI:
                fluxes.append(float(line.strip()))
            picrustFluxesArr.append(fluxes)
            inFI.close()
            exps = []
            inFI = open(inputdir+sampname+'/normalized_otus.modelexpr')
            for line in inFI:
                exps.append(float(line.strip()))
            acArr.append(exps)
            inFI.close()
            if useDiabetes:
                if sampname in ctrlsamps:
                    NormalIdxs.append(len(picrustFluxesArr)-1)
                if sampname in type2samps:
                    IBDIdxs.append(len(picrustFluxesArr)-1)
            elif useNormObese:
                if sampname[6:] in ctrlsamps:
                    NormalIdxs.append(len(picrustFluxesArr)-1)
                if sampname[6:] in type2samps:
                    IBDIdxs.append(len(picrustFluxesArr)-1)

    corr1 = []
    corr2 = []
    for i in range(len(picrustFluxesArr)):
        for j in range(len(picrustFluxesArr[i])):
            corr1.append(picrustFluxesArr[i][j])
            corr2.append(acArr[i][j])
    corr1 = abs(np.array(corr1))
    corr2 = abs(np.array(corr2))
    corr2 = corr2/sum(corr2)
    [rho, pval] = scipy.stats.spearmanr(corr1,corr2)
    corr1temp = corr1[np.logical_and(corr1!=0,corr2!=0)]
    corr2temp = corr2[np.logical_and(corr1!=0,corr2!=0)]
    [rhononzero, pvalnonzero] = scipy.stats.spearmanr(corr1temp,corr2temp)
    print(rho)
    print(pval)
    print(rhononzero)
    print(pvalnonzero)

    writeData([corr1,corr2],'/mnt/vdb/home/ubuntu2/exprVsFlux'+outputsuffix+'.txt',headers=['flux','expr'])
    writeData([corr1temp,corr2temp],'/mnt/vdb/home/ubuntu2/exprVsFluxNonzero'+outputsuffix+'.txt',headers=['flux','expr'])
    
    nonsense = nonsense+1
                    
    tobemergedrxns = []
    tobemergedrxnnames = []
    tobemergedsubsystems = []
    tobemergedFI = open('/mnt/vdb/home/ubuntu2/tobemerged2.txt')
    tobemergedFI.readline()
    for line in tobemergedFI:
        words = line.strip().split('\t')
        tobemergedrxns.append(words[0])
        tobemergedrxnnames.append(words[1])
        if len(words)>2:
            tobemergedsubsystems.append(words[2])
        else:
            tobemergedsubsystems.append('')
    tobemergedFI.close()

    [expSubDiffArr,expSubDiffArrNum,expSubDiffMapIndividual,expRxnDiffArr,expRxnDiffArrNum,expRxnDiffMapIndividual] = segmentFluxBySubsystemGroup(NormalIdxs,IBDIdxs,acArr,tobemergedrxns,tobemergedrxnnames,tobemergedsubsystems)
    [subDiffArr,subDiffArrNum,subDiffMapIndividual,rxnDiffArr,rxnDiffArrNum,rxnDiffMapIndividual] = segmentFluxBySubsystemGroup(NormalIdxs,IBDIdxs,picrustFluxesArr,tobemergedrxns,tobemergedrxnnames,tobemergedsubsystems)
    
    randIters = 200
    expRxnDiffKthMap = {}
    rxnMetricKthMap = {}
    rxnDiffKthMap = {}
    subDiffKthMap = {}
    for k in range(randIters):
        currentTime1 = time.time()
        print(k)
        NormalAndIBDIdxs = []
        for i in range(len(NormalIdxs)):
            NormalAndIBDIdxs.append(NormalIdxs[i])
        for i in range(len(IBDIdxs)):
            NormalAndIBDIdxs.append(IBDIdxs[i])
        randidxs1 = np.random.choice(len(NormalIdxs)+len(IBDIdxs),len(NormalIdxs),replace=False)
        randidxs1temp = []
        for i in range(len(randidxs1)):
            randidxs1temp.append(NormalAndIBDIdxs[randidxs1[i]])
        randidxs1 = randidxs1temp
        randidxs2 = []
        for k1 in range(len(picrustFluxesArr)):
            if k1 not in randidxs1 and (k1 in NormalIdxs or k1 in IBDIdxs):
                randidxs2.append(k1)
        #print(randidxs1)
        #print(randidxs2)
        #print(NormalAndIBDIdxs)
        #print(NormalIdxs)
        #print(IBDIdxs)
        #print(randidxs1)
        #print(randidxs2)
        [expSubDiffArrKth,expSubDiffArrNumKth,expSubDiffMapIndividualKth,expRxnDiffArrKth,expRxnDiffArrNumKth,expRxnDiffMapIndividualKth] = segmentFluxBySubsystemGroup(randidxs1,randidxs2,acArr,tobemergedrxns,tobemergedrxnnames,tobemergedsubsystems)
        for k1 in range(len(expRxnDiffArrKth)):
            if expRxnDiffArrKth[k1] not in expRxnDiffKthMap:
                expRxnDiffKthMap[expRxnDiffArrKth[k1]] = []
                expRxnDiffKthMap[expRxnDiffArrKth[k1]].append(expRxnDiffArrNumKth[k1])
            else:
                expRxnDiffKthMap[expRxnDiffArrKth[k1]].append(expRxnDiffArrNumKth[k1])
                
        [subDiffArrKth,subDiffArrNumKth,subDiffMapIndividualKth,rxnDiffArrKth,rxnDiffArrNumKth,rxnDiffMapIndividualKth] = segmentFluxBySubsystemGroup(randidxs1,randidxs2,picrustFluxesArr,tobemergedrxns,tobemergedrxnnames,tobemergedsubsystems)
        for k1 in range(len(rxnDiffArrKth)):
            k1thfluxes = rxnDiffMapIndividualKth[rxnDiffArrKth[k1]]
            normalfluxes1 = []
            for i in range(len(randidxs1)):
                normalfluxes1.append(k1thfluxes[randidxs1[i]])
            obesefluxes1 = []
            for i in range(len(randidxs2)):
                obesefluxes1.append(k1thfluxes[randidxs2[i]])
            metric = varianceMetric(normalfluxes1,obesefluxes1)
            if rxnDiffArrKth[k1] not in rxnMetricKthMap:
                rxnMetricKthMap[rxnDiffArrKth[k1]] = []
                rxnMetricKthMap[rxnDiffArrKth[k1]].append(metric)
            else:
                rxnMetricKthMap[rxnDiffArrKth[k1]].append(metric)
        for k1 in range(len(rxnDiffArrKth)):
            if rxnDiffArrKth[k1] not in rxnDiffKthMap:
                rxnDiffKthMap[rxnDiffArrKth[k1]] = []
                rxnDiffKthMap[rxnDiffArrKth[k1]].append(rxnDiffArrNumKth[k1])
            else:
                rxnDiffKthMap[rxnDiffArrKth[k1]].append(rxnDiffArrNumKth[k1])
        for k1 in range(len(subDiffArrKth)):
            if subDiffArrKth[k1] not in subDiffKthMap:
                subDiffKthMap[subDiffArrKth[k1]] = []
                subDiffKthMap[subDiffArrKth[k1]].append(subDiffArrNumKth[k1])
            else:
                subDiffKthMap[subDiffArrKth[k1]].append(subDiffArrNumKth[k1])
        currentTime2 = time.time()
        print(currentTime2-currentTime1)

    for k1 in range(len(rxnDiffArr)):
        rxnDiffKthMap[rxnDiffArr[k1]].sort(reverse=True)
    for k1 in range(len(rxnDiffArr)):
        expRxnDiffKthMap[rxnDiffArr[k1]].sort(reverse=True)
    for k1 in range(len(rxnDiffArr)):
        rxnMetricKthMap[rxnDiffArr[k1]].sort()
    for k1 in range(len(subDiffArr)):
        subDiffKthMap[subDiffArr[k1]].sort(reverse=True)
                
    topSubs = []
    bottomSubs = []
    topPsSub = []
    bottomPsSub = []
    topRxns = []
    bottomRxns = []
    topPsRxn = []
    bottomPsRxn = []
    topRxnsT = []
    topPsRxnT = []
    topRxnsMetric = []
    topNumsRxnMetric = []
    topPsRxnMetric = []
    topPsRxnMetricNames = []
    topExpRxns = []
    bottomExpRxns = []
    topPsExpRxn = []
    bottomPsExpRxn = []
    topExpRxnsT = []
    topPsExpRxnT = []
    topExpRxnsMetric = []
    topNumsExpRxnMetric = []
    for k1 in range(len(expRxnDiffArr)):
        individualMap = expRxnDiffMapIndividual[expRxnDiffArr[k1]]
        [stat, p] = scipy.stats.ranksums(np.array(individualMap)[IBDIdxs],np.array(individualMap)[NormalIdxs])
        topExpRxnsT.append(expRxnDiffArr[k1])
        topPsExpRxnT.append(p)
        expRxnDiffArrNumKth = expRxnDiffArrNum[k1]
        expRxnDiffKthMapVals = expRxnDiffKthMap[expRxnDiffArr[k1]]
        top = -1
        for k2 in range(len(expRxnDiffKthMapVals)):
            if expRxnDiffKthMapVals[k2] >= expRxnDiffArrNumKth:
                top=k2
        bottom = len(expRxnDiffKthMapVals)
        for k2 in range(len(expRxnDiffKthMapVals)-1,-1,-1):
            if expRxnDiffKthMapVals[k2] <= expRxnDiffArrNumKth:
                bottom=k2
        if True:
            topExpRxns.append(expRxnDiffArr[k1])
            topPsExpRxn.append(float(top+1)/randIters)
            bottomExpRxns.append(expRxnDiffArr[k1])
            bottomPsExpRxn.append(float(randIters-(bottom-1))/randIters)
    for k1 in range(len(rxnDiffArr)):
        individualMap = rxnDiffMapIndividual[rxnDiffArr[k1]]
        [stat, p] = scipy.stats.ranksums(np.array(individualMap)[IBDIdxs],np.array(individualMap)[NormalIdxs])
        topRxnsT.append(rxnDiffArr[k1])
        topPsRxnT.append(p)
        rxnDiffArrNumKth = rxnDiffArrNum[k1]
        rxnDiffKthMapVals = rxnDiffKthMap[rxnDiffArr[k1]]
        top = -1
        for k2 in range(len(rxnDiffKthMapVals)):
            if rxnDiffKthMapVals[k2] >= rxnDiffArrNumKth:
                top=k2
        bottom = len(rxnDiffKthMapVals)
        for k2 in range(len(rxnDiffKthMapVals)-1,-1,-1):
            if rxnDiffKthMapVals[k2] <= rxnDiffArrNumKth:
                bottom=k2
        if True:
            topRxns.append(rxnDiffArr[k1])
            topPsRxn.append(float(top+1)/randIters)
            bottomRxns.append(rxnDiffArr[k1])
            bottomPsRxn.append(float(randIters-(bottom-1))/randIters)
    for k1 in range(len(subDiffArr)):
        individualMap = subDiffMapIndividual[subDiffArr[k1]]
        subDiffArrNumKth = subDiffArrNum[k1]
        subDiffKthMapVals = subDiffKthMap[subDiffArr[k1]]
        top = -1
        for k2 in range(len(subDiffKthMapVals)):
            if subDiffKthMapVals[k2] >= subDiffArrNumKth:
                top=k2
        bottom = len(subDiffKthMapVals)
        for k2 in range(len(subDiffKthMapVals)-1,-1,-1):
            if subDiffKthMapVals[k2] <= subDiffArrNumKth:
                bottom=k2
        if True:
            topSubs.append(subDiffArr[k1])
            topPsSub.append(float(top+1)/randIters)
            bottomSubs.append(subDiffArr[k1])
            bottomPsSub.append(float(randIters-(bottom-1))/randIters)

    sortIdxs = zip(absArr(rxnDiffArrNum),range(len(rxnDiffArrNum)))
    sortIdxs.sort(reverse=True)
    sortIdxsTemp = []
    for i in range(len(sortIdxs)):
        sortIdxsTemp.append(sortIdxs[i][1])
    sortIdxs = sortIdxsTemp
    rxnDiffArr = sortByIdx(rxnDiffArr,sortIdxs)
    rxnDiffArrNum = sortByIdx(rxnDiffArrNum,sortIdxs)
    rxnSubArr = []
    rxnTsArr = []
    rxnNameDiffArr = []
    rxnPermPValsArr = []
    rxnExpDiffArr = []
    rxnExpPValsArr = []
    for i in range(len(rxnDiffArr)):
        rxnSubArr.append(tobemergedsubsystems[tobemergedrxns.index(rxnDiffArr[i])])
        rxnNameDiffArr.append(tobemergedrxnnames[tobemergedrxns.index(rxnDiffArr[i])])
        rxnTsArr.append(topPsRxnT[topRxnsT.index(rxnDiffArr[i])])
        rxnPermPValsArr.append(min(topPsRxn[topRxns.index(rxnDiffArr[i])],bottomPsRxn[bottomRxns.index(rxnDiffArr[i])]))
        rxnExpDiffArr.append(expRxnDiffArrNum[expRxnDiffArr.index(rxnDiffArr[i])])
        rxnExpPValsArr.append(min(topPsExpRxn[topExpRxns.index(rxnDiffArr[i])],bottomPsExpRxn[bottomExpRxns.index(rxnDiffArr[i])]))
    writeData([addIdxStrings(rxnDiffArr),rxnDiffArrNum,rxnSubArr,rxnNameDiffArr,rxnPermPValsArr,rxnExpDiffArr,rxnExpPValsArr,rxnTsArr],'/mnt/vdb/home/ubuntu2/mergedModel'+outputsuffix+'RxnDiffs.txt',headers=['rxn','diffflux','subsystem','rxnname','bootstrap p-val','expression difference','expression bootstrap p-val','wilcoson rank-sum p-value'])

    sortIdxs = zip(absArr(subDiffArrNum),range(len(subDiffArrNum)))
    sortIdxs.sort(reverse=True)
    sortIdxsTemp = []
    for i in range(len(sortIdxs)):
        sortIdxsTemp.append(sortIdxs[i][1])
    sortIdxs = sortIdxsTemp
    subDiffArr = sortByIdx(subDiffArr,sortIdxs)
    subDiffArrNum = sortByIdx(subDiffArrNum,sortIdxs)
    rxnNumSubMap = {}
    rxnNumSigSubMap = {}
    for i in range(len(tobemergedsubsystems)):
        if tobemergedsubsystems[i] not in rxnNumSubMap:
            rxnNumSubMap[tobemergedsubsystems[i]] = 1
        else:
            rxnNumSubMap[tobemergedsubsystems[i]] += 1
        if rxnPermPValsArr[rxnDiffArr.index(tobemergedrxns[i])] < .05:
            if tobemergedsubsystems[i] not in rxnNumSigSubMap:
                rxnNumSigSubMap[tobemergedsubsystems[i]] = 1
            else:
                rxnNumSigSubMap[tobemergedsubsystems[i]] += 1
    rxnNumSubArr = []
    rxnNumSigSubArr = []
    hypergeomparr = []
    for sub in subDiffArr:
        rxnNumSubArr.append(rxnNumSubMap[sub])
        if sub in rxnNumSigSubMap:
            rxnNumSigSubArr.append(rxnNumSigSubMap[sub])
        else:
            rxnNumSigSubArr.append(0)
        M = rxnNumSubMap[sub]
        N = len(rxnDiffArr)
        K = sum([1 for i in range(len(rxnPermPValsArr)) if rxnPermPValsArr[i] < .05])
        if sub in rxnNumSigSubMap:
            x = rxnNumSigSubMap[sub]
        else:
            x = 0
        hypergeomparr.append(1-hypergeom.cdf(x-1,N,M,K))
    writeData([addIdxStrings(subDiffArr),subDiffArrNum,rxnNumSubArr,rxnNumSigSubArr,hypergeomparr],'/mnt/vdb/home/ubuntu2/mergedModel'+outputsuffix+'SubDiffs.txt',headers=['sub','diffflux','total number of rxns','number of significant rxns','Hypergeometric p-val'])

    for i in range(len(rxnDiffArr)):
        if rxnPermPValsArr[i] < .05:
            strrepString = rxnDiffArr[i].replace('/','_').replace(' ','_')
            rxnDiffKthMapVals = rxnDiffKthMap[rxnDiffArr[i]]
            writeData([rxnDiffKthMapVals],'/mnt/vdb/home/ubuntu2/mergedModelRxnDists'+outputsuffix+'/'+strrepString+'KthMap.txt',headers=['fluxdiff'])
            rxnDiffMapIndividualVals = rxnDiffMapIndividual[rxnDiffArr[i]]
            highlightArr = []
            rxnDiffMapIndividualValsTemp = []
            for k1 in range(len(rxnDiffMapIndividualVals)):
                if k1 in NormalIdxs:
                    highlightArr.append('normal')
                    rxnDiffMapIndividualValsTemp.append(rxnDiffMapIndividualVals[k1])
                elif k1 in IBDIdxs:
                    highlightArr.append('diabetes')
                    rxnDiffMapIndividualValsTemp.append(rxnDiffMapIndividualVals[k1])
            rxnDiffMapIndividualVals = rxnDiffMapIndividualValsTemp
            writeData([rxnDiffMapIndividualVals,highlightArr],'/mnt/vdb/home/ubuntu2/mergedModelRxnDists'+outputsuffix+'/'+strrepString+'DiffMap.txt',headers=['flux','highlight'])
