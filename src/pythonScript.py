import os
from scipy import stats
import math
import numpy as np
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

formatespeciesmap = {}
formatematrix = []
if True:
    files = os.listdir('/mnt/vdb/home/ubuntu2/MATLAB/SEED/output/simulateSmallModelsSeparatePicrustAGORA')

    obeselist = []
    normallist = []
    alllist = []
    inFI = open('/mnt/vdb/home/ubuntu2/MHobese.txt')
    for line in inFI:
        obeselist.append(line.strip())
        alllist.append(line.strip())
    inFI.close()
    inFI = open('/mnt/vdb/home/ubuntu2/MHnormal.txt')
    for line in inFI:
        normallist.append(line.strip())
        alllist.append(line.strip())
    inFI.close()

    for i in range(len(alllist)):
        formatematrix.append([])

    wordsMap = {}
    massiveArr = []
    for i in range(len(alllist)):
        massiveArr1 = []
        for j in range(773):
            massiveArr1.append([])
        massiveArr.append(massiveArr1)
    for i in range(len(files)):
        if files[i].endswith('.flux'):
            words = files[i].split('_')
            firstSpecies = int(words[0])
            samplename = words[3]
            samplename = samplename[samplename.find('ucrC97')+6:samplename.find('.flux')]
            if words[2]=='1':
                inFI = open('/mnt/vdb/home/ubuntu2/MATLAB/SEED/output/simulateSmallModelsSeparatePicrustAGORA/'+files[i])
                fluxlist = []
                for line in inFI:
                    fluxlist.append(float(line.strip()))
                inFI.close()
                #print(len(massiveArr))
                #print(firstSpecies-1)
                #print(alllist.index(samplename))
                if samplename in alllist:
                    massiveArr[alllist.index(samplename)][firstSpecies-1] = fluxlist
                print(i)
    for j in range(773):
        print(j)
        rxns = []
        rxnnames = []
        subsystems = []
        inFI = open('/mnt/vdb/home/ubuntu2/speciesSep'+str(j+1)+'Rxns.txt')
        for line in inFI:
            words = line.strip().split('\t')
            rxns.append(words[0])
            rxnnames.append(words[1])
            if len(words)==3:
                subsystems.append(words[2])
            else:
                subsystems.append('None')
        inFI.close()
        obesefluxes = []
        normalfluxes = []
        outFI = open('/mnt/vdb/home/ubuntu2/speciesSep'+str(j+1)+'RxnsFluxes.txt','w')
        outFI.write('rxn\tnormmean\tobesemean\tt-test pval\n')
        for i in range(len(alllist)):
            if i < len(obeselist) and len(massiveArr[i][j])!=0:
                obesefluxes.append(massiveArr[i][j])
            if i >= len(obeselist) and len(massiveArr[i][j])!=0:
                normalfluxes.append(massiveArr[i][j])
        if len(obesefluxes)!=0 and len(normalfluxes)!=0 and len(obesefluxes[0])!=0 and len(normalfluxes[0])!=0:
            for k in range(len(obesefluxes[0])):
                obesefluxes1 = []
                normalfluxes1 = []
                for i in range(len(obesefluxes)):
                    obesefluxes1.append(obesefluxes[i][k])
                for i in range(len(normalfluxes)):
                    normalfluxes1.append(normalfluxes[i][k])
                if rxnnames[k]=='Formate exchange':
                    formatespeciesmap[j] = np.mean(normalfluxes1) - np.mean(obesefluxes1)
                    for i in range(len(obesefluxes)):
                        formatematrix[i].append(obesefluxes[i][k])
                    for i in range(len(normalfluxes)):
                        formatematrix[i+len(obesefluxes)].append(normalfluxes[i][k])
                sumnormal = 0
                sumsquarednormal = 0
                for i in range(len(normalfluxes1)):
                    sumnormal += normalfluxes1[i]
                    sumsquarednormal += normalfluxes1[i]*normalfluxes1[i]
                sumobese = 0
                sumsquaredobese = 0
                for i in range(len(obesefluxes1)):
                    sumobese += obesefluxes1[i]
                    sumsquaredobese += obesefluxes1[i]*obesefluxes1[i]
                numer = (sumsquarednormal - sumnormal*sumnormal/len(normalfluxes1) + sumsquaredobese - sumobese*sumobese/len(obesefluxes1))
                denom = (sumsquarednormal + sumsquaredobese - (sumnormal+sumobese)*(sumnormal+sumobese)/(len(normalfluxes1)+len(obesefluxes1)))
                if denom==0:
                    metric = 'nan'
                else:
                    metric = str(numer/denom)
                if float(metric) < .05:
                    print(normalfluxes1)
                    print(obesefluxes1)
                    #nonsense = nonsense+1
                [tval,pval] = stats.ttest_ind(obesefluxes1,normalfluxes1,equal_var=False)
                #print(rxns[k])
                outFI.write('\t'.join([rxns[k],rxnnames[k],subsystems[k],str(np.mean(normalfluxes1)),str(np.mean(obesefluxes1)),str(pval),metric])+'\n')
        for i in range(len(formatematrix)):
            while len(formatematrix[i])<=j+1:
                formatematrix[i].append(0)
        outFI.close()
    #nonsense = nonsense+1
                #if math.isnan(tval):
                #    print(obesefluxes1)
                #    print(normalfluxes1)
                #    nonsense = nonsense+1
        #if len(words) >= 2:
            #wordsMap[words[0]] = ''
            #wordsMap[words[1]] = ''
    #print(len(wordsMap))

    #print(formatematrix)
    print(len(formatematrix))
    print(len(formatematrix[0]))
    matrixFI = open('/mnt/vdb/home/ubuntu2/formatematrix.txt','w')
    for i in range(len(formatematrix)):
        print(len(formatematrix[i]))
        for j in range(len(formatematrix[0])-1):
            matrixFI.write(str(formatematrix[i][j])+'\t')
        matrixFI.write(str(formatematrix[i][len(formatematrix[0])-1])+'\n')
    matrixFI.close()

    X = np.array(formatematrix)
    y = []
    for i in range(len(obeselist)):
        y.append(1)
    for i in range(len(normallist)):
        y.append(2)
    y = np.array(y)
    clf = LinearDiscriminantAnalysis()
    clf.fit(X,y)
    coefs2 = clf.coef_[0]
    coefs2.sort()
    sumcoef2 = 0
    for i in range(len(coefs2)):
        sumcoef2 += coefs2[i]*coefs2[i]
    for i in range(len(coefs2)):
        coefs2[i] = coefs2[i]/math.sqrt(sumcoef2)
    print(clf.score(X,y))

    nonsense = nonsense+1
    
    formatespecies = formatespeciesmap.keys()
    formatespecies.sort(key=lambda x: formatespeciesmap[x])
    AGORASpecies = []
    inFI = open('/mnt/vdb/home/ubuntu2/AGORASpecies.txt')
    for line in inFI:
        AGORASpecies.append(line.strip())
    inFI.close()
    outFI = open('/mnt/vdb/home/ubuntu2/formatecontributions.txt','w')
    outFI.write('xvals\tyvals\txlabels\tgrouplabels\n')
    formatespeciesmaptemp = {}
    graysum1 = 0
    graysum2 = 0
    for i in range(len(formatespecies)):
        if i < 5 or i > len(formatespecies)-6:
            formatespeciesmaptemp[formatespecies[i]] = formatespeciesmap[formatespecies[i]]
        else:
            if formatespeciesmap[formatespecies[i]] >= 0:
                graysum1 += formatespeciesmap[formatespecies[i]]
            else:
                graysum2 += formatespeciesmap[formatespecies[i]]
    formatespeciesmap = formatespeciesmaptemp
    formatespecies = formatespeciesmap.keys()
    formatespecies.sort(key=lambda x: formatespeciesmap[x])
    for formateithspecies in formatespecies:
        if formatespeciesmap[formateithspecies] < 0:
            outFI.write('2\t'+str(formatespeciesmap[formateithspecies])+'\t'+'Negative difference'+'\t'+AGORASpecies[formateithspecies]+'\n')
        else:
            outFI.write('1\t'+str(formatespeciesmap[formateithspecies])+'\t'+'Positive difference'+'\t'+AGORASpecies[formateithspecies]+'\n')
    outFI.write('1\t'+str(graysum1)+'\t'+'Positive difference'+'\t'+'Remaining Positive Difference Species'+'\n')
    outFI.write('2\t'+str(graysum2)+'\t'+'Negative difference'+'\t'+'Remaining Negative Difference Species'+'\n')
    outFI.close()
if False:
    rxnsToMeanPVal = {}
    rxnsToSubs = {}
    for i in range(773):
        inFI = open('/mnt/vdb/home/ubuntu2/speciesSep'+str(i+1)+'RxnsFluxes.txt')
        inFI.readline()
        for line in inFI:
            words = line.strip().split('\t')
            rxnsToSubs[words[0]] = words[2]
            if words[5]!='nan':
                if words[0] not in rxnsToMeanPVal:
                    rxnsToMeanPVal[words[0]] = [float(words[5])]
                else:
                    rxnsToMeanPVal[words[0]].append(float(words[5]))
    for rxn in rxnsToMeanPVal:
        rxnsToMeanPVal[rxn] = np.mean(rxnsToMeanPVal[rxn])
    subsToMeanPVal = {}
    for rxn in rxnsToMeanPVal:
        if rxnsToSubs[rxn] not in subsToMeanPVal:
            subsToMeanPVal[rxnsToSubs[rxn]] = [rxnsToMeanPVal[rxn]]
        else:
            subsToMeanPVal[rxnsToSubs[rxn]].append(rxnsToMeanPVal[rxn])
    for sub in subsToMeanPVal:
        subsToMeanPVal[sub] = np.mean(subsToMeanPVal[sub])
    subsSorted = subsToMeanPVal.keys()
    subsSorted.sort(key=lambda x: subsToMeanPVal[x])
    rxnsSorted = rxnsToMeanPVal.keys()
    rxnsSorted.sort(key=lambda x: rxnsToMeanPVal[x])
    outFI = open('/mnt/vdb/home/ubuntu2/sepRxnsByPVal.txt','w')
    outFI.write('rxn\tt-test pval\n')
    for rxn in rxnsSorted:
        outFI.write(rxn+'\t'+str(rxnsToMeanPVal[rxn])+'\n')
    outFI.close()
    outFI = open('/mnt/vdb/home/ubuntu2/sepSubsByPVal.txt','w')
    outFI.write('sub\tmean t-test pval\n')
    for sub in subsSorted:
        outFI.write(sub+'\t'+str(subsToMeanPVal[sub])+'\n')
    outFI.close()
