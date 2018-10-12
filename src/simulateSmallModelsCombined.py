import cobra
import cobra.test
import time
from bootstrapAGORAIBD import writeData
from computeMinDisj import computeMinDisj
from falconMulti import falconMulti
import numpy as np
import scipy.stats
from runFALCONStripped import convertIrrevFluxDistribution
from runFALCONStripped import runFALCONStripped
import os
import copy
import multiprocessing
from simulateSmallModelsSeparatePicrustAGORA4 import runFunc

class simulateSmallModelsCombinedClass(object):

    def __init__(self):
        self.manager = multiprocessing.Manager()
        self.AGORAModelArr = []#self.manager.list()
        self.AGORAModelArrTemp = []
        files = os.listdir('/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/AGORAModels/Western-Diet-Paper')
        files.sort()
        for i in range(10):#len(files)):
            if files[i].endswith('.xml'):
                print(files[i])
                print(i)
                ithModel = cobra.io.read_sbml_model('/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/AGORAModels/Western-Diet-Paper/'+files[i])
                for k in range(len(ithModel.reactions)):
                    ithModel.reactions[k].gene_reaction_rule = ithModel.reactions[k].id
                selectgenenames = [rxn.id for rxn in ithModel.reactions]
                unselectgenes = []
                for k in range(len(ithModel.genes)):
                    if ithModel.genes[k].id not in selectgenenames:
                        unselectgenes.append(ithModel.genes[k])
                cobra.manipulation.delete.remove_genes(ithModel,unselectgenes)
                self.AGORAModelArrTemp.append(ithModel)
        self.AGORAModelArrTemp.sort(key=lambda x: x.description)
        for i in range(len(self.AGORAModelArrTemp)):
            self.AGORAModelArr.append(self.AGORAModelArrTemp[i])
        self.subsystemsArr = []
        for i in range(len(self.AGORAModelArr)):
            inFI = open('/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData/CommonSubsystems/species'+str(i+1)+'Subsystems.txt')
            ithSubsystems = []
            for line in inFI:
                ithSubsystems.append(line.strip())
            self.subsystemsArr.append(ithSubsystems)
            inFI.close()
        self.allSamples = []
        allSpecies = []
        inFI = open('/mnt/vdb/home/ubuntu2/diabetesVec.txt')
        inFI.readline()
        for line in inFI:
            line = line.strip()
            words = line.split('\t')
            self.allSamples.append(words[0])
            allSpecies.append(words[1])
        inFI.close()
        self.allSamples = np.unique(self.allSamples)
        self.allSamples.sort()
        allSamples1 = {}
        for i in range(len(self.allSamples)):
            allSamples1[self.allSamples[i]] = i
        self.allSamples = allSamples1
        self.allSamplesArr = self.allSamples.keys()
        self.allSamplesArr.sort()
        allSpecies = np.unique(allSpecies)
        allSpecies.sort()
        allSpecies1 = {}
        for i in range(len(allSpecies)):
            allSpecies1[allSpecies[i]] = i
        allSpecies = allSpecies1
        self.AGORAModelAbundMat = np.zeros([len(self.allSamples),len(allSpecies)])
        inFI = open('/mnt/vdb/home/ubuntu2/diabetesVec.txt')
        inFI.readline()
        for line in inFI:
            line = line.strip()
            words = line.split('\t')
            self.AGORAModelAbundMat[self.allSamples[words[0]],allSpecies[words[1]]] = float(words[2])
        inFI.close()

    def runFunc2(self,i,j):
        print(str(i)+' '+str(j))
    def runFunc(self,i,j):
        for z in [0]:
            if i<j:# and self.AGORAModelAbundMat[z][i]!=0 and self.AGORAModelAbundMat[z][j]!=0:
                for checkIndex in [i,j]:
                    if True:#not os.path.exists('/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData/IBDCommon/speciesSep'+str(checkIndex+1)+'.modelexpr'):
                        eclist = []
                        for k in range(len(self.AGORAModelArr[i].reactions)):
                            if 'ec-code' in self.AGORAModelArr[i].reactions[k].annotation:
                                eclist.append(self.AGORAModelArr[i].reactions[k].annotation['ec-code'])
                            else:
                                eclist.append('')
                        writeData([eclist],'/mnt/vdb/home/ubuntu2/speciesSep'+str(checkIndex+1)+'.modelec',delimiter='\t')
                        runFunc(self.allSamplesArr[z],self.AGORAModelArr[checkIndex].description,str(i+1))

                inFI = open('/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData/IBDCommon/speciesSep'+str(i+1)+'.modelexpr')
                ithExpr = []
                #print('continue')
                for line in inFI:
                    ithExpr.append(self.AGORAModelAbundMat[z][i]*float(line.strip()))
                inFI.close()
                inFI = open('/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData/IBDCommon/speciesSep'+str(j+1)+'.modelexpr')
                jthExpr = []
                for line in inFI:
                    jthExpr.append(self.AGORAModelAbundMat[z][j]*float(line.strip()))
                inFI.close()
                pairmodelith = copy.deepcopy(self.AGORAModelArr[i])
                for k in range(len(pairmodelith.reactions)):
                    if not self.subsystemsArr[i][k]=='Exchange/demand reaction' or not pairmodelith.reactions[k].id.startswith('EX'):
                        pairmodelith.reactions[k].id = '1_'+pairmodelith.reactions[k].id
                        pass
                for k in range(len(pairmodelith.metabolites)):
                    if not pairmodelith.metabolites[k].id.startswith('_e'):
                        pairmodelith.metabolites[k].id = '1_'+pairmodelith.metabolites[k].id
                        pass
                pairmodeljth = copy.deepcopy(self.AGORAModelArr[j])
                for k in range(len(pairmodeljth.reactions)):
                    if not self.subsystemsArr[j][k]=='Exchange/demand reaction' or not pairmodeljth.reactions[k].id.startswith('EX'):
                        pairmodeljth.reactions[k].id = '2_'+pairmodeljth.reactions[k].id
                        pass
                for k in range(len(pairmodeljth.metabolites)):
                    if not pairmodeljth.metabolites[k].id.startswith('_e'):
                        pairmodeljth.metabolites[k].id = '2_'+pairmodeljth.metabolites[k].id
                        pass
                pairmodel = copy.deepcopy(pairmodelith)
                pairexpr = copy.deepcopy(ithExpr)
                for k in range(len(pairmodeljth.reactions)):
                    if not self.subsystemsArr[j][k]=='Exchange/demand reaction' or not pairmodeljth.reactions[k].id.startswith('EX'):
                        pairmodel.add_reactions([pairmodeljth.reactions[k]])
                        pairexpr.append(jthExpr[k])
                        pass
                # for k in range(len(self.subsystemsArr[j])):
                #     if not self.subsystemsArr[j][k]=='Exchange/demand reaction' or not self.AGORAModelArr[j].reactions[k].id.startswith('EX'):
                #         kthReaction = self.AGORAModelArr[j].reactions[k]
                #         kthReaction.id = '2_'+kthReaction.id
                #         pairexpr.append(jthExpr[k])
                #         pairmodel.add_reactions([kthReaction])
                expressionIDs = [rxn.id for rxn in pairmodel.reactions]
                expressionSDs = [1 for rxn in pairmodel.reactions]
                for k in range(len(pairmodel.reactions)):
                    pairmodel.reactions[k].gene_reaction_rule = pairmodel.reactions[k].id
                selectgenenames = [rxn.id for rxn in pairmodel.reactions]
                unselectgenes = []
                for k in range(len(pairmodel.genes)):
                    if pairmodel.genes[k].id not in selectgenenames:
                        unselectgenes.append(pairmodel.genes[k])
                cobra.manipulation.delete.remove_genes(pairmodel,unselectgenes)
                nReps = 1
                #try:
                [v_falcon, objValue, cost_rev, corr_rho, rxn_exp_md_rev] = runFALCONStripped(pairmodel,expressionIDs,pairexpr,expressionSDs,nReps,ani=i,aj=j)
                rxnids = [rxn.id for rxn in pairmodel.reactions]
                writeData([v_falcon,pairexpr,rxnids],'/mnt/vdb/home/ubuntu2/MATLAB/SEED/output/simulateSmallModelsCombined/outputfile_'+str(i)+'_'+str(j)+'.txt',delimiter='\t',headers=['flux','expr','rxnid'])
                #except:
                    #pass
                    #self.twofluxesarr[i].append([])
                    #self.twoexprarr[i].append([])
                    #self.tworxnsarr[i].append([])

    def printSummary(self):
        corr1 = []
        corr2 = []
        intakeArr = []
        outputArr = []
        wholeArr = []
        for i in range(len(self.twofluxesarr)):
            for j in range(len(self.twofluxesarr[i])):
                if sum(self.twofluxesarr[i][j])!=0 and len(self.twofluxesarr[i][j])!=0:
                    for k in range(len(self.twofluxesarr[i][j])):
                        if self.tworxnsarr[i][j][k].startswith('EX'):
                            wholeArr.append(self.twofluxesarr[i][j][k])
                            if self.twofluxesarr[i][j][k]<0:
                                intakeArr.append(self.twofluxesarr[i][j][k])
                            elif self.twofluxesarr[i][j][k]>0:
                                outputArr.append(self.twofluxesarr[i][j][k])
                    for k in range(len(self.twofluxesarr[i][j])):
                        corr1.append(self.twofluxesarr[i][j][k])
                        corr2.append(self.twoexprarr[i][j][k])
        writeData([corr1,corr2],'/mnt/vdb/home/ubuntu2/pythonSmallModelsFluxExprCorr.txt',headers=['flux','expr'])
        writeData([wholeArr],'/mnt/vdb/home/ubuntu2/pythonInteractionDist.txt',headers=['interaction'])
        #nonsense = nonsense+1

    def readResults(self):
        self.twofluxesarr = []
        self.twoexprarr = []
        self.tworxnsarr = []
        for i in range(len(self.AGORAModelArr)):
            self.twofluxesarr.append([])
            self.twoexprarr.append([])
            self.tworxnsarr.append([])
            for j in range(len(self.AGORAModelArr)):
                self.twofluxesarr[i].append([])
                self.twoexprarr[i].append([])
                self.tworxnsarr[i].append([])

        for i in range(len(self.AGORAModelArr)):
            for j in range(len(self.AGORAModelArr)):
                datafile = '/mnt/vdb/home/ubuntu2/MATLAB/SEED/output/simulateSmallModelsCombined/outputfile_'+str(i)+'_'+str(j)+'.txt'
                twofluxes = []
                twoexpr = []
                tworxns = []
                if os.path.exists(datafile):
                    dataFI = open(datafile)
                    dataFI.readline()
                    for line in dataFI:
                        words = line.strip().split('\t')
                        twofluxes.append(float(words[0]))
                        twoexpr.append(float(words[1]))
                        tworxns.append(words[2])
                self.twofluxesarr[i][j] = twofluxes
                self.twoexprarr[i][j] = twoexpr
                self.tworxnsarr[i][j] = tworxns
                
        
if __name__ == '__main__':
    obj = simulateSmallModelsCombinedClass()
    prlist = []
    obj.runFunc(3,4)
    # for i in range(len(obj.AGORAModelArr)):
    #     for j in range(len(obj.AGORAModelArr)):
    #         pr = multiprocessing.Process(target=obj.runFunc, args=(i,j,))
    #         prlist.append(pr)
    # for pr in prlist:
    #     print('start')
    #     pr.start()
    # obj.readResults()
    # nonsense = nonsense+1
