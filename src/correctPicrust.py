import os
import subprocess

inputDir1 = '/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData/'
inputList = os.listdir(inputDir1)
abunds1 = []
abunds2 = []
abundFI = open('/mnt/vdb/home/ubuntu2/compareSepSpeciesVsGroup.txt','w')
for inputItem in inputList:
    if inputItem.startswith('ucrC97'):
        inputDir = inputDir1+inputItem+'/'
        if False:
            inputFiles = os.listdir(inputDir)
            onlyModelFI = open(inputDir+'onlyModel.tsv','w')
            abundMap = {}
            for inputFile in inputFiles:
                if inputFile.startswith('one_otu') and inputFile.endswith('.tsv') and inputFile.find('metagenome')==-1:
                    speciesNumber = inputFile[7:len(inputFile)-4]
                    if speciesNumber!='':
                        inFI = open(inputDir+inputFile)
                        inFI.readline()
                        inFI.readline()
                        dataLine = inFI.readline()
                        words = dataLine.strip().split()
                        if words[0] not in abundMap:
                            abundMap[words[0]] = float(words[1])
                        else:
                            abundMap[words[0]] += float(words[1])
                        inFI.close()
            for speciesNum in abundMap:
                onlyModelFI.write(speciesNum+'\t'+str(abundMap[speciesNum])+'\n')
            onlyModelFI.close()
            proc = subprocess.Popen(['biom','convert','-i',inputDir+'onlyModel.tsv','-o',inputDir+'onlyModel.biom','--to-hdf5'])
            proc.wait()
            proc = subprocess.Popen(['predict_metagenomes.py','-i',inputDir+'onlyModel.biom','-o',inputDir+'onlyModel_predicted_metagenome.biom'])
            proc.wait()
            proc = subprocess.Popen(['biom','convert','-i',inputDir+'onlyModel_predicted_metagenome.biom','-o',inputDir+'onlyModel_predicted_metagenome.tsv','--to-tsv'])
            proc.wait()

        if True and os.path.exists(inputDir+'onlyModel_predicted_metagenome.tsv'):
            print(inputDir)
            inFI = open(inputDir+'onlyModel_predicted_metagenome.tsv')
            KOAbund1 = {}
            inFI.readline()
            inFI.readline()
            for line in inFI:
                words = line.split()
                KOAbund1[words[0]] = float(words[1])
            inFI.close()

            inputFiles = os.listdir(inputDir)
            KOAbund2 = {}
            for inputFile in inputFiles:
                if inputFile.startswith('one_otu') and inputFile.endswith('.tsv') and inputFile.find('metagenome')!=-1 and inputFile.find('test')==-1:
                    speciesNumber = inputFile[7:len(inputFile)-26]
                    if speciesNumber!='':
                        inFI = open(inputDir+inputFile)
                        inFI.readline()
                        inFI.readline()
                        for line in inFI:
                            words = line.strip().split()
                            if words[0] not in KOAbund2:
                                KOAbund2[words[0]] = float(words[1])
                            else:
                                KOAbund2[words[0]] += float(words[1])
                        inFI.close()

            abund1 = 0
            abund2 = 0
            for KO in KOAbund1:
                abund1 += KOAbund1[KO]
                abund2 += KOAbund2[KO]
                #print(str(KOAbund1[KO])+' 1')
                #print(str(KOAbund2[KO])+' 2')
            abunds1.append(abund1)
            abunds2.append(abund2)
    print(abunds1)
    print(abunds2)
    abundFI.write(str(abund1)+'\t'+str(abund2)+'\n')
abundFI.close()
