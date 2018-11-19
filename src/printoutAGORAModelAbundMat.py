import sys
import os
import copy
from bootstrapAGORAIBD import writeData

FI1 = open('/mnt/vdb/home/ubuntu2/justAGORADistsManual.txt')
numbers1 = []
numToSpec = {}
for line in FI1:
    words = line.split('|')
    if words[1]!='':
        numbers1.append(int(words[1]))
    numToSpec[int(words[1])] = words[0]
FI1.close()

ggToFI = open('/mnt/vdb/home/ubuntu2/ggIDsToTax.txt')
ggTo = {}
for line in ggToFI:
    words = line.split()
    ggTo[int(words[0])] = int(words[1])
ggToFI.close()

xmlfiles = []
allmodelfiles = os.listdir('/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/AGORAModels/Western-Diet-Paper')
for i in range(len(allmodelfiles)):
    if allmodelfiles[i].endswith('.xml'):
        xmlfiles.append(allmodelfiles[i])

filenames = os.listdir('/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData')
samplenamearrall = []
abundsarrall = []
speciesnamearrall = []
for filename in filenames:
    if filename.startswith('ucrC97'):
        normOTU = '/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData/'+filename+'/normalized_otus.tsv'
        if os.path.exists(normOTU):
            samplenamearr = [filename[6:] for ark in range(len(xmlfiles))]
            abundsarr = [0 for ark in range(len(xmlfiles))]
            speciesnamearr = copy.deepcopy(xmlfiles)
            FI2 = open(normOTU)
            FI2.readline()
            FI2.readline()
            taxonids = []
            taxonabunds = []
            for line in FI2:
                taxonid = int(line.split()[0])
                taxonabund = float(line.split()[1])
                if taxonid in ggTo and ggTo[taxonid] in numbers1:
                    speciesid = numToSpec[ggTo[taxonid]]
                    for i in range(len(speciesnamearr)):
                        specieswords1 = speciesnamearr[i].split('_')
                        specieswords2 = speciesid.split(' ')
                        matchesenough = True
                        for j in range(min(len(specieswords1),len(specieswords2))):
                            if specieswords1[j]!=specieswords2[j]:
                                matchesenough = False
                        if matchesenough:
                            abundsarr[i] += taxonabund
            FI2.close()
            speciesnamearrall.extend(speciesnamearr)
            abundsarrall.extend(abundsarr)
            samplenamearrall.extend(samplenamearr)
            print(filename)

pickFI = open('/mnt/vdb/home/ubuntu2/normObeseVec.txt','w')
writeData([samplenamearrall,speciesnamearrall,abundsarrall],'/mnt/vdb/home/ubuntu2/normObeseVec.txt',delimiter='\t',headers=['sample','species','abund'])
