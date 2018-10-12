import sys
import subprocess
from biom import load_table

ucrFolder = sys.argv[1]
modelname1 = sys.argv[2]
modelname1 = ' '.join(modelname1.split('_')[0:2])
modelname2 = sys.argv[3]
modelname2 = ' '.join(modelname2.split('_')[0:2])
isSpecies3 = sys.argv[4]
if isSpecies3=='False':
    isSpecies3 = False
else:
    isSpecies3 = True
statusFI = open('/mnt/xvdf/home/ubuntu2/'+sys.argv[1]+sys.argv[2]+sys.argv[3]+sys.argv[4]+'.txt','w')
statusFI.write('1\n')
statusFI.close()

if isSpecies3=='True':
    isSpecies3 = True

print(modelname1)
print(modelname2)
justFI = open('/mnt/xvdf/home/ubuntu2/justAGORADistsManual.txt')
taxid1 = 0
taxid2 = 0
for line in justFI:
    words = line.split('|')
    if words[0] == modelname1:
        taxid1 = int(words[1])
    if words[0] == modelname2:
        taxid2 = int(words[1])
justFI.close()

print(taxid1)
print(taxid2)
if taxid1==0 or taxid2==0:
    nonsense = nonsense+1;

#taxid1 = 472834#590458
#taxid2 = 136703#314824
ggids1 = []
ggids2 = []
ggIDsToTaxFI = open('/mnt/xvdf/home/ubuntu2/ggIDsToTax.txt')
for line in ggIDsToTaxFI:
    words = line.split()
    if int(words[1])==taxid1:
        ggids1.append(int(words[0]))
    if int(words[1])==taxid2:
        ggids2.append(int(words[0]))

origTsvFile = '/mnt/xvdf/home/ubuntu2/MATLAB/SEED/input/MGMData/'+ucrFolder+'/normalized_otus.tsv'
origTsvFI = open(origTsvFile)
tsvFile = '/mnt/xvdf/home/ubuntu2/MATLAB/SEED/input/MGMData/'+ucrFolder+'/two_otus.tsv'
tsvFI = open(tsvFile,'w')
origTsvFI.readline()
origTsvFI.readline()
found1 = False
found2 = False
for line in origTsvFI:
    words = line.split()
    if int(words[0]) in ggids1:
        found1 = True
        tsvFI.write(line)
    if int(words[0]) in ggids2:
        found2 = True
        tsvFI.write(line)
origTsvFI.close()
tsvFI.close()

if not found1 or not found2:
    pass
    #nonsense = nonsense+1

biomFile = '/mnt/xvdf/home/ubuntu2/MATLAB/SEED/input/MGMData/'+ucrFolder+'/two_otus.biom'
#table = load_table(tsvFile)
#biomFI = open(biomFile,'w')
#biomFI.write(str(table))
proc = subprocess.Popen(['biom','convert','-i',tsvFile,'-o',biomFile,'--table-type=OTU table','--to-hdf5'])
proc.wait()
metaFile1 = '/mnt/xvdf/home/ubuntu2/MATLAB/SEED/input/MGMData/'+ucrFolder+'/two_otus_predicted_metagenome.biom'
metaFile2 = '/mnt/xvdf/home/ubuntu2/MATLAB/SEED/input/MGMData/'+ucrFolder+'/two_otus_predicted_metagenome.tsv'
proc = subprocess.Popen(['predict_metagenomes.py','-i',biomFile,'-o',metaFile1])
proc.wait()
proc = subprocess.Popen(['biom','convert','-i',metaFile1,'-o',metaFile2,'--to-tsv'])
proc.wait()

kegg2ECFI = open('/mnt/xvdf/home/ubuntu2/kegg2EC.txt')
kegg2EC = {}
for line in kegg2ECFI:
    words = line.strip().split()
    kegg2EC[words[0]] = words[1:]
kegg2ECFI.close()

if isSpecies3:
    speciesArr = ['species3']
else:
    speciesArr = ['species1','species2']
for species in speciesArr:
    metaFile2 = '/mnt/xvdf/home/ubuntu2/MATLAB/SEED/input/MGMData/'+ucrFolder+'/two_otus_predicted_metagenome.tsv'
    metaFI2 = open(metaFile2)
    ECFile = '/mnt/xvdf/home/ubuntu2/MATLAB/SEED/input/MGMData/'+ucrFolder+'/'+sys.argv[2]+sys.argv[3]+species+'.modelexpr'
    ECFI = open(ECFile,'w')
    EC2Expr = {}
    for line in metaFI2:
        words = line.split()
        if words[0] in kegg2EC:
            ECs = kegg2EC[words[0]]
            for EC in ECs:
                expr = float(words[1])/len(ECs)
                if EC in EC2Expr:
                    EC2Expr[EC] = EC2Expr[EC]+expr
                else:
                    EC2Expr[EC] = expr
    metaFI2.close()

    speciesFI = open('/mnt/xvdf/home/ubuntu2/MATLAB/SEED/input/MGMData/'+ucrFolder+'/'+sys.argv[2]+sys.argv[3]+species+'.modelec')
    specieslines = []
    speciesECCounts = {}
    for line in speciesFI:
        line = line.strip()
        specieslines.append(line)
        ECsToAdd = []
        if line.find(',')!=-1:
            words = line.split(',')
            for word in words:
                word = word.strip()
                ECsToAdd.append(word)
        elif line != '':
            ECsToAdd.append(line)
        for EC in ECsToAdd:
            if EC in speciesECCounts:
                speciesECCounts[EC] = speciesECCounts[EC]+1
            else:
                speciesECCounts[EC] = 1

    for line in specieslines:
        ECsToAdd = []
        if line.find(',')!=-1:
            words = line.split(',')
            for word in words:
                word = word.strip()
                ECsToAdd.append(word)
        elif line != '':
            ECsToAdd.append(line)

        totalexp = 0
        #print('HERE')
        for EC in ECsToAdd:
            if EC in EC2Expr:
                #print(EC2Expr[EC])
                totalexp += EC2Expr[EC]
        #print(totalexp)
        totalexp /= speciesECCounts[EC]
        ECFI.write(str(totalexp)+'\n')
    ECFI.close()
statusFI = open('/mnt/xvdf/home/ubuntu2/'+sys.argv[1]+sys.argv[2]+sys.argv[3]+sys.argv[4]+'.txt','w')
statusFI.write('0\n')
statusFI.close()
