import sys
import subprocess
from biom import load_table

ucrFolder = sys.argv[1]
modelnames = sys.argv[2].split(',')
abundance = 0
totalabundance = 0
for modelname in modelnames:
    
    modelname = ' '.join(modelname.split('_')[0:2])
    justFI = open('/mnt/vdb/home/ubuntu2/justAGORADistsManual.txt')
    taxid = 0
    for line in justFI:
        words = line.split('|')
        if words[0] == modelname:
            taxid = int(words[1])
    justFI.close()

    ggids = []
    ggIDsToTaxFI = open('/mnt/vdb/home/ubuntu2/ggIDsToTax.txt')
    for line in ggIDsToTaxFI:
        words = line.split()
        if int(words[1])==taxid:
            ggids.append(int(words[0]))

    origTsvFile = '/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData/'+ucrFolder+'/normalized_otus.tsv'
    origTsvFI = open(origTsvFile)
    origTsvFI.readline()
    origTsvFI.readline()
    for line in origTsvFI:
        words = line.split()
        totalabundance += float(words[1])
        if int(words[0]) in ggids:
            abundance += float(words[1])
    origTsvFI.close()
print('Result ' + str(abundance) + ' ' + str(totalabundance))
