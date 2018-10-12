import sys
import subprocess
from biom import load_table
from optparse import OptionParser

def taxidFunc(modelname):
    justFI = open('/mnt/vdb/home/ubuntu2/justAGORADistsManual.txt')
    taxid = 0
    for line in justFI:
        words = line.split('|')
        if words[0] == modelname:
            taxid = int(words[1])
    justFI.close()
    return(taxid)

def ggidsFunc(taxid):
    ggIDsToTaxFI = open('/mnt/vdb/home/ubuntu2/ggIDsToTax.txt')
    ggids = []
    for line in ggIDsToTaxFI:
        words = line.split()
        if int(words[1])==taxid:
            ggids.append(int(words[0]))
    return(ggids)
            
def makeSepOTUsTSV(ggidsArr,idx,justTwoSpecies):
    origTsvFile = '/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData/'+ucrFolder+'/normalized_otus.tsv'
    origTsvFI = open(origTsvFile)
    if justTwoSpecies:
        tsvFile = '/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData/'+ucrFolder+'/two_otus_test'+str(idx)+'.tsv'
    else:
        if len(ggidsArr)==2:
            tsvFile = '/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData/'+ucrFolder+'/two_otus'+str(idx)+'.tsv'
        else:
            tsvFile = '/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData/'+ucrFolder+'/one_otu'+str(idx)+'.tsv'
    tsvFI = open(tsvFile,'w')
    #origTsvFI.readline()
    #origTsvFI.readline()
    foundsArr = [0 for i in range(len(ggidsArr))]
    abundance = 0
    firstline = True
    for line in origTsvFI:
        if line.startswith('#'):
            tsvFI.write(line)
        else:
            words = line.split()
            if firstline:
                tsvFI.write(line)
                firstline = False
            else:
                for i in range(len(ggidsArr)):
                    if int(words[0]) in ggidsArr[i]:
                        foundsArr[i] = 1
                        tsvFI.write(line)
                        abundance += float(words[1])
    print('Result ' + str(abundance))
    origTsvFI.close()
    tsvFI.close()
    return(foundsArr)

def makeSepBiom(ucrFolder,foundsArr,idx,justTwoSpecies):
    if justTwoSpecies:
        prefix = 'two_otus_test'+str(idx)
    else:
        if len(foundsArr)==2:
            prefix = 'two_otus'+str(idx)
        elif len(foundsArr)==1:
            prefix = 'one_otu'+str(idx)
        else:
            prefix = 'normalized_otus'
    biomFile = '/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData/'+ucrFolder+'/'+prefix+'.biom'
    #table = load_table(tsvFile)
    #biomFI = open(biomFile,'w')
    #biomFI.write(str(table))
    tsvFile = '/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData/'+ucrFolder+'/'+prefix+'.tsv'
    proc = subprocess.Popen(['biom','convert','-i',tsvFile,'-o',biomFile,'--table-type=OTU table','--to-hdf5'])
    proc.wait()
    metaFile1 = '/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData/'+ucrFolder+'/'+prefix+'_predicted_metagenome.biom'
    metaFile2 = '/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData/'+ucrFolder+'/'+prefix+'_predicted_metagenome.tsv'
    proc = subprocess.Popen(['predict_metagenomes.py','-i',biomFile,'-o',metaFile1])
    proc.wait()
    proc = subprocess.Popen(['biom','convert','-i',metaFile1,'-o',metaFile2,'--to-tsv'])
    proc.wait()

def makeSepExprAndEC(species,ucrFolder,modelnamesArr,idx,justTwoSpecies):
    kegg2ECFI = open('/mnt/vdb/home/ubuntu2/kegg2EC.txt')
    kegg2EC = {}
    for line in kegg2ECFI:
        words = line.strip().split()
        kegg2EC[words[0]] = words[1:]
    kegg2ECFI.close()

    if justTwoSpecies:
        prefix = 'two_otus_test'+str(idx)
    else:
        if len(modelnamesArr)==2:
            prefix = 'two_otus'+str(idx)
        elif len(modelnamesArr)==1:
            prefix = 'one_otu'+str(idx)
        else:
            prefix = 'normalized_otus'
    metaFile2 = '/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData/'+ucrFolder+'/'+prefix+'_predicted_metagenome.tsv'
    metaFI2 = open(metaFile2)
    if len(modelnamesArr)==2:
        ECFile = '/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData/'+ucrFolder+'/'+modelnamesArr[0]+modelnamesArr[1]+species+'.modelexpr'
    elif len(modelnamesArr)==1:
        ECFile = '/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData/'+ucrFolder+'/'+species+'.modelexpr'
    else:
        ECFile = '/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData/'+ucrFolder+'/'+'normalized_otus.modelexpr'
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

    if len(modelnamesArr)==2:
        speciesFI = open('/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData/'+ucrFolder+'/'+modelnamesArr[0]+modelnamesArr[1]+species+'.modelec')
    elif len(modelnamesArr)==1:
        speciesFI = open('/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData/'+ucrFolder+'/'+species+'.modelec')
    else:
        speciesFI = open('/mnt/vdb/home/ubuntu2/allOxECs.txt')
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
        for EC in ECsToAdd:
            if EC in EC2Expr:
                totalexp += EC2Expr[EC]
        totalexp /= speciesECCounts[EC]
        ECFI.write(str(totalexp)+'\n')
    ECFI.close()

if __name__=='__main__':
    usage = "usage %prog [options] \n"
    prog = "this prog"

    parser = OptionParser(usage=usage)
    parser.add_option("--isTwoSpecies",help="")
    parser.add_option("--ucrFolder",help="")
    parser.add_option("--modelname1",help="")
    parser.add_option("--modelname2",help="")
    parser.add_option("--isSpecies3",help="")
    parser.add_option("--IBDFlag",default="False",help="")
    parser.add_option("--modelname",help="")
    parser.add_option("--modelnumber",help="")
    parser.add_option("--justTwoSpecies",help="")

    (opts,args) = parser.parse_args()
    print(opts)
    print(args)
    
    isTwoSpecies = opts.isTwoSpecies
    ucrFolder = opts.ucrFolder
    if isTwoSpecies=='False':
        isTwoSpecies = False
    else:
        isTwoSpecies = True
    if isTwoSpecies:
        modelname1 = opts.modelname1
        origmodelname1 = modelname1
        modelname1 = ' '.join(modelname1.split('_')[0:2])
        modelname2 = opts.modelname2
        origmodelname2 = modelname2
        modelname2 = ' '.join(modelname2.split('_')[0:2])
        isSpecies3 = opts.isSpecies3
        origisSpecies3 = isSpecies3
        if isSpecies3=='False':
            isSpecies3 = False
        else:
            isSpecies3 = True
        IBDFlag = opts.IBDFlag
        if IBDFlag=='True':
            IBDFlag = True
        statusFI = open('/mnt/vdb/home/ubuntu2/'+ucrFolder+origmodelname1+origmodelname2+origisSpecies3+'.txt','w')
        statusFI.write('1\n')
        statusFI.close()

        taxid1 = taxidFunc(modelname1)
        taxid2 = taxidFunc(modelname2)

        if taxid1==0 or taxid2==0:
            nonsense = nonsense+1;

        ggids1 = ggidsFunc(taxid1)
        ggids2 = ggidsFunc(taxid2)

        if len(set(ggids1))!=len(ggids1) or len(set(ggids2))!=len(ggids2):
            nonsense = nonsense+1
        duplicates = False
        for ggid1 in ggids1:
            if ggid1 in ggids2:
                duplicates = True
        if duplicates:
            nonsense = nonsense+1
        
        foundsArr = makeSepOTUsTSV([ggids1,ggids2],0,False)
        found1 = foundsArr[0]==1
        found2 = foundsArr[1]==1

        if not found1 or not found2:
            pass
            #nonsense = nonsense+1

        makeSepBiom(ucrFolder,foundsArr,0,False)

        if isSpecies3:
            speciesArr = ['species3']
        else:
            speciesArr = ['species1','species2']
        for species in speciesArr:
            if IBDFlag:
                makeSepExprAndEC(species,ucrFolder,[origmodelname1,origmodelname2],0,False)
            else:
                makeSepExprAndEC(species,ucrFolder,[modelname1,modelname2],0,False)
        statusFI = open('/mnt/vdb/home/ubuntu2/'+ucrFolder+origmodelname1+origmodelname2+origisSpecies3+'.txt','w')
        statusFI.write('0\n')
        statusFI.close()
    else:
        modelname = opts.modelname
        modelnumber = int(opts.modelnumber)
        justTwoSpecies = opts.justTwoSpecies
        if justTwoSpecies == 'False':
            justTwoSpecies = False
        else:
            justTwoSpecies = True
        IBDFlag = opts.IBDFlag
        if IBDFlag=='True':
            IBDFlag = True
        origmodelname = modelname
        modelname = ' '.join(modelname.split('_')[0:2])
        abundance = 0
        statusFI = open('/mnt/vdb/home/ubuntu2/'+ucrFolder+origmodelname+'.txt','w')
        statusFI.write('1\n')
        statusFI.close()

        taxid = taxidFunc(modelname)
        if taxid==0:
            nonsense = nonsense+1;

        print(taxid)
        ggids = ggidsFunc(taxid)

        print(justTwoSpecies)
        print(ucrFolder)
        foundsArr = makeSepOTUsTSV([ggids],modelnumber,justTwoSpecies)
        found = foundsArr[0]==1

        if not found:
            pass
            #nonsense = nonsense+1

        makeSepBiom(ucrFolder,foundsArr,modelnumber,justTwoSpecies)

        speciesArr = ['speciesSep'+str(modelnumber)]
        for species in speciesArr:
            makeSepExprAndEC(species,ucrFolder,[modelname],modelnumber,justTwoSpecies)
        statusFI = open('/mnt/vdb/home/ubuntu2/'+ucrFolder+origmodelname+'.txt','w')
        statusFI.write('0\n')
        statusFI.close()
