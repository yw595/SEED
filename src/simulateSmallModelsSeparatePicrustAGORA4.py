import sys
import subprocess
from biom import load_table
from optparse import OptionParser
import os

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
            
def makeSepOTUsTSV(ggidsArr,idx,justTwoSpecies,useCommon):
    origTsvFile = '/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData/'+ucrFolder+'/normalized_otus.tsv'
    if not os.path.exists(origTsvFile):
        foundsArr = [0 for i in range(len(ggidsArr))]
        tsvFile = '/mnt/vdb/home/ubuntu2/commonone_otu'+str(idx)+'.tsv'
        tsvFI = open(tsvFile,'w')
        for i in range(len(ggidsArr)):
            for k in range(len(ggidsArr[i])):
                foundsArr[i] = 1
                line = '#\n'
                tsvFI.write(line)
                line = '4453773\t.01\n'
                tsvFI.write(line)
                line = '\t'.join([str(ggidsArr[i][k]),'1.0\n'])
                tsvFI.write(line)
        tsvFI.close()
        return(foundsArr)

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
            if firstline and int(words[0]) not in ggidsArr[i]:
                if useCommon:
                    line = '\t'.join([words[0],'0.01\n'])
                tsvFI.write(line)
                firstline = False
            else:
                for i in range(len(ggidsArr)):
                    if int(words[0]) in ggidsArr[i]:
                        foundsArr[i] = 1
                        if useCommon:
                            line = '\t'.join([words[0],'1.0\n'])
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
    metaFile1 = '/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData/'+ucrFolder+'/'+prefix+'_predicted_metagenome.biom'
    metaFile2 = '/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData/'+ucrFolder+'/'+prefix+'_predicted_metagenome.tsv'
    if not os.path.exists(tsvFile):
        homedir = '/mnt/vdb/home/ubuntu2/'
        prefix = 'commonone_otu'+str(idx)
        tsvFile = homedir+prefix+'.tsv'
        biomFile = homedir+prefix+'.biom'
        metaFile1 = homedir+prefix+'_predicted_metagenome.biom'
        metaFile2 = homedir+prefix+'_predicted_metagenome.tsv'
        print(metaFile2)
    proc = subprocess.Popen(['biom','convert','-i',tsvFile,'-o',biomFile,'--table-type=OTU table','--to-hdf5'])
    proc.wait()
    print(['predict_metagenomes.py','-i',biomFile,'-o',metaFile1])
    proc = subprocess.Popen(['predict_metagenomes.py','-i',biomFile,'-o',metaFile1])
    proc.wait()
    proc = subprocess.Popen(['biom','convert','-i',metaFile1,'-o',metaFile2,'--to-tsv'])
    proc.wait()

def makeSepExprAndEC(species,ucrFolder,modelnamesArr,idx,justTwoSpecies,useCommon,forGiantModel=False):
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
    if not os.path.exists(metaFile2):
        homedir = '/mnt/vdb/home/ubuntu2/'
        prefix = 'commonone_otu'+str(idx)
        metaFile2 = homedir+prefix+'_predicted_metagenome.tsv'
    metaFI2 = open(metaFile2)
    if len(modelnamesArr)==2:
        ECFile = '/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData/'+ucrFolder+'/'+modelnamesArr[0]+modelnamesArr[1]+species+'.modelexpr'
    elif len(modelnamesArr)==1:
        ECFile = '/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData/'+ucrFolder+'/'+species+'.modelexpr'
    else:
        ECFile = '/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData/'+ucrFolder+'/'+'normalized_otus.modelexpr'
    if useCommon:
        ECFile = '/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData/IBDCommon/'+species+'.modelexpr'
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

    if not os.path.exists('/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData/'+ucrFolder+'/'):
        speciesFI = open(homedir+'speciesSep'+str(idx)+'.modelec')
    else:
        if len(modelnamesArr)==2:
            speciesFI = open('/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData/'+ucrFolder+'/'+modelnamesArr[0]+modelnamesArr[1]+species+'.modelec')
        elif len(modelnamesArr)==1:
            speciesFI = open('/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/MGMData/'+ucrFolder+'/'+species+'.modelec')
        else:
            if forGiantModel:
                speciesFI = open('/mnt/vdb/home/ubuntu2/forGiantModelECs.modelec')
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
        print(line)
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
        totalSpeciesCount = 0
        for EC in ECsToAdd:
            totalSpeciesCount += speciesECCounts[EC]
        if totalSpeciesCount!=0:
            totalexp /= totalSpeciesCount
        ECFI.write(str(totalexp)+'\n')
    #nonsense = nonsense+1
    ECFI.close()

def runFunc(aucrFolder,modelnameslist,modelnumberslist,isSpeciesMerged=None,justTwoSpecies=None,useCommon=None,IBDFlag=None):

    global ucrFolder
    ucrFolder = aucrFolder
    modelnameslist = modelnameslist.split(',')
    modelnumberslist = modelnumberslist.split(',')
    origmodelnameslist = modelnameslist
    origisSpeciesMerged = isSpeciesMerged
    if origisSpeciesMerged==None:
        origisSpeciesMerged = 'False'
    isSpeciesMerged = origisSpeciesMerged
    if isSpeciesMerged=='False':
        isSpeciesMerged = False
    else:
        isSpeciesMerged = True
    justTwoSpecies = justTwoSpecies
    if justTwoSpecies == 'False':
        justTwoSpecies = False
    else:
        justTwoSpecies = True
    useCommon = useCommon
    if useCommon == 'False':
        useCommon = False
    else:
        useCommon = True
    IBDFlag = IBDFlag
    if IBDFlag=='True':
        IBDFlag = True
    else:
        IBDFlag = False

    statusFile = '/mnt/vdb/home/ubuntu2/'+ucrFolder
    for i in range(len(modelnameslist)):
        origmodelname = modelnameslist[i]
        statusFile += origmodelname
    statusFile += origisSpeciesMerged
    statusFile += '.txt'
    statusFI = open(statusFile,'w')
    statusFI.write('1\n')
    statusFI.close()

    ggidsArr = []
    for i in range(len(modelnameslist)):
        origmodelname = modelnameslist[i]
        modelname = origmodelname
        if modelname.find('_')!=-1:
            modelname = ' '.join(modelname.split('_')[0:2])
        else:
            modelwords = modelname.split(' ')[0:2]
            if modelwords[1]=='sp.':
                modelwords = modelwords[:1]
            modelname = ' '.join(modelwords)
        modelnameslist[i] = modelname
        modelnumber = int(modelnumberslist[i])
        abundance = 0

        taxid = taxidFunc(modelname)
        if taxid==0:
            nonsense = nonsense+1;
            
        ggids = ggidsFunc(taxid)
        if len(set(ggids))!=len(ggids):
            nonsense = nonsense+1
        ggidsArr.append(ggids)

    print(ggidsArr)
    duplicates = False
    for i in range(len(ggidsArr)):
        ggids1 = ggidsArr[i]
        for j in range(len(ggidsArr)):
            if i!=j:
                ggids2 = ggidsArr[j]
                for ggid1 in ggids1:
                    if ggid1 in ggids2:
                        duplicates = True
    if duplicates:
        nonsense = nonsense+1

    foundsArr = []
    print(ggidsArr)
    for i in range(len(modelnameslist)):
        found = makeSepOTUsTSV([ggidsArr[i]],modelnumber,justTwoSpecies,useCommon)
        foundsArr.append(found)
        if not found==1:
            pass
            #nonsense = nonsense+1

    makeSepBiom(ucrFolder,foundsArr,modelnumber,justTwoSpecies)

    if isSpeciesMerged:
        speciesArr = ['speciesMerged']
    else:
        speciesArr = []
        for i in range(len(modelnameslist)):
            speciesArr.append('speciesSep'+str(modelnumberslist[i]))
    for i in range(len(speciesArr)):
        if IBDFlag:
            makeSepExprAndEC(speciesArr[i],ucrFolder,[origmodelnameslist[i]],modelnumberslist[i],justTwoSpecies,useCommon)
        else:
            makeSepExprAndEC(speciesArr[i],ucrFolder,[modelnameslist[i]],modelnumberslist[i],justTwoSpecies,useCommon)
    statusFI = open(statusFile,'w')
    statusFI.write('0\n')
    statusFI.close()


if __name__=='__main__':
    usage = "usage %prog [options] \n"
    prog = "this prog"

    parser = OptionParser(usage=usage)
    parser.add_option("--ucrFolder",help="")
    parser.add_option("--modelnameslist",help="")
    parser.add_option("--modelnumberslist",help="")
    parser.add_option("--isSpeciesMerged",help="")
    parser.add_option("--IBDFlag",default="False",help="")
    parser.add_option("--justTwoSpecies",help="")
    parser.add_option("--useCommon",default="False",help="")    

    (opts,args) = parser.parse_args()
    print(opts)
    print(args)
    ucrFolder = opts.ucrFolder
    modelnameslist = opts.modelnameslist
    modelnumberslist = opts.modelnumberslist
    isSpeciesMerged = opts.isSpeciesMerged
    justTwoSpecies = opts.justTwoSpecies
    useCommon = opts.useCommon
    IBDFlag = opts.IBDFlag
    runFunc(ucrFolder,modelnameslist,modelnumberslist,isSpeciesMerged=isSpeciesMerged,justTwoSpecies=justTwoSpecies,useCommon=useCommon,ooIBDFlag=IBDFlag)
