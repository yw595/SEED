import cobra
import cobra.test
import time
from bootstrapAGORAIBD import writeData
from computeMinDisj import computeMinDisj
from falconMulti import falconMulti
import numpy as np
import scipy.stats
import subprocess

def convertIrrevFluxDistribution(fluxdist,modelIrrev):
    revfluxdist = []
    irrev2rev = [0 for i in range(len(modelIrrev.reactions))]
    rev2irrev = []
    for i in range(len(modelIrrev.reactions)):
        if not modelIrrev.reactions[i].id.endswith('_reverse'):
            revfluxsum = fluxdist[i]
            revfluxmatchidxs = [i]
            for j in range(len(modelIrrev.reactions)):
                if modelIrrev.reactions[j].id==modelIrrev.reactions[i].id+'_reverse':
                    revfluxsum += -fluxdist[j]
                    revfluxmatchidxs.append(j)
                    irrev2rev[j] = len(revfluxdist)
            revfluxdist.append(revfluxsum)
            rev2irrev.append(revfluxmatchidxs)
            irrev2rev[i] = len(revfluxdist)-1

    return [np.array(revfluxdist), irrev2rev, rev2irrev]

def runFALCONStripped(model,expressionIDs,expData,expressionSDs,nReps,overwriteSDs={},overwriteMeans={},ani='',aj=''):
    useMinDisj = True
    expCon = False
    minFit = 0.0
    regC = 0.0
    FDEBUG = False

    nrxns = len(model.reactions)
    cobra.manipulation.modify.convert_to_irreversible(model)
    modelIrrev = model
    modelIrrev.S = cobra.util.array.create_stoichiometric_matrix(modelIrrev)
    modelIrrev.lb = np.array([rxn.lower_bound for rxn in modelIrrev.reactions])
    modelIrrev.ub = np.array([rxn.upper_bound for rxn in modelIrrev.reactions])
    modelIrrev.rxnGeneMat = np.zeros([len(modelIrrev.reactions),len(modelIrrev.genes)])
    modelIrrev.c = np.array([rxn.objective_coefficient for rxn in modelIrrev.reactions])
    for i in range(len(modelIrrev.reactions)):
        for j in range(len(modelIrrev.genes)):
            for gene in modelIrrev.reactions[i].genes:
                if modelIrrev.genes[j].id==gene.id:
                    modelIrrev.rxnGeneMat[i,j]=1
    expTime = time.time()
    genedata_filename = '/mnt/vdb/home/ubuntu2/temp.csv'
    if ani!='' and aj!='':
        genedata_filename = '/mnt/vdb/home/ubuntu2/temp_'+str(i)+'_'+str(j)+'.csv'
    writeData([expressionIDs,expData,expressionSDs],genedata_filename)
    [rxn_exp_md, rxn_exp_sd_md, rxn_rule_group] = computeMinDisj(modelIrrev,genedata_filename,0,True,overwriteSDs,overwriteMeans)
    minDisjTime = time.time()-expTime
    minDisjTime = 0

    [v_falconIrr,temp1,temp2,temp3,temp4,temp5,v_falconIrr_s,temp6,temp7,temp8,temp9,fOpt,f_easyLP,v_easyLP,cost_irrev] = falconMulti(modelIrrev, nReps, rxn_exp_md, rxn_exp_sd_md, rxn_rule_group,rc=regC,minFit=minFit,EXPCON=expCon,FDEBUG=FDEBUG,ani=ani,aj=aj)
    print(v_falconIrr)
    proc = subprocess.Popen(['rm',genedata_filename])
    proc.wait()
    [v_falcon, irrev2rev, rev2irrev] = convertIrrevFluxDistribution(v_falconIrr, modelIrrev)
    print(v_falcon)
    [v_falcon_s, irrev2rev, rev2irrev] = convertIrrevFluxDistribution(v_falconIrr_s, modelIrrev)
    validIdxs = (rxn_exp_md!=float('nan'))
    objValue = fOpt
    objValueOld = sum(abs(v_falconIrr[validIdxs]-rxn_exp_md[validIdxs])/rxn_exp_sd_md[validIdxs])
    fluxSolSum = sum(abs(v_falconIrr[validIdxs]))
    rxn_exp_md_corr = rxn_exp_md[validIdxs]
    v_falconIrr_corr = abs(v_falconIrr[validIdxs])
    if len(rxn_exp_md_corr)>0:
        [corr_rho, corr_pval] = scipy.stats.spearmanr(rxn_exp_md_corr,v_falconIrr_corr)
    else:
        corr_rho = 0
    
    rxn_exp_md_rev = np.zeros([max(irrev2rev)+1,1])
    for i in range(len(rxn_exp_md_rev)):
        rxn_exp_md_rev[i] = sum(rxn_exp_md[rev2irrev[i]])
    #for i in range(len(irrev2rev)):
    #    rxn_exp_md_rev(irrev2rev(i)) = rxn_exp_md_rev(irrev2rev(i)) + rxn_exp_md(i)
    #cost_irrev = columnVector(f_easyLP).*columnVector(v_easyLP);
    cost_rev = np.zeros([max(irrev2rev)+1,1]);
    for i in range(len(irrev2rev)):
        cost_rev[irrev2rev[i]] = cost_rev[irrev2rev[i]] + cost_irrev[i];

    return [v_falcon, objValue, cost_rev, corr_rho, rxn_exp_md_rev]

if __name__=='__main__':
    whatthehell = cobra.io.read_sbml_model('/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/AGORAModels/Western-Diet-Paper/Yokenella_regensburgei_ATCC_43003.xml')
    #whatthehell = cobra.io.read_sbml_model('/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/AGORAModels/Western-Diet-Paper/Abiotrophia_defectiva_ATCC_49176.xml')
    #whatthehell.genes = []
    for i in range(len(whatthehell.reactions)):
        whatthehell.reactions[i].gene_reaction_rule = whatthehell.reactions[i].id
    selectgenenames = [rxn.id for rxn in whatthehell.reactions]
    unselectgenes = []
    for k in range(len(whatthehell.genes)):
        if whatthehell.genes[k].id not in selectgenenames:
            unselectgenes.append(whatthehell.genes[k])
    cobra.manipulation.delete.remove_genes(whatthehell,unselectgenes)
    expressionIDs = [rxn.id for rxn in whatthehell.reactions]
    expData = [1 for rxn in whatthehell.reactions]
    #expData = []
    #for i in range(len(whatthehell.reactions)):
    #    expData.append(np.random.random())
    expressionSDs = [1 for rxn in whatthehell.reactions]
    nReps = 1
    [v_falcon, objValue, cost_rev, corr_rho, rxn_exp_md_rev] = runFALCONStripped(whatthehell,expressionIDs,expData,expressionSDs,nReps)
    outFI = open('/mnt/vdb/home/ubuntu2/compareFALCONs.txt','w')
    for i in range(len(v_falcon)):
        outFI.write(str(v_falcon[i])+'\n')
    outFI.close()
