import time
import numpy as np
from bootstrapAGORAIBD import writeData
import subprocess

def computeMinDisj(model, genedata_filename, sigma=0,FDEBUG=True,overwriteSDs={},overwriteMeans={}):
    ztol = 1e-4

    mdesc = model.id.replace(' ','_')
    rfid = str(np.random.randint(low=1,high=10e10))
    rfname = genedata_filename+'_'+mdesc+rfid
    rfname = rfname.replace(' ','')
    model.grRules = [rxn.gene_reaction_rule for rxn in model.reactions]
    writeData([model.grRules],rfname,delimiter=',')
    rfout = rfname+'_out'
    nrxns = len(model.reactions)
    if FDEBUG:
        print('Filebase: ')
        print(rfname)

    print('minDisj HERE')
    #print(rfout)
    #print(genedata_filename)
    rfoutFI = open(rfout,'w')
    proc = subprocess.Popen(['/mnt/vdb/home/ubuntu2/minDisj',genedata_filename,rfname],stdout=rfoutFI)
    status = proc.wait()
    rfoutFI.close()
    print('minDisj THERE')
    if status != 0:
        time.sleep(0.03)
        print('try #2...')
        rfoutFI = open(rfout,'w')
        proc = subprocess.Popen(['/mnt/vdb/home/ubuntu2/minDisj',genedata_filename,rfname],stdout=rfoutFI)
        status = proc.wait()
        rfoutFI.close()
        if status != 0:
            time.sleep(3)
            print('try #3...')
            rfoutFI = open(rfout,'w')
            proc = subprocess.Popen(['/mnt/vdb/home/ubuntu2/minDisj',genedata_filename,rfname],stdout=rfoutFI)
            status = proc.wait()
            rfoutFI.close()
            if status != 0: 
                minDisjCmd = '/mnt/vdb/home/ubuntu2/minDisj '+genedata_filename+' '+rfname+' > '+rfout
                print('minDisj failed with return code '+str(status))
                return

    if not FDEBUG:
        proc = subprocess.Popen(['rm',rfname])
        proc.wait()
    if sigma > 0 and not FDEBUG:
        proc = subprocess.Popen(['rm',genedata_filename])
        proc.wait()

    ruleFirstIdx = {}
    cidx = -1
    rxn_rule_group = [0 for nrxn in range(nrxns)]
    #print(len(rxn_rule_group))
    #print(len(model.grRules))
    while cidx < nrxns-1:
        cidx = cidx + 1
        key = model.grRules[cidx].strip()
        if len(key) > 0:
            if key in ruleFirstIdx:
                rxn_rule_group[cidx] = ruleFirstIdx[key]
            else:
                rxn_rule_group[cidx] = cidx
                ruleFirstIdx[key] = cidx
        else:
            rxn_rule_group[cidx] = cidx

    tempFI = open(rfout)
    print(rfout)
    rxnData = []
    for line in tempFI:
        words = line.strip().split()
        rxnData.append([float(words[0]), float(words[1])])
    tempFI.close()
    if not FDEBUG:
        proc = subprocess.Popen(['rm',rfout])
        proc.wait()
    print(nrxns)
    if len(rxnData) == nrxns-1:
        rxnData.append([float('nan'),float('nan')])

    rxnData = np.array(rxnData)
    rxn_exp = rxnData[:,0]
    rxn_exp_sd = rxnData[:,1]
    if FDEBUG:
        print('min max (rxn_exp then rxn_exp_sd)')
        print([min(rxn_exp), max(rxn_exp)])
        print([min(rxn_exp_sd), max(rxn_exp_sd)])
    if max(rxn_exp_sd) > ztol:
        rxn_exp_sd[np.logical_or(rxn_exp_sd < ztol,rxn_exp_sd==float('nan'))] = min(rxn_exp_sd[rxn_exp_sd >= ztol])/2
    else:
        rxn_exp_sd[np.logical_or(rxn_exp_sd < ztol,rxn_exp_sd==float('nan'))] = 1
    if FDEBUG:
        print('New min rxn_exp_sd')
        print(min(rxn_exp_sd))

    return [rxn_exp,rxn_exp_sd,rxn_rule_group]
