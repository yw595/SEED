import numpy as np
import time
from falcon import falcon

def falconMulti(m, nReps, rxn_exp_md=[], rxn_exp_sd_md=[], rxn_rule_group=[], rc=0.0,minFit=0.0,EXPCON=False,FDEBUG=False,ani='',aj=''):

    nRxns = len(m.reactions)
    v_sol_Dist   = np.zeros([nReps, nRxns])
    corrval_Dist = np.zeros([nReps, 1])
    nvar_Dist    = np.zeros([nReps, 1])
    fTime_Dist   = np.zeros([nReps, 1])
    fIter_Dist   = np.zeros([nReps, 1])
    v_all = []

    timeInit = str(time.time())

    for i in range(nReps):
        print(nReps)
        [v_sol, corrval, nvar, temp1, fTime, fIter, fOpt, f_easyLP, v_easyLP, cost_irrev] = falcon(m, r=rxn_exp_md, r_sd=rxn_exp_sd_md, r_group=rxn_rule_group,rc=rc,minFit=minFit,EXPCON=EXPCON,FDEBUG=FDEBUG,ani=ani,aj=aj)
        v_sol_Dist[i,:] = np.reshape(v_sol,[1,nRxns])
        corrval_Dist[i] = corrval
        fTime_Dist[i] = fTime
        fIter_Dist[i] = fIter

    v_sol   = np.mean(v_sol_Dist,axis=0)
    nvar    = np.mean(nvar_Dist)
    corrval = np.mean(corrval_Dist)
    fTime   = np.mean(fTime_Dist)
    fIter   = np.mean(fIter_Dist)

    v_sol_s   = np.std(v_sol_Dist,axis=0)
    nvar_s    = np.std(nvar_Dist)
    corrval_s = np.std(corrval_Dist)
    fTime_s   = np.std(fTime_Dist)
    fIter_s   = np.std(fIter_Dist)

    return [v_sol, corrval, nvar, v_all, fTime, fIter, v_sol_s, corrval_s, nvar_s, fTime_s, fIter_s, fOpt, f_easyLP, v_easyLP, cost_irrev]
