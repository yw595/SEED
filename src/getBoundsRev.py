import cobra
import numpy as np

def getBoundsRev(m):
    nrxns = len(m.reactions)
    nmets = len(m.metabolites)
    boundsRev = np.array([False for i in range(nrxns)])
    i=0
    while i < nrxns-1:
        if np.all(m.S[:,i]==-m.S[:,i+1]):
            if m.ub[i]>0 and m.ub[i+1]<0:
                boundsRev[i] = True
                boundsRev[i+1] = True
                i += 1
        i += 1
    return boundsRev
