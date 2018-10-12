import time
import numpy as np
import scipy.optimize

def runCobraLPWithScipy(A,b,c,anLb,aUb,osense,csense):

    ctemp = [c[i] for i in range(len(c))]

    Aeq = []
    beq = []
    Aub = []
    bub = []
    for i in range(len(csense)):
        print(i)
        if isinstance(b[i],float):
            bcoeff = b[i]
        else:
            bcoeff = b[i][0]
        if csense[i]=='L':
            Aub.append(list(A[i,:]))
            bub.append(bcoeff)
        elif csense[i]=='G':
            Aub.append(list(A[i,:]))
            bub.append(-bcoeff)
        elif csense[i]=='E':
            Aeq.append(list(A[i,:]))
            beq.append(bcoeff)

    abounds = []
    for i in range(len(aUb)):
        print(i)
        ithLb = None
        if anLb[i]!=float('nan') and anLb[i]!=float('inf'):
            ithLb = anLb[i]
        ithUb = None
        if aUb[i]!=float('nan') and aUb[i]!=float('inf'):
            ithUb = aUb[i]
        abounds.append((ithLb,ithUb))

    res = scipy.optimize.linprog(ctemp, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=abounds,options={"disp": True})
    nonsense = nonsense+1
