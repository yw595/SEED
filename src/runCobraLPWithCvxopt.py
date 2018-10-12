from cvxopt import matrix, solvers
import cvxopt
import time

def runCobraLPWithCvxopt(A,b,c,anLb,aUb,osense,csense):

    Atemp = []
    for j in range(len(A[0,:])):
        Atemp.append([])
    btemp = []
    for i in range(len(csense)):
        for j in range(len(A[0,:])):
            if csense[i]=='L':
                Atemp[j].append(A[i,j])
            elif csense[i]=='G':
                Atemp[j].append(-A[i,j])
            elif csense[i]=='E':
                Atemp[j].append(A[i,j])
                Atemp[j].append(-A[i,j])
        if csense[i]=='E':
            if isinstance(b[i],list):
                btemp.append(b[i][0])
                btemp.append(b[i][0])
            else:
                btemp.append(b[i])
                btemp.append(b[i])
        elif csense[i]=='L' or csense[i]=='U':
            if isinstance(b[i],list):
                btemp.append(b[i][0])
            else:
                btemp.append(b[i])

    for i in range(len(anLb)):
        for j in range(len(A[0,:])):
            if i==j:
                Atemp[j].append(1)
                Atemp[j].append(-1)
            else:
                Atemp[j].append(0)
                Atemp[j].append(0)
        btemp.append(aUb[i])
        btemp.append(anLb[i])
                
    Atemp = matrix(Atemp)
    btemp2 = []
    for i in range(len(btemp)):
        if isinstance(btemp[i],(float,int)):
            btemp2.append(btemp[i])
        else:
            btemp2.append(btemp[i][0])
    btemp = matrix(btemp2)
    ctemp = []
    for i in range(len(c)):
        t1 = time.time()
        if isinstance(c[i],float):
            ctemp.append(c[i])
        else:
            ctemp.append(c[i][0])
    ctemp = matrix(ctemp)
    #sol = solvers.lp(ctemp,Atemp,btemp)
    sol = cvxopt.glpk.lp(ctemp,Atemp,btemp)
    stat = 0
    if sol['status']=='optimal':
        stat = 1
    return [stat, [sol['x'][i] for i in range(len(sol['x']))], sol['primal objective']]
