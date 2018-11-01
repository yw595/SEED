import time
import numpy as np
import scipy.sparse
import subprocess
from testpython import *

def runCobraLPWithGLPK(A,b,c,anLb,aUb,osense,csense,ani='',aj=''):

    LPfile = '/mnt/vdb/home/ubuntu2/glpk.txt'
    if ani!='' and aj!='':
        LPfile = '/mnt/vdb/home/ubuntu2/glpk_'+str(ani)+'_'+str(aj)+'.txt'
    LPFI = open(LPfile,'w')
    LPFI.write('Maximize\n')
    objLine = '   obj:'
    for i in range(len(c)):
        if isinstance(c[i],(float,int)):
            coeff = c[i]
        else:
            coeff = c[i][0]
        if coeff>=0:
            if i==0:
                objLine += ' '
            else:
                objLine += ' + '
        else:
            objLine += ' - '
        objLine += str(abs(coeff))+' '+'x'+str(i)
    LPFI.write(objLine+'\n')
    
    LPFI.write('Subject To\n')
    Asparse = scipy.sparse.dok_matrix(A)
    Akeys = Asparse.keys()
    ithToX = {}
    for i in range(len(Akeys)):
        if Akeys[i][0] not in ithToX:
            ithToX[Akeys[i][0]] = [[],[]]
        ithToX[Akeys[i][0]][0].append('x'+str(Akeys[i][1]))
        ithToX[Akeys[i][0]][1].append(A[Akeys[i][0],Akeys[i][1]])
    for i in range(len(csense)):
        consLine = '   c'+str(i)+':'
        if i not in ithToX:
            #print('CONTINUE')
            #print(i)
            continue
        xNames = ithToX[i][0]
        xCoeffs = ithToX[i][1]
        #xNames = np.array(['x'+str(x) for x in range(len(A[0,:]))])
        #xNames = xNames[A[i,:]!=0]
        #xCoeffs = A[i,A[i,:]!=0]
        #print(len(xNames))
        #print(len(xCoeffs))
        #print(i)
        for j in range(len(xNames)):
            coeff = xCoeffs[j]
            xName = xNames[j]
            if coeff!=0:
                if coeff>0:
                    if j==0:
                        consLine += ' '
                    else:
                        consLine += ' + '
                else:
                    consLine += ' - '
                consLine += str(abs(coeff))+' '+xName
        # for j in range(len(A[0,:])):
        #     coeff = A[i,j]
        #     if coeff!=0:
        #         if coeff>0:
        #             if j==0:
        #                 consLine += ' '
        #             else:
        #                 consLine += ' + '
        #         else:
        #             consLine += ' - '
        #         consLine += str(abs(coeff))+' '+'x'+str(j)
        if csense[i]=='L':
            consLine += ' <= '
        elif csense[i]=='G':
            consLine += ' >= '
        elif csense[i]=='E':
            consLine += ' = '
        if isinstance(b[i],float):
            consLine += str(b[i])
        else:
            consLine += str(b[i][0])
        if csense[i]=='E' or csense[i]=='U' or csense[i]=='L':
            LPFI.write(consLine+'\n')

    LPFI.write('Bounds\n')
    for i in range(len(A[0,:])):
        boundsLine = '   '
        setLB = False
        if anLb[i]!=float('nan') and anLb[i]!=float('inf'):
            boundsLine += str(anLb[i])+' <= '+'x'+str(i)
            setLB = True
        if aUb[i]!=float('nan') and aUb[i]!=float('inf'):
            if not setLB:
                boundsLine += 'x'+str(i)
            boundsLine += ' <= '+str(aUb[i])
        LPFI.write(boundsLine+'\n')
    LPFI.write('End\n')

    LPFI.close()

    print('etsi')
    outfile = '/mnt/vdb/home/ubuntu2/glpkout.txt'
    if ani!='' and aj!='':
        outfile = '/mnt/vdb/home/ubuntu2/glpkout_'+str(ani)+'_'+str(aj)+'.txt'
    myFunction(LPfile,outfile)
    #proc = subprocess.Popen(['glpsol','--lp',LPfile,'--tmlim','100','-o',outfile])
    #proc.wait()
    xCount = 0
    xFlux = []
    status = 0
    #return [1, [0 for cidx in range(len(c))], 0]
    inFI = open(outfile)
    for line in inFI:
        words = line.strip().split()
        #print(words)
        if len(words)>=2 and words[0]=='Status:':
            if words[1]=='OPTIMAL':
                status = 1
        if len(words)>=2 and words[0]=='Objective:':
            objval = float(words[3])
        if len(words)>=5 and words[1]==('x'+str(xCount)):
            xFlux.append(float(words[3]))
            xCount += 1
    inFI.close()

    #proc = subprocess.Popen(['rm',LPfile])
    #proc.wait()
    #proc = subprocess.Popen(['rm',outfile])
    #proc.wait()
    return [status, xFlux, objval]
