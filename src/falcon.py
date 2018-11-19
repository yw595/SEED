import time
import math
from getBoundsRev import getBoundsRev
import copy
from runCobraLPWithOptlang import runCobraLPWithOptlang
#from runCobraLPWithCvxopt import runCobraLPWithCvxopt
from runCobraLPWithGLPK import runCobraLPWithGLPK
from runCobraLPWithScipy import runCobraLPWithScipy
from runCobraLPWithSWIGLPK import runCobraLPWithSWIGLPK
from runCobraLPWithCobra import runCobraLPWithCobra
import sys
import numpy as np
import scipy.sparse as sp

def insertFunc(arr,index,insertval,nonsense=False):
    while index > len(arr)-1:
        if isinstance(arr,list):
            if isinstance(insertval,(int,long,float)):
                arr.append(0)
            elif isinstance(insertval,str):
                arr.append('')
        else:
            temp = len(arr)+1
            arr = np.resize(arr,temp)
            if isinstance(insertval,(int,long,float)):
                arr[temp-1] = 0
            elif isinstance(insertval,str):
                arr[temp-1] = ''
    if nonsense:
        nonsense2 = nonsense2+1
    arr[index]=insertval
    return arr

def clipper(inputArr):
    endidx = len(inputArr)-1
    while endidx>=0:
        if inputArr[endidx]!='ark':
            break
        else:
            endidx -= 1
    inputArr = inputArr[:endidx+1]
    for i in range(len(inputArr)):
        if inputArr[i]=='ark':
            if isinstance(inputArr[endidx],str):
                inputArr[i]=''
            else:
                inputArr[i]=0
    return inputArr
    
def isnanArr(arr,reverse=False):
    boolArr = []
    for i in range(len(arr)):
        if arr[i]==float('nan'):
            if reverse:
                boolArr.append(False)
            else:
                boolArr.append(True)
        else:
            if reverse:
                boolArr.append(True)
            else:
                boolArr.append(False)
    return(boolArr)

def falcon(m,r = [],r_sd = [],r_group = [],rc = 0,minFit = 0,LPmeth = 1,LPseed = float('inf'),flux_sum = -1,FDEBUG = False,EXPCON = False,TESTING = False,MASSPROD = False,ani='',aj=''):

    pFunctionName = 'falcon'
    pStructExpand = False
    pCaseSensitive = True
    pPartialMatching = False

    t_falcon = time.time()

    boundsRev = getBoundsRev(m)

    nrxns = len(m.reactions)
    nmets = len(m.metabolites)
    notnan_r = isnanArr(r,reverse=True)

    minUB = np.min(m.ub[m.ub > 0])
    if flux_sum < 0:
        flux_sum = sum(np.logical_and(np.logical_not(boundsRev),notnan_r))*minUB
        if flux_sum == 0:
            flux_sum = np.mean(m.ub[m.ub > 0])/2
    if FDEBUG:
        print(flux_sum)

    ZMIN = 0

    r_group = np.array(r_group)
    rgrp_notnan = r_group[notnan_r]
    vbasN = []
    cbasN = []
    vbasS = []
    cbasS = []
    nSnz = sum(sum(abs(np.sign(m.S))))
    ngroups = copy.deepcopy(r_group)
    if 1 not in ngroups:
        ngroups.append(1)
    ngroups = np.unique(ngroups)
    v_all = []

    ecrxns = []
    for i in range(len(m.rxnGeneMat)):
        if np.any(m.rxnGeneMat[i,:]!=0):
            ecrxns.append(True)
        else:
            ecrxns.append(False)

    r_sum = sum(r[notnan_r]) 
    r_pri_max = max(r)
    if r_sum!=0:
        r = flux_sum * r / r_sum
        r_sd = flux_sum * r_sd / r_sum
    if FDEBUG:
        print('Sum new r:')
        print(sum(r[notnan_r]))
        print('Max r (scaled r) is:')
        print([r_pri_max,max(r)])

    if TESTING:
        expZtol = 2*flux_sum/nrxns
    if FDEBUG:
        r_med = np.median(r)
        r_min = min(r[r>0])

    fOpt = 0
    v_orig = np.zeros([nrxns,1])
    m.lb[m.c == 1] = minFit

    rxnhasgene = ecrxns #rxnhasgene = (sum(m.rxnGeneMat')~=0)????

    corrval = float('nan')
    nR_old = 0
    v_sol = np.zeros([len(m.reactions),1])
    fIter = 0
    conv = 0
    nvar = float('nan')
    fUpdate = 0
    rGrpsUsed = 0
    flux_sum_pri = flux_sum
    v_pri = []
    print('FDEBUG')
    print(FDEBUG)
    while sum(np.logical_not(boundsRev)) > nR_old:
        fIter = fIter + 1
        nR_old = sum(np.logical_not(boundsRev))
        # nnnanrev = sum((notnan_r) & boundsRev) / 2
        nnan_irr = r_group[np.logical_and(notnan_r,np.logical_not(boundsRev))]
        nnnan_irr = len(np.unique(nnan_irr))

        r_group_cons = np.zeros([nrxns,1])
        if FDEBUG:
            NColLab = m.reactions
            NRowLab = m.metabolites

        #N = sp.lil_matrix([nmets+1+2*nrxns+1+2*nnnan_irr+1,nrxns+1+1+nnnan_irr,int(math.floor(2.3*nSnz))])
        N = np.zeros([nmets+1+2*nrxns+1+2*nnnan_irr+1,nrxns+1+1+nnnan_irr])
        dimFail = False
        if FDEBUG:
            sz_N = np.shape(N)
        sz_N = np.shape(N)

        N[0:nmets,0:nrxns] = m.S
        L = m.lb
        U = m.ub
        f = np.zeros([len(m.reactions),1])
        b = np.zeros([len(m.metabolites),1])
        csense = np.array(['' for blen in range(len(b))])
        if MASSPROD:
            csense[0:len(b)] = 'G'
        else:
            csense[0:len(b)] = 'E'

        s1 = nmets
        s2 = nrxns
        
        N[s1-1, s2] = 0 #n
        N[s1-1, s2 + 1] = 0 #z
        if FDEBUG:
            insertFunc(NColLab,s2,'n')
            insertFunc(NColLab,s2+1,'z')
        L = insertFunc(L,s2,0)
        U = insertFunc(U,s2,float('inf'))
        L = insertFunc(L,s2+1,ZMIN)
        U = insertFunc(U,s2+1,float('inf'))
        f = insertFunc(f,s2,0)
        f = insertFunc(f,s2+1,0)
        s2 = s2 + 2

        csense = insertFunc(csense,s1,'E')
        N[s1,nrxns] = 1
        if FDEBUG:
            insertFunc(NRowLab,s1,'LFP unitary')
        b = insertFunc(b,s1,1,nonsense=False)
        s1 = s1 + 1

        b=list(b)
        csense=list(csense)
        if FDEBUG:
            NRowLab=list(NRowLab)
            NColLab=list(NColLab)
        for k in range(nrxns):
            #print(k)
            f[k] = -rc #regularization constant 
            insertFunc(b,s1,0)
            insertFunc(b,s1+1,0)
            if not EXPCON or r(k)==float('nan') or U[k] == 0:
                N[s1,k] = 1 
                N[s1+1,k] = -1 
                N[s1,nrxns+1] = -U[k]
                N[s1+1,nrxns+1] = L[k]
            else:
                N[s1,k] = 1 
                N[s1+1,k] = 1 
                N[s1,nrxns] = -r[k]
                N[s1+1,nrxns] = -r[k]
            L[k] = 0
            U[k] = float('inf')
            if m.ub[k] == 0:
                U[k] = 0
            insertFunc(csense,s1,'L')
            insertFunc(csense,s1+1,'L') 
            if FDEBUG:
                insertFunc(NRowLab,s1,m.reactions[k].id+':U')
                insertFunc(NRowLab,s1+1,m.reactions[k].id+':L')
            s1+=2
            #print('HARSH')
        b=np.array(b)
        csense=np.array(csense)
        if FDEBUG:
            NRowLab=np.array(NRowLab)
            NColLab=np.array(NColLab)
        
        for k in range(len(ecrxns)):
            if ecrxns[k]:
                N[s1,k] = -1
        if np.ma.size(v_pri) > 0:
            flux_sum_pri = sum(v_pri[ecrxns])

        N[s1, nrxns + 1] = flux_sum
        b = insertFunc(b,s1,0) 
        csense = insertFunc(csense,s1,'L')
        if FDEBUG:
            insertFunc(NRowLab,s1,'FlxSum')
        s1+=1

        objPriorRow = s1
        csense = insertFunc(csense,objPriorRow,'L')
        insertFunc(b,objPriorRow,0)
        if FDEBUG:
            insertFunc(NRowLab,s1,'ObjPrior')
        s1+=1

        k = -1
        r_group_visited = np.array([False for nrxn in range(nrxns)])
        first_r_group_visited = -1
        rGrpsPrev = rGrpsUsed
        rGrpsUsed = 0
        benchMarkS1=s1
        b=list(b)
        b2 = copy.deepcopy(b)
        b.extend(['ark' for arkidx in range(len(b),sz_N[0],1)])
        f=list(f)
        f2 = copy.deepcopy(f)
        f.extend(['ark' for arkidx in range(len(f),sz_N[1],1)])
        L=list(L)
        L2 = copy.deepcopy(L)
        L.extend(['ark' for arkidx in range(len(L),sz_N[1],1)])
        U=list(U)
        U2 = copy.deepcopy(U)
        U.extend(['ark' for arkidx in range(len(U),sz_N[1],1)])
        csense=list(csense)
        csense2 = copy.deepcopy(csense)
        csense.extend(['ark' for arkidx in range(len(csense),sz_N[0],1)])
        if FDEBUG:
            NRowLab=list(NRowLab)
            NColLab=list(NColLab)
        #if sum(N[6153,0:4147]!=0)>3000:
        #    nonsense = nonsense+1
        #nonsense2 = nonsense2+1
        while k < nrxns-1:
            #print(k)
            time1 = time.time()
            k+=1
            s = r_sd[k]
            objDenom = s
            cons1 = 0
            if not r_group_visited[r_group[k]]:
                #print('HERE')
                cons1 = s1
                r_group_cons[r_group[k]] = cons1
            else:
                #print('THERE')
                cons1 = r_group_cons[r_group[k]]
            if isinstance(cons1,np.ndarray) and len(cons1)==1:
                cons1=int(cons1[0])
            if not boundsRev[k] and not r[k]==float('nan') and s > 0:
                if first_r_group_visited < 0:
                    first_r_group_visited = r_group[k]
                if not r_group_visited[r_group[k]]:
                    r_group_visited[r_group[k]] = True
                    s1+=2
                    rGrpsUsed+=1 
                    if r_group[k] != first_r_group_visited:
                        s2+=1
                #print(cons1)
                N[cons1,nrxns] = -r[k]
                N[cons1,k] = 1  
                N[cons1,s2] = -1
                b[cons1] = 0
                insertFunc(b2,cons1,0)
                N[cons1+1,nrxns] = r[k]
                N[cons1+1,k] = -1  
                N[cons1+1,s2] = -1
                b[cons1+1] = 0
                insertFunc(b2,cons1+1,0)
                f[s2] = -1/objDenom
                insertFunc(f2,s2,-1/objDenom)
                L[s2] = 0
                insertFunc(L2,s2,0)
                U[s2] = float('inf')
                insertFunc(U2,s2,float('inf'))
                csense[cons1] = 'L'
                csense[cons1+1] = 'L'
                insertFunc(csense2,cons1,'L')
                insertFunc(csense2,cons1+1,'L')
                fUpdate -= abs(v_orig[k] - r[k])/objDenom
                if FDEBUG:
                    insertFunc(NRowLab,cons1,'RG_'+str(r_group[k]))
                    insertFunc(NRowLab,cons1+1,'RG_'+str(r_group[k]))
                    insertFunc(NColLab,s2,'t_'+str(r_group[k]))
                if TESTING and fIter > 1 and abs(corrval) > 0:
                    N[objPriorRow,s2-1] = 1/objDenom

            time2 = time.time()
            #print(time2-time1)
        b=clipper(b)
        b=np.array(b)
        b=b[b!='ark']
        b2=np.array(b2)
        f=clipper(f)
        f=np.array(f)
        f=f[f!='ark']
        f2=np.array(f2)
        L=clipper(L)
        L=np.array(L)
        L=L[L!='ark']
        if len(L)==1:
            L=L[0]
        L2=np.array(L2)
        U=clipper(U)
        U=np.array(U)
        U=U[U!='ark']
        if len(U)==1:
            U=U[0]
        U2=np.array(U2)
        csense=clipper(csense)
        csense=np.array(csense)
        csense=csense[csense!='ark']
        csense2=np.array(csense)
        if FDEBUG:
            NRowLab=np.array(NRowLab)
            NColLab=np.array(NColLab)
        if TESTING and fIter > 1 and abs(corrval) > 0:
            N[objPriorRow,nrxns+1] = rGrpsUsed*(corrval/rGrpsPrev)

        for checkidx in range(len(r_group_visited)):
            if not r_group_visited[checkidx]:
                pass
                #nonsense = nonsense+1
    if not np.all(sz_N == np.shape(N)) or sz_N[0] != np.size(b) or sz_N[1] != np.size(L):
        print('WARNING: mismatch in estimated and actual dimension detected!!!')
        dimFail = True
        print(sz_N)
        sz_N_new = np.shape(N)
        blen = np.size(b)
        Llen = np.size(L)
        save('size_debug.mat', 'r_group', 'r_sd', 'r')

    if not dimFail:
        if True:
            if FDEBUG:
                print('Not Reversible: '+str(sum(not boundsRev)))
            #b[1406] = 1
            #for UTest in range(len(U)):
            #    U[UTest]=float('inf')
            #for senseTest in range(1407):
            #    csense[senseTest] = 'E'
            #for senseTest in range(1407,9699):
            #    csense[senseTest] = 'L'
            if len(b)==1 and len(b[0])>1:
                b = b[0]
            if len(f)==1 and len(f[0])>1:
                f = f[0]
            if len(L)==1 and len(L[0])>1:
                L = L[0]
            if len(U)==1 and len(U[0])>1:
                U = U[0]
            [v,fOpt,conv,vbasN,cbasN] = easyLP(f,N,b,L,U,csense,vbasN,cbasN,ani=ani,aj=aj)
            f_rxn = f[0:nrxns]
            v_rxn = v[0:nrxns]
            f_easyLP = f[nrxns+2:]
            v_easyLP = v[nrxns+2:]
            cost_easyLP = f_easyLP*v_easyLP
            cost_irrev = np.zeros([nrxns,1])
            r_group_cons_bycol = (r_group_cons-1-benchMarkS1)/2
            cost_naive = abs(r-abs(v_rxn))
            cost_naive[cost_naive==float('nan')]=0
            if True:
                cost_proxy=cost_naive
            else:
                cost_proxy=abs(v_rxn)
            cost_proxy2 = np.zeros([len(cost_proxy),1])
            for i in range(len(cost_proxy)):
                cost_proxy2[i] = cost_proxy[i]
            cost_proxy = cost_proxy2
            for i in range(len(f_easyLP)):
                cost_irrev[r_group_cons_bycol==i] = cost_easyLP[i]/sum(r_group_cons_bycol==i)
                if sum(cost_proxy[r_group_cons_bycol==i])!=0:
                    cost_irrev[r_group_cons_bycol==i] = np.inner(cost_irrev[r_group_cons_bycol==i],cost_proxy[r_group_cons_bycol==i])/sum(cost_proxy[r_group_cons_bycol==i])
                else:
                    cost_irrev[r_group_cons_bycol==i] = cost_irrev[r_group_cons_bycol==i]/sum(r_group_cons_bycol==i)

    if FDEBUG:
        print('fOpt, n, z:')
        print(str(fOpt)+' '+str(v[nrxns+1])+' '+str(v[nrxns+2]))
    if conv:
        v_pri = v
        v_orig = v
        if v[nrxns+1]!=0:
            v_orig = v/v[nrxns+1]
            corrval = fOpt/v[nrxns+1]
        v_sol = v_orig[0:nrxns]
        v_all.append(v[0:nrxns+2])
        nvar = v_orig[nrxns+1]
        [m.lb,m.ub,boundsRev] = setRxnDirection(v[0:nrxns],m.lb,m.ub,boundsRev,nrxns)
        if FDEBUG:
            print('New nvar, zvar is:')
            print(str(nvar)+' '+str(v[nrxns+1]))
            print('First 15 fluxes:')
            print(v_sol[0:15])
            print('Num Irrev, Previous Num Irrev:')
            print(str(sum(not boundsRev))+' '+str(nR_old))

    fTime = time.time() - t_falcon
    if conv:
        print('FALCON: solver converged in '+str(fTime)+' seconds and '+str(fIter)+' iterations.')
    else:
        print('FALCON: solver did NOT converge in '+str(fTime)+' seconds and '+str(fIter)+' iterations.')

    return [v_sol, corrval, nvar, v_all, fTime, fIter, fOpt, f_easyLP, v_easyLP, cost_irrev]


def easyLP(f,a,b,vlb,vub,csense,vbas,cbas,FDEBUG=False,ani='',aj=''):

    if np.any(f==float('nan')) or np.any(np.any(a==float('nan'))) or np.any(b==float('nan')) or np.any(vlb==float('nan')) or np.any(vub==float('nan')) or np.any(csense==float('nan')): 
        print('nan inputs not allowed')
        return

    print('EASY')
    v = np.zeros(np.shape(vlb))
    #v = v(:)
    #f = full(f(:))
    #vlb = vlb(:)
    #vub = vub(:)

    j1 = (vlb != vub)
    j2 = (vlb == vub)

    v[j2] = vlb[j2]
    b = b-np.dot(a,v)
    a = a[:,np.logical_not(j2)]
    vlb = vlb[np.logical_not(j2)] 
    vub = vub[np.logical_not(j2)]
    f0 = f
    f = f[np.logical_not(j2)]
    fOpt = float('nan')

    if np.any(b==float('nan')):
        print('nan inputs not allowed: something went wrong')
        return

    LPmeth = None
    LPseed = None
    paramsmethod = LPmeth
    if LPseed==float('inf'):
        paramsseed = np.random.randint(low=math.floor(sys.maxint/10))
    else: 
        paramsseed = LPseed
    if len(vbas) > 0:
        paramsvbasis = vbas
        paramscbasis = np.zeros([1, len(b)])
        paramscbasis[0:len(cbas)] = cbas

    if FDEBUG:
        t_easy = time.time()
    #[solutionstat,solutionfull,solutionobj] = runCobraLPWithOptlang(a,b,f,vlb,vub,-1,csense)
    #[solutionstat,solutionfull,solutionobj] = runCobraLPWithCvxopt(a,b,f,vlb,vub,-1,csense)
    #[solutionstat,solutionfull,solutionobj] = runCobraLPWithGLPK(a,b,f,vlb,vub,-1,csense)
    #[solutionstat,solutionfull,solutionobj] = runCobraLPWithScipy(a,b,f,vlb,vub,-1,csense)
    #nonsense = nonsense+1
    #nonsensepeasy = nonsensepeasy+1
    [solutionstat,solutionfull,solutionobj] = runCobraLPWithGLPK(a,b,f,vlb,vub,-1,csense,ani=ani,aj=aj)
    #[solutionstat,solutionfull,solutionobj] = runCobraLPWithCobra(a,b,f,vlb,vub,-1,csense,ani=ani,aj=aj)
    if FDEBUG:
        t_easy = time.time()-t_easy

    conv = solutionstat == 1
    svbas = []
    scbas = []

    if conv:
        v0 = solutionfull
        v[j1] = v0
        fOpt = np.inner(f0,v)
        if np.inner(f,v0)!=0:
            #print(f)
            #print(v0)
            #print(np.inner(f,v0))
            pass
        if FDEBUG:
            print('Convergent optimum is: '+str(solutionobj))
            if fOpt==float('nan'):
                print('Converged, but fOpt still nan!')
            NColLab = NColLab[np.logical_not(j2)]
    else:
        solstat = solutionstat
        sol = solutionfull   
    return [v, fOpt, conv, svbas, scbas]

def setRxnDirection(vI, iLB, iUB, isRev, nrxns):
    rthresh = 0.5

    tol = 1e-9

    k = -1
    while k < nrxns-1:
        k+=1
        if isRev[k]:
            vSum = vI[k] + vI[k+1]
            if vI[k]/vSum > rthresh:
                iLB[k+1] = 0
                iUB[k+1] = 0

                isRev[k] = False
                isRev[k+1] = False
            elif vI[k+1]/vSum > rthresh:
                iLB[k] = 0
                iUB[k] = 0

                isRev[k] = False
                isRev[k+1] = False      
            k+=1
    return [iLB, iUB, isRev]
    
def setFBRxnDirection(vI, iLB, iUB, isRev, nrxns, m):

    rthresh = 0.5
    tol = 1e-8
    k = 0
    while k < nrxns-1:
        k+=1
        if isRev[k]:
            vSum = vI[k] + vI[k+1]
            if vI[k] >= tol and vI[k] - vI[k+1] <= tol:
                if FDEBUG:
                    print('set rxns to zero:')
                    print(m.reactions[k:k+2])
                    print(vI[k:k+2])
                    print('END set rxns to zero')
                iLB[k+1] = 0
                iUB[k+1] = 0
                isRev[k] = 0
                isRev[k+1] = 0
                iLB[k] = 0
                iUB[k] = 0
                isRev[k] = 0
                isRev[k+1] = 0      
            k+=1
    return [iLB, iUB, isRev]

def countNonZeroEq(vI, isRev, nrxns):
    tol = 1e-8
    k = 0
    nNZE = 0
    while k < nrxns: 
        k+=1
        if isRev[k]:
            if vI[k] >= tol and vI[k] - vI[k+1] <= tol:
                nNZE+=1
            k+=1
    return nNZE
