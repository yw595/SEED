from swiglpk import *
import scipy.sparse

def runCobraLPWithSWIGLPK(A,b,c,anLb,aUb,osense,csense):
    #ia = intArray((len(csense)+10)*(len(A[0,:])+10));
    #ja = intArray((len(csense)+10)*(len(A[0,:])+10));
    #ar = doubleArray((len(csense)+10)*(len(A[0,:])+10));
    ia = intArray(100000);
    ja = intArray(100000);
    ar = doubleArray(100000);
    lp = glp_create_prob();
    glp_set_prob_name(lp, "sample");
    glp_set_obj_dir(lp, GLP_MAX);
    glp_add_rows(lp, len(csense));
    for i in range(len(csense)):
        glp_set_row_name(lp, i+1, "cons"+str(i));
        bcoeff = b[i]
        if not isinstance(bcoeff,float):
            bcoeff = b[i][0]
        if csense[i]=='U':
            glp_set_row_bnds(lp, i+1, GLP_UP, 0.0, bcoeff);
        elif csense[i]=='L':
            glp_set_row_bnds(lp, i+1, GLP_LO, 0.0, bcoeff);
        elif csense[i]=='E':
            glp_set_row_bnds(lp, i+1, GLP_FX, 0.0, bcoeff);
    glp_add_cols(lp, len(aUb));
    for i in range(len(aUb)):
        glp_set_col_name(lp, i+1, "x"+str(i));
        if anLb[i]!=float('nan') and anLb[i]!=float('inf'):
            glp_set_col_bnds(lp, i+1, GLP_LO, 0.0, anLb[i]);
        if aUb[i]!=float('nan') and aUb[i]!=float('inf'):
            glp_set_col_bnds(lp, i+1, GLP_HI, 0.0, aUb[i]);
        cCoeff = c[i]
        if not isinstance(cCoeff,(float,int)):
            cCoeff = c[i][0]
        glp_set_obj_coef(lp, i+1, cCoeff);

    Asparse = scipy.sparse.dok_matrix(A)
    Akeys = Asparse.keys()
    flatidx = 0
    for i in range(len(Akeys)):
        flatidx += 1
        ia[flatidx] = Akeys[i][0]+1
        ja[flatidx] = Akeys[i][1]+1
        #nonsense = nonsense+1
        ar[flatidx] = A[Akeys[i][0],Akeys[i][1]]

    # flatidx = 0
    # for i in range(len(csense)):
    #     for j in range(len(A[i,:])):
    #         if A[i,j]!=0:
    #             flatidx += 1
    #             #flatidx = i*(len(A[i,:]))+(j+1)
    #             print(str(i)+' '+str(j))
    #             print(flatidx)
    #             ia[flatidx] = i+1
    #             ja[flatidx] = j+1
    #             ar[flatidx] = A[i,j]
    #glp_load_matrix(lp, len(csense)*len(A[0,:]), ia, ja, ar);
#    nonsense = nonsense+1
    glp_load_matrix(lp, flatidx, ia, ja, ar);
    glp_simplex(lp, None);
    Z = glp_get_obj_val(lp);
    xArr = []
    for i in range(len(A[0,:])):
        xArr.append(glp_get_col_prim(lp,i+1))
    glp_delete_prob(lp);
    nonsense += 1
