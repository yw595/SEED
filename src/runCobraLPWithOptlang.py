from optlang import Model, Variable, Constraint, Objective
import time

def runCobraLPWithOptlang(A,b,c,anLb,aUb,osense,csense):
    xArr = []
    for i in range(len(anLb)):
        t1 = time.time()
        xArr.append(Variable('x'+str(i), lb=anLb[i], ub=aUb[i]))
        t2 = time.time()
        print(str(i)+'Var'+str(t2-t1))

    constraintArr = []
    for i in range(len(csense)):
        aconstraint = 0
        t1 = time.time()
        summa = 0
        for j in range(len(A[i,:])):
            if A[i,j]!=0:
                summa += 1
        if summa < 100:
            for j in range(len(A[i,:])):
                if A[i,j]!=0:
                    aconstraint += A[i,j]*xArr[j]
            if csense[i]=='L':
                constraintArr.append(Constraint(aconstraint, ub=b[i]))
            elif csense[i]=='G':
                constraintArr.append(Constraint(aconstraint, lb=b[i]))
            elif csense[i]=='E':
                constraintArr.append(Constraint(aconstraint, ub=b[i], lb=b[i]))
        t2 = time.time()
        print(str(i)+'Con'+str(t2-t1))

    anObj = 0
    for i in range(len(c)):
        t1 = time.time()
        if isinstance(c[i],(float,int)):
            if c[i]!=0:
                anObj += c[i]*xArr[i]
        else:
            if c[i][0]!=0:
                anObj += c[i][0]*xArr[i]
        t2 = time.time()
        print(str(i)+'Obj'+str(t2-t1))
    obj = Objective(anObj,direction='min')

    model = Model(name='Simple model')
    model.objective = obj
    model.add(constraintArr)

    model.optimize()
    stat = 0
    if model.status=='optimal':
        stat = 1
    return [stat, [var.primal for var in model.variables], model.objective.value]
