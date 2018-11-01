import time
import numpy as np
import scipy.sparse
import subprocess
from cobra import Model, Reaction, Metabolite

def runCobraLPWithCobra(A,b,c,anLb,aUb,osense,csense,ani='',aj=''):

    model = Model('example_model')
    reactions = []
    for i in range(len(c)):
        reaction = Reaction('r'+str(i))
        reaction.name = ''
        reaction.subsystem = ''
        if anLb[i]!=float('nan') and anLb[i]!=float('inf'):
            reaction.lower_bound = anLb[i]
        if aUb[i]!=float('nan') and aUb[i]!=float('inf'):
            reaction.upper_bound = aUb[i]
        reactions.append(reaction)
    model.add_reactions(reactions)

    Asparse = scipy.sparse.dok_matrix(A)
    Akeys = Asparse.keys()
    ithToX = {}
    for i in range(len(Akeys)):
        if Akeys[i][0] not in ithToX:
            ithToX[Akeys[i][0]] = [[],[]]
        ithToX[Akeys[i][0]][0].append('x'+str(Akeys[i][1]))
        ithToX[Akeys[i][0]][1].append(A[Akeys[i][0],Akeys[i][1]])
    for i in range(7000):#len(csense)):
        if i not in ithToX:
            continue
        xNames = ithToX[i][0]
        xCoeffs = ithToX[i][1]
        coefficients = dict()
        print(i)
        for j in range(len(xNames)):
            coeff = xCoeffs[j]
            xName = xNames[j]
            if coeff!=0:
                coefficients[reactions[int(xName[1:])].forward_variable] = float(coeff)
                coefficients[reactions[int(xName[1:])].reverse_variable] = -float(coeff)
        bVar = 0
        if isinstance(b[i],float):
            bVar = b[i]
        else:
            bVar = b[i][0]
        if csense[i]=='L':
            constraint = model.problem.Constraint(0, lb=bVar)
        elif csense[i]=='U':
            constraint = model.problem.Constraint(0, ub=bVar)
        elif csense[i]=='E':
            constraint = model.problem.Constraint(0, lb=bVar, ub=bVar)
        model.add_cons_vars(constraint)
        model.solver.update()
        constraint.set_linear_coefficients(coefficients=coefficients)

    coefficients = dict()
    for i in range(len(c)):
        if isinstance(c[i],(float,int)):
            coeff = c[i]
        else:
            coeff = c[i][0]
        coefficients[reactions[i].forward_variable] = float(coeff)
        coefficients[reactions[i].reverse_variable] = -float(coeff)

    test_objective = model.problem.Objective(0,direction='max')
    model.objective = test_objective
    test_objective.set_linear_coefficients(coefficients=coefficients)
    solution = model.optimize(objective_sense=None)
    nonsense = nonsense+1
