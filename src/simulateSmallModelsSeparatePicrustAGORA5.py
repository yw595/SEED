from optparse import OptionParser
from biom import load_table
import subprocess
import sys
import simulateSmallModelsSeparatePicrustAGORA4 as sa

usage = "usage %prog [options] \n"
prog = "this prog"
parser = OptionParser(usage=usage)
parser.add_option("--ucrFolder",help="")
(opts,args) = parser.parse_args()
ucrFolder = opts.ucrFolder

sa.makeSepBiom(ucrFolder,[True,True,True],0,False)
sa.makeSepExprAndEC('',ucrFolder,['','',''],0,False,False,True)
