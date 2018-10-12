import sys
import subprocess
from biom import load_table
from simulateSmallModelsSeparatePicrustAGORA import makeSepBiom
from simulateSmallModelsSeparatePicrustAGORA import makeSepExprAndEC

ucrFolder = sys.argv[1]
#makeSepBiom(ucrFolder,[],800)
makeSepExprAndEC('',ucrFolder,[],800,False)
