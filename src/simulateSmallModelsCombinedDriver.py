import subprocess
import time
from optparse import OptionParser

diabetesctrl = ['N044A','SZEY-75A'] #from diabetesctrl.txt
for sampleid in diabetesctrl:
    for i in range(10):
        proc = subprocess.Popen(['python','/mnt/vdb/home/ubuntu2/MATLAB/SEED/src/simulateSmallModelsCombined.py','--diseaseType','diabetes','--sampleid',sampleid])
        time.sleep(900)
        proc.kill()

diabetestype2nomet = ['NG-5636_551','DOM024','DLM001'] #from diabetestype2nomet.txt
for sampleid in diabetestype2nomet:
    for i in range(10):
        proc = subprocess.Popen(['python','/mnt/vdb/home/ubuntu2/MATLAB/SEED/src/simulateSmallModelsCombined.py','--diseaseType','diabetes','--sampleid',sampleid])
        time.sleep(900)
        proc.kill()

HMP2Normal = ['206703','206704','206700'] #from HMP2Normal.txt
for sampleid in HMP2Normal:
    for i in range(10):
        proc = subprocess.Popen(['python','/mnt/vdb/home/ubuntu2/MATLAB/SEED/src/simulateSmallModelsCombined.py','--diseaseType','IBD','--sampleid',sampleid])
        time.sleep(900)
        proc.kill()

HMP2IBD = ['206701','206708','206709'] #from HMP2IBD.txt
for sampleid in HMP2IBD:
    for i in range(10):
        proc = subprocess.Popen(['python','/mnt/vdb/home/ubuntu2/MATLAB/SEED/src/simulateSmallModelsCombined.py','--diseaseType','IBD','--sampleid',sampleid])
        time.sleep(900)
        proc.kill()

# MHnormal = ['MH0005','MH0006','MH0008'] #from MHnormal.txt
# for sampleid in MHnormal:
#     for i in range(10):
#         proc = subprocess.Popen(['python','/mnt/vdb/home/ubuntu2/MATLAB/SEED/src/simulateSmallModelsCombined.py','--diseaseType','normObese','--sampleid',sampleid])
#         time.sleep(900)
#         proc.kill()

# MHobese = ['MH0001','MH0002','MH0003'] #from MHobese.txt
# for sampleid in MHobese:
#     for i in range(10):
#         proc = subprocess.Popen(['python','/mnt/vdb/home/ubuntu2/MATLAB/SEED/src/simulateSmallModelsCombined.py','--diseaseType','normObese','--sampleid',sampleid])
#         time.sleep(900)
#         proc.kill()
