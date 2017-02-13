import subprocess

metIn = open('/home/fs01/yw595/MATLAB/SEED/input/fecalmicrobiome.txt')
outFI = open('/home/fs01/yw595/MATLAB/SEED/input/fecalmicrobiomeauto.txt','w')
inLines = []
count = 0
for line in metIn:
    if count > 100:
        #break
        pass
    count = count+1
    line = line.strip()
    words = line.split('\t')
    outline = line
    if len(words)==1:
        proc = subprocess.Popen(['wget','www.hmdb.ca/metabolites/HMDB' + words[0],'-O','HMDB.html'])
        proc.wait()
        shFI = open('downloadMetabolome.sh','w')
        shFI.write('grep -Po "KEGG Compound ID(.)*C(\d)(\d)(\d)(\d)(\d)" HMDB.html > KEGG.txt')
        shFI.close()
        proc = subprocess.Popen(['chmod','u=rwx','downloadMetabolome.sh'])
        proc.wait()
        proc = subprocess.Popen(['sh','downloadMetabolome.sh'])
        proc.wait()
        keggIn = open('KEGG.txt')
        keggID = keggIn.readline().strip()
        keggID = keggID[len(keggID)-5:len(keggID)]
        keggIn.close()
        outline = words[0] + '\t' + keggID
    print(outline)
    outFI.write(outline + '\n')
outFI.close()
        #proc = subprocess.Popen(['grep','-Po','">C(\d)(\d)(\d)(\d)(\d)"','HMDB.html'])
        #proc.wait()
