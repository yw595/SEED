import subprocess

#C:parse fecalmicrobiome.txt, first column is HMDB, second is KEGG, so if just one entry with HMDB ID, download html page for a HMDB ID from hmdb 
metIn = open('/home/fs01/yw595/MATLAB/SEED/input/fecalmicrobiome.txt')
outFI = open('/home/fs01/yw595/MATLAB/SEED/input/fecalmicrobiomeauto.txt','w')
inLines = []
count = 0
for line in metIn:
    count = count+1
    line = line.strip()
    words = line.split('\t')
    outline = line
    if len(words)==1:
        proc = subprocess.Popen(['wget','www.hmdb.ca/metabolites/HMDB' + words[0],'-O','HMDB.html'])
        proc.wait()

        #C:use downloadMetabolome.sh as grep wrapper to find whole line with KEGG ID, write line to KEGG.txt, get last five digits as KEGG ID, write out HMDB and KEGG IDs to fecalmicrobiomeauto.txt
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



