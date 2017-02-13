import re

inputDir = '/home/fs01/yw595/MATLAB/SEED/input/MGMData'

p = re.compile('[0-9|-]+\.[0-9|-]+\.[0-9|-]+\.[0-9|-]+')
KOFI = open(inputDir + '/KEGG/ko.id2hierachy')
KOsToECs = {}
for line in KOFI:
    words = line.strip().split('\t')
    m = p.findall(words[4])
    if m:
        KOsToECs[words[5]] = list(set(m))
KOFI.close()

ontologyFI = open(inputDir + '/ontology_map')
ontoIDsToKOs = {}
for line in ontologyFI:
    words = line.strip().split('\t')
    if words[1].startswith('K'):
        ontoIDsToKOs[words[0]] = words[1]
ontologyFI.close()

md5FI = open(inputDir + '/md5_ontology_map')
md5sToOntoIDs = {}
for line in md5FI:
    words = line.strip().split('\t')
    if int(words[1])==12:
        md5sToOntoIDs[words[0]] = words[3]
md5FI.close()

count = 0
numKOs = {}
md5ToECOnto = {}
for md5 in md5sToOntoIDs:
    KO = ontoIDsToKOs[md5sToOntoIDs[md5]]
    if KO in KOsToECs:
        numKOs[KO] = ''
        count = count+1
        md5ToECOnto[md5] = KOsToECs[KO]
        if (count % 10000)==0:
            print(KOsToECs[KO])
            print(count)
print(len(numKOs.keys()))

mapFI = open(inputDir + '/md5ToECOnto.map','w')
for md5 in md5ToECOnto:
    mapFI.write(md5 + '\t' + str(md5ToECOnto[md5]) + '\n')
mapFI.close()
