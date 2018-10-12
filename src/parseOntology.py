import re

inputDir = '/home/fs01/yw595/MATLAB/SEED/input/MGMData'

#C:parse ko.id2hierachy, make map KOsToECs of all ECs for a KO
p = re.compile('[0-9|-]+\.[0-9|-]+\.[0-9|-]+\.[0-9|-]+')
KOFI = open(inputDir + '/KEGG/ko.id2hierachy')
KOsToECs = {}
for line in KOFI:
    words = line.strip().split('\t')
    m = p.findall(words[4])
    if m:
        KOsToECs[words[5]] = list(set(m))
KOFI.close()

#C:parse ontology_map, for K KEGG entries, make map ontoIDsToKOs of all KOs for an ontoID
ontologyFI = open(inputDir + '/ontology_map')
ontoIDsToKOs = {}
for line in ontologyFI:
    words = line.strip().split('\t')
    if words[1].startswith('K'):
        ontoIDsToKOs[words[0]] = words[1]
ontologyFI.close()

#C:parse md5_ontology_map, for type 12 entries, make map md5sToOntoIDs of one ontoID for each md5
md5FI = open(inputDir + '/md5_ontology_map')
md5sToOntoIDs = {}
for line in md5FI:
    words = line.strip().split('\t')
    if int(words[1])==12:
        md5sToOntoIDs[words[0]] = words[3]
md5FI.close()

#C:using previous maps, make map md5ToECOnto of all ECs for an md5, using JUST FIRST KO FROM ontoIDsToKOs???, use counting map numKOs, print diagnostic of every 10000 KOsToECs entry
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

#C:write tab-delimited md5ToECOnto.map, each md5 has JUST ONE EC ASSOCIATED?
mapFI = open(inputDir + '/md5ToECOnto.map','w')
for md5 in md5ToECOnto:
    mapFI.write(md5 + '\t' + str(md5ToECOnto[md5]) + '\n')
mapFI.close()
