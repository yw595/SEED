from ete3 import NCBITaxa
import itertools

SEEDSpecies = []
with open('/home/ubuntu/SEEDSpecies.txt') as f:
    for line in f:
        SEEDSpecies.append(line.strip())
    f.close()

HMPFamilies = []
with open('/home/ubuntu/HMPFamilies.txt') as f:
    for line in f:
        SEEDSpecies.append(line.strip())
    f.close()

ncbi = NCBITaxa()
name2taxid = ncbi.get_name_translator(list(set(SEEDSpecies+HMPFamilies)))

tree = ncbi.get_topology(list(itertools.chain.from_iterable(list(name2taxid.values()))),intermediate_nodes=True)
#print(tree.get_ascii(attributes=['sci_name']), file=open('/home/ubuntu/taxonomy.txt','w'))

#print(tree.name)
#fh = open('/home/ubuntu/closestFamilies.txt','w')
for species in SEEDSpecies:
    print(species)
    minDist = -1
    minDistSpecies = ''
    for family in HMPFamilies:
        familyNode = tree.search_nodes(name=str(name2taxid[family][0]))[0]
        speciesNode = tree.search_nodes(name=str(name2taxid[species][0]))[0]
        dist = tree.get_distance(speciesNode, familyNode)
        if minDist == -1:
            minDist = dist
            minDistFamily = family
        elif minDist < dist:
            minDist = dist
            minDistFamily = family
    print(species+" "+minDistFamily+" "+str(minDist))
#     print(genus+"\t"+minDistSpecies,file=fh)
# fh.close()



























