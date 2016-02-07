from Bio.Seq import Seq
from Bio.KEGG import Enzyme
from Bio.KEGG import REST
from Bio.KEGG.KGML import KGML_parser
from Bio.KEGG import Map
#request = REST.kegg_get("ec:5.4.2.2")
#open("ec_5.4.2.2.txt",'w').write(request.read())
#records = Enzyme.parse(open("ec_5.4.2.2.txt"))
#record = list(records)[0]
#print(record.classname)
#print(record.entry)
organisms = REST.kegg_list("organism").read()
organismlist = []
for line in organisms.rstrip().split("\n"):
    #print(line)
    code = line.split("\t")[1]
    organismlist.append(code)

#print(organismlist)

#parser = KGML_parser.KGMLparser()
#open("human_map.xml",'w').write(REST.kegg_get("hsa05130",option="kgml").read())
human_map = KGML_parser.read(REST.kegg_get("hsa01100",option="kgml"))
cpds = human_map.compounds
for cpd in cpds:
    print(cpd.name)
    graphics = cpd.graphics
    for graphic in graphics:
        print(graphic.x)

rxns = human_map.reaction_entries
for rxn in rxns:
    print(rxn.name)
    graphics = rxn.graphics
    for graphic in graphics:
        print(graphic.x)
