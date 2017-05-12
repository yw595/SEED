import subprocess

#proc = subprocess.Popen(['wget','ftp://public-ftp.hmpdacc.org/HMQCP/otu_table_v13.txt.gz'])
#proc.wait()
#proc = subprocess.Popen(['gunzip','otu_table_v13.txt.gz'])
#proc.wait()

inputFI = open('otu_table_v13.txt')
lineNo = 0
allEnds = []
fh = open('HMPFamilies.txt','w')
fh2 = open('HMPFamiliesAbundsAndOccurs.txt','w')
fh3 = open('HMPAllAbunds.txt','w')
familiesToAbunds = {}
familiesToOccurs = {}
familiesToAllAbunds = {}
for line in inputFI:
    lineNo = lineNo+1
    print(lineNo)
    if lineNo > 3000:
        #break
        pass
    words = line.split('\t')
    family = words[len(words)-1].split(';')
    family2 = family[len(family)-1]
    #allEnds.append(family2)
    #allEnds = list(set(allEnds))
    if family2[0:3]=='f__':
        allEnds.append(family2[3:len(family2)])
        allEnds = list(set(allEnds))
    if lineNo >= 3:#len(words) > 1000:
        if family2[0:3]=='f__':
            family3 = family2[3:len(family2)]
        firstTime = False
        for i in range(1,len(words)-1):
            #print(words[1000])
            #print(family2[3:len(family2)])
            #print(family2[0:3])
            if family3 not in familiesToAbunds:
                familiesToAbunds[family3] = int(words[i])
                familiesToOccurs[family3] = 0
                if int(words[i])!=0:
                    familiesToOccurs[family3] = familiesToOccurs[family3] + 1
                familiesToAllAbunds[family3] = [words[i]]
                firstTime = True
            else:
                familiesToAbunds[family3] = familiesToAbunds[family3] + int(words[i])
                if int(words[i])!=0:
                    familiesToOccurs[family3] = familiesToOccurs[family3] + 1
                if firstTime:
                    familiesToAllAbunds[family3].append(words[i])     
                else:
                    #print(lineNo)
                    #print(familiesToAllAbunds[family3])
                    familiesToAllAbunds[family3][i-1] = str(int(familiesToAllAbunds[family3][i-1]) + int(words[i]))

for family in familiesToAllAbunds:
    abunds = familiesToAllAbunds[family]
    fh3.write(family.strip() + '\t')
    fh3.write('\t'.join(abunds) + '\n')
fh3.close()
for family in familiesToAbunds:
    fh2.write('\t'.join([family.strip(),str(familiesToAbunds[family]),str(familiesToOccurs[family])]) + '\n')
fh2.close()
for family in allEnds:
    fh.write(family)
fh.close()






