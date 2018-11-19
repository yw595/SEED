import sys

normOTU = sys.argv[1]
if len(sys.argv)>2:
    z = int(sys.argv[2])

FI1 = open('/mnt/vdb/home/ubuntu2/justAGORADistsManual.txt')
numbers1 = []
numToSpec = {}
for line in FI1:
    words = line.split('|')
    if words[1]!='':
        numbers1.append(int(words[1]))
    numToSpec[int(words[1])] = words[0]
FI1.close()

FI2 = open(normOTU)
numbers2 = []
FI2.readline()
FI2.readline()
for line in FI2:
    numbers2.append(int(line.split()[0]))
FI2.close()

ggToFI = open('/mnt/vdb/home/ubuntu2/ggIDsToTax.txt')
ggTo = {}
for line in ggToFI:
    words = line.split()
    ggTo[int(words[0])] = int(words[1])
ggToFI.close()

if len(sys.argv)>2:
    pickFI = open('/mnt/vdb/home/ubuntu2/pickPresent'+str(z)+'.txt','w')
else:
    pickFI = open('/mnt/vdb/home/ubuntu2/pickPresent.txt','w')
count = 0
for number in numbers2:
    if number in ggTo:
        count+=1
    if number in ggTo and ggTo[number] in numbers1:
        pickFI.write(numToSpec[ggTo[number]].replace(' ','_')+'\n')
        #print(number)
        #print(numToSpec[ggTo[number]])
pickFI.close()
print(count)
