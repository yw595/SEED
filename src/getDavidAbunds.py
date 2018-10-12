import subprocess

#C:parse DavidS9, based on FIRST AND SECOND COORDS OF PCA??, assign cluster numbers to foundClusts 1 or 2
inputFI2 = open('DavidS9.txt')
foundTaxa = []
foundClusts1 = []
foundClusts2 = []
foundTaxa1 = []
foundTaxa2 = []
for line in inputFI2:
    words = line.split('\t')
    if float(words[1]) < 0 and float(words[2]) > 0:
        foundClusts1.append(words[0])
    elif float(words[1]) > 0 and float(words[2]) < 0:
        foundClusts2.append(words[0])
inputFI2.close()

#C:parse DavidS8, find dot in taxonomic name, maybe single letter taxonomic rank abbrev, so strip off, or species name after, so append. Record currentClust number, tracking taxa in it and abunds, at switch to new currentClust, append maxTaxon to foundTaxa representatives of clusters, put in 1 or 2 if match foundClusts from above
inputFI1 = open('DavidS8.txt')
currentClust = ''
currentTaxa = []
currentAbunds = []
for line in inputFI1:
    words = line.split('\t')
    dotIdx = words[3].find('.')
    #if words[3][len(words[3])-1]=='f':
    foundTaxon = ''
    if dotIdx==len(words[3])-2:
        foundTaxon = words[3][0:dotIdx]
    else:
        foundTaxon = words[3][0:dotIdx] + ' ' + words[3][dotIdx+1:len(words[3])]
    if currentClust=='':
        currentClust = words[0]
    #print(currentClust)
    if currentClust==words[0]:
        currentTaxa.append(foundTaxon)
        currentAbunds.append(float(words[4]))
    else:
        print(currentClust)
        currentClust = words[0]
        print(words)
        maxAbund = max(currentAbunds)
        for i in range(len(currentAbunds)):
            if currentAbunds[i]==maxAbund:
                maxTaxon = currentTaxa[i]
        currentAbunds = []
        currentTaxa = []
        currentTaxa.append(foundTaxon)
        currentAbunds.append(float(words[4]))
        foundTaxa.append(maxTaxon)
        foundTaxa = list(set(foundTaxa))
        if words[0][8:len(words[0])] in foundClusts1:
            foundTaxa1.append(foundTaxon)
        if words[0][8:len(words[0])] in foundClusts2:
            foundTaxa2.append(foundTaxon)
    #print(currentClust)
inputFI1.close()

#C:print all foundTaxa representatives of clusters, and those matching foundClusts 1 or 2
outputFI = open('DavidTaxa.txt','w')
for taxon in foundTaxa:
    outputFI.write(taxon + '\n')
outputFI.close()
outputFI = open('DavidTaxa1.txt','w')
for taxon in foundTaxa1:
    outputFI.write(taxon + '\n')
outputFI.close()
outputFI = open('DavidTaxa2.txt','w')
for taxon in foundTaxa2:
    outputFI.write(taxon + '\n')
outputFI.close()
