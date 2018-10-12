inFI = open('/mnt/vdb/home/ubuntu2/data_and_code/processed_data_tables/data.r')
inFI.readline()
samplesToData = {}
samplesToStatus = {}
species = []
for line in inFI:
    words = line.strip().split('\t')
    if words[0] not in samplesToStatus:
        samplesToStatus[words[0]] = words[2][1:len(words[2])-1]
    if words[0] not in samplesToData:
        samplesToData[words[0]] = {}
    if words[5]!='NA' and float(words[5])!=0 and (words[3]=='Genus' or words[3]=='Family'):
        speciesname = words[4]
        speciesname = speciesname[1:len(speciesname)-1]
        if speciesname not in species:
            species.append(speciesname)
        if speciesname not in samplesToData[words[0]]:
            samplesToData[words[0]][speciesname] = float(words[5])
        else:
            samplesToData[words[0]][speciesname] = samplesToData[words[0]][speciesname] + float(words[5])
inFI.close()
outFI = open('/mnt/vdb/home/ubuntu2/taxonomyDiabetes.tsv','w')
outFI2 = open('/mnt/vdb/home/ubuntu2/diabetesctrl.txt','w')
outFI3 = open('/mnt/vdb/home/ubuntu2/diabetestype2nomet.txt','w')
outFI4 = open('/mnt/vdb/home/ubuntu2/diabetestype2withmet.txt','w')
samples = samplesToData.keys()
outFI.write('#OTU ID\t')
print(len(species))
for i in range(len(samples)):
    sample = samples[i]
    if samplesToStatus[sample]=='ND CTRL':
        outFI2.write(sample+'\n')
    if samplesToStatus[sample]=='T2D metformin-':
        outFI3.write(sample+'\n')
    if samplesToStatus[sample]=='T2D metformin+':
        outFI4.write(sample+'\n')
    if i!=len(samples)-1:
        outFI.write(sample+'\t')
    else:
        outFI.write(sample+'\n')
for ithspecies in species:
    outFI.write(ithspecies+'\t')
    for i in range(len(samples)):
        sample = samples[i]
        delimiter = ''
        if i!=len(samples)-1:
            delimiter='\t'
        else:
            delimiter='\n'
        if ithspecies not in samplesToData[sample]:
            outFI.write('0.0'+delimiter)
        else:
            outFI.write(str(samplesToData[sample][ithspecies])+delimiter)
outFI.close()
outFI2.close()
outFI3.close()
outFI4.close()
