inFI = open('/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/taxonomic_profiles.tsv')
line = inFI.readline().strip('\t')
taxIDs = line.split()[1:-1]
print(taxIDs)
inFI = open('/mnt/vdb/home/ubuntu2/MATLAB/SEED/input/hmp2_metadata.csv')
for line in inFI:
    words = line.strip().split(',')
    print(words)
    nonsense += 1
    if words[1] in taxIDs:
        pass
        #print(line)
inFI.close()
