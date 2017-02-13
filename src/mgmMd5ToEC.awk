#!/usr/bin/gawk -f
BEGIN {
    FS="\t";
    OFS="\t";
}

FILENAME=="input/GreenblumObesities.txt" {
    GreenblumObesities[$1]=$3;
}

FILENAME=="md5ToEC.map2" {
    if (NR%1000000==0) {
	print NR > "/dev/stderr";
    }
    if (FNR<1000) {
	#print FILENAME;
    }
    md5ToEC[$1]=$2;
    next
}
FILENAME=="/home/user/Downloads/SEED/MGMData/mgm4448044.3.650.protein.sims" {
    if (NR%1000000==0) {
	print NR > "/dev/stderr";
    }
    if (md5ToEC[$2] ~ /[0-9|-][0-9|-]*\.[0-9|-][0-9|-]*\.[0-9|-][0-9|-]*\.[0-9|-][0-9|-]*/) {
	split(md5ToEC[$2], a, " ");
	delete hasEC;
	for (a1 in a) {
	    hasEC[a[a1]]=1;
	}
	for (a2 in hasEC) {
	    if (a2 ~ /[0-9|-][0-9|-]*\.[0-9|-][0-9|-]*\.[0-9|-][0-9|-]*\.[0-9|-][0-9|-]*/) {
		match($1,"MH[0-9]+",a);
		MHID=a[0];
		#print MHID;
		if (MHID=="") {
		    expressNorm[a2]++;
		}
		else {
		    if (GreenblumObesities[MHID]=="Y") {
			expressObese[a2]++;
		    }
		    if (GreenblumObesities[MHID]=="N") {
			expressNorm[a2]++;
		    }
		}
	    }
	}
    }
}
END {	
    for (i in GreenblumObesities) {
	print i, GreenblumObesities[i];
    }
    for (i in expressObese) {
	print i, expressObese[i]
	print i, expressObese[i] > "/home/user/Downloads/homeBackup/MATLAB/SEED/testObese.txt";
    }
    for (i in expressNorm) {
	print i, expressNorm[i] > "/home/user/Downloads/homeBackup/MATLAB/SEED/testNorm.txt";
    }
}
