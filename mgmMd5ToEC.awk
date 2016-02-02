#!/usr/bin/awk -f
BEGIN {
    FS="\t";
    OFS="\t";
}
FNR==NR {
    if (NR%1000000==0) {
	print NR > "/dev/stderr";
    }
    if (FNR<1000) {
	#print $1 > "/dev/stderr";
	#print md5ToEC[$2] > "/dev/stderr";
	#print md5ToEC[a50317a9d9ad88a9700900ed34fc4815] > "/dev/stderr";
    }
    md5ToEC[$1]=$2;
    next
}
{
    if (NR%1000000==0) {
	print NR > "/dev/stderr";
    }
    if (FNR<1000) {
	#print $2 > "/dev/stderr";
	#print md5ToEC[$2] > "/dev/stderr";
	#print md5ToEC[a50317a9d9ad88a9700900ed34fc4815] > "/dev/stderr";
    }
    if (md5ToEC[$2] ~ /[0-9|-][0-9|-]*\.[0-9|-][0-9|-]*\.[0-9|-][0-9|-]*\.[0-9|-][0-9|-]*/) {
	match(md5ToEC[$2],"[0-9|-][0-9|-]*\.[0-9|-][0-9|-]*\.[0-9|-][0-9|-]*\.[0-9|-][0-9|-]*");
	#express[md5ToEC[$2]]++;
	express[substr(md5ToEC[$2],RSTART,RLENGTH)]++;
    }
}
END {	
    for (i in express) {
	print i, express[i];
    }
}
