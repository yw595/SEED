#!/usr/bin/gawk -f
BEGIN {
    FS="\t";
    OFS="\t";
    command="echo";
}
FNR==NR {
    if ($2 ~ /.*[0-9|-][0-9|-]*\.[0-9|-][0-9|-]*\.[0-9|-][0-9|-]*\.[0-9|-][0-9|-]*.*/) {
	funcToEC[$1]=$2;
    }
    if (NR%1000000==0) {
	print NR > "/dev/stderr";
    }
    next
}
{
    if (funcToEC[$3] ~ /.*[0-9|-][0-9|-]*\.[0-9|-][0-9|-]*\.[0-9|-][0-9|-]*\.[0-9|-][0-9|-]*.*/) {
	md5ToEC[$1]=md5ToEC[$1] " " funcToEC[$3];
    }
    if (NR%1000000==0) {
	print NR > "/dev/stderr";
    }
}
END {

    #for (i in funcToEC) {
	#print i, funcToEC[i];
    #}

    for (i in md5ToEC) {
	print i, md5ToEC[i];
    }
}
