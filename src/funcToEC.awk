#!/usr/bin/awk -f
BEGIN {
    FS="\t";
    OFS="\t";
}
{
    if ($2 ~ /.*[0-9|-][0-9|-]*\.[0-9|-][0-9|-]*\.[0-9|-][0-9|-]*\.[0-9|-][0-9|-]*.*/) {
	funcToEC[$1]=funcToEC[$1] " " $2;
    }
    #if ($2 ~ /3\./) {
	#funcToEC[$1]=funcToEC[$1] FS $2;
    #}
}
END {
    for (i in funcToEC) {
	print i, funcToEC[i];
    }
}
