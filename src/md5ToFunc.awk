#!/usr/bin/awk -f
BEGIN {
    FS="\t";
    OFS="\t";
}
{
    md5ToFunc[$1]=md5ToFunc[$1] " " $3;
}
END {
    for (i in md5ToFunc) {
	print i, md5ToFunc[i];
    }
}
