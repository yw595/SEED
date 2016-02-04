#!/usr/bin/awk -f
BEGIN {
    FS="\t";
    OFS="\t";
}
{
    expressSum[$2]++;
}
END {
    for (i in expressSum) {
	print i, expressSum[i];
    }
}
