#!/usr/bin/awk -f
BEGIN { FS=OFS="\t" }
{
    coverage = $(NF-5) + 0
    similarity = $(NF-4) + 0
    if (coverage >= 90 && similarity >= 70) {
        print $0
    }
}
