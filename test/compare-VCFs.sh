#!/bin/bash
# compare-VCFs.sh

vcf1=$1
vcf2=$2

process_vcf() {
    zcat "$1" | grep -v '^##' | sort -k1,1 -k2,2n
}

# Run diff and check if it fails
if ! diff -q <(process_vcf "$vcf1") <(process_vcf "$vcf2"); then
    echo "Warning: The VCF files are different!"
    exit 1
fi
