#!/usr/bin/bash

# this script takes results from the hmmsearch and runs a collecting script
# to build tables of counts for each domain

mkdir -p tables
TARGETDIR=proteins
TABLEFOLDER=tables
RESULTFOLDER=results
DOMAINSEQ=domains_seq

if [ -f config.txt ]; then
 source config.txt
fi

perl scripts/gather_domaincounts.pl -o $TABLEFOLDER/domain_counts.csv -i $RESULTFOLDER -db $TARGETDIR -c 1e-5  -d $DOMAINSEQ

