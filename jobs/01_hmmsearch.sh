#!/usr/bin/bash
#SBATCH -p short --time 2:00:00 --mem 2G --ntasks 4

module load hmmer/3
CPUS=2
if [ $SLURM_CPUS_ON_NODE ]; then
 CPUS=$SLURM_CPUS_ON_NODE
fi

# query are a folder full of .hmm files
EVALUE=1e-4
QUERYDIR=HMM
TARGETDIR=proteins
OUT=results
mkdir -p $OUT
for hmmfile in $QUERYDIR/*.hmm
do
 stem=$(basename $hmmfile .hmm)
 for db in $TARGETDIR/*.fasta
 do
  targetname=$(basename $db .fasta)
 # could add GA_CUT for 
  if [ ! -f $OUT/${stem}__$targetname.domtbl ]; then
   hmmsearch -E $EVALUE --cpu $CPUS --domtblout $OUT/${stem}__$targetname.domtbl $hmmfile $db > $OUT/${stem}__$targetname.hmmsearch
  fi
 done
done
 
