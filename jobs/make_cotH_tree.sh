#!/usr/bin/bash

#SBATCH --mem 16G --ntasks 1 --nodes 1 
module load hmmer/3
module load trimal
module load fasttree

cd domainseq
if [ ! -f CotH_genes.hmmalign ]; then
hmmalign ../HMM/CotH.hmm CotH_genes.fas  > CotH_genes.hmmalign
fi

if [ ! -f CotH_genes.msa ]; then
esl-reformat --replace=\*:-  --gapsym=- clustal CotH_genes.hmmalign > CotH_genes.msa
fi

if [ ! -f CotH_genes.msa2 ]; then
  esl-reformat --replace=x:- clustal CotH_genes.msa > CotH_genes.msa2
fi
if [ ! -f CotH_genes.trim1 ]; then
 trimal -resoverlap 0.50 -seqoverlap 60 -in CotH_genes.msa2 -out CotH_genes.trim1
fi
if [ ! -f CotH_genes.trim ]; then
  trimal -automated1 -in CotH_genes.msa2 -out CotH_genes.trim -fasta
fi

FastTreeMP -gamma -wag < CotH_genes.trim > CotH_genes.nj.tre
