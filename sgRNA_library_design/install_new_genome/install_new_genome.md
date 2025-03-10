This file documents the steps to installing a new genome

1) Download GFF3/GTF/GFF and Genome Fasta file 
ex. from https://ftp.ensembl.org/pub/release-113/


2) Create the gffutils database

python3 install_gff3.py -G <assembly>.gff.gz

3) Create fai index for header

gunzip <assembly>.fa.gz
samtools faidx <assembly>.fa

4) Creade Auxiliairy info (NN bed file / 2bit / sizes)
faToTwoBit <assembly>.fa <assembly>.2bit
twoBitInfo <assembly>.2bit -nBed <assembly>_NN.bed # used to determine validity
awk '!seen[$2]++ {print $1"\t"$2}' <assembly>.fa.fai > <assembly>.sizes # used to determine validity

5) Create a R BsGenome Should one not already exist
forge_bsgenome.R <path> <genus> <species> <assembly> <mt> <library_dir>

6) Download indexed vep cache availlable at ensembl

for release 113 -> https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/
tar xzf <vep cache>.tar.gz
rm <vep cache>.tar.gz

7) Data Organisation

starting from the parameter pipeline_gen

pipeline_gen -> genomes -> <assembly> -> <assembly>._NN.bed
                                      -> <assembly>.db
                                      -> <assembly>.sizes
                               ...
             -> bowtie_index -> <assembly>.1.ebwt
                                    ...

only .db .NN.bed .sizes and bowtie index are used.

starting form VEP_dir

VEP_dir -> vep.sif
        -> vep_cache -> <assembly>
