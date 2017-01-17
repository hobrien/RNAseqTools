#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

BASEDIR=~/Benchmark
# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
if [ ! -d $BASEDIR/Tophat2 ]
then
    echo "Making Tophat2 folder"
    mkdir $BASEDIR/Tophat2
    if [ $? -eq 0 ]
    then
        echo "Finished making Tophat2 folder"
    else
        echo "Could not create Tophat2 folder"
        exit 1
    fi    
fi
    
if [ ! -f $BASEDIR/Tophat2/accepted_hits.bam  ] || [ ! -f $BASEDIR/Tophat2/unmapped.bam ]   
then
    echo "Running Tophat2 mapping"
    time tophat --keep-fasta-order --library-type fr-secondstrand --num-threads 8 \
      --transcriptome-index /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.inx \
      --output-dir $BASEDIR/$SampleID \
      /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Sequence/Bowtie2Index/genome \
      $BASEDIR/FastQ/test_R1.fastq.gz $BASEDIR/FastQ/test_R2.fastq.gz > $BASEDIR/Tophat2/time.txt
    if [ $? -eq 0 ]
    then
        echo "Finished Tophat mapping "
    else
        echo "Tophat mapping failed"
        exit 1
    fi    
fi
