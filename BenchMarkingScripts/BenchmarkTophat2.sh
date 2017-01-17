#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

BASEDIR=~/Benchmark/Tophat2
# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
if [ ! -d $BASEDIR ]
then
    echo "Making Tophat2 folder"
    mkdir $BASEDIR
    if [ $? -eq 0 ]
    then
        echo "Finished making Tophat2 folder"
    else
        echo "Could not create Tophat2 folder"
        exit 1
    fi    
fi
    
if [ ! -f ~/RNAseqTools/BenchMarks/Tophat2_time.txt ]   
then
    echo "Running Tophat2 mapping"
    /usr/bin/time -v -o ~/RNAseqTools/BenchMarks/Tophat2_time.txt tophat --keep-fasta-order --library-type fr-secondstrand --num-threads 8 \
      --transcriptome-index /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.inx \
      --output-dir $BASEDIR \
      /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Sequence/Bowtie2Index/genome \
      ~/Benchmark/FastQ/test_R1.fastq.gz ~/Benchmark/FastQ/test_R2.fastq.gz 
    if [ $? -eq 0 ]
    then
        echo "Finished Tophat mapping "
    else
        echo "Tophat mapping failed"
        exit 1
    fi    
fi

if [ ! -f ~/RNAseqTools/BenchMarks/Tophat2_flagstat.txt ]   
then
    echo "Running samtools flagstat Tophat2 mapping"
    samtools flagstat -n $BASEDIR/accepted_hits.bam \
      > ~/RNAseqTools/BenchMarks/Tophat2_flagstat.txt
    if [ $? -eq 0 ]
    then
        echo "Finished running samtools flagstat on Tophat mapping "
    else
        echo "samtools flagstat could not be run on tophat mapping"
        exit 1
    fi    
fi


if [ ! -f $BASEDIR/accepted_hits_nsort.bam ]   
then
    echo "Sorting Tophat2 mapping"
    samtools sort -n $BASEDIR/accepted_hits.bam \
      $BASEDIR/accepted_hits_nsort 
    if [ $? -eq 0 ]
    then
        echo "Finished sorting Tophat mapping "
    else
        echo "Could not sort tophat mapping"
        exit 1
    fi    
fi

if [ ! -f ~/RNAseqTools/BenchMarks/Tophat2_counts.txt ]   
then
    echo "Running htseq-count on sorted Tophat2 mapping"
    htseq-count -f bam -s reverse -t exon -i gene_id -m intersection-strict \
      $BASEDIR/accepted_hits_nsort.bam \
      /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.gtf \
      > ~/RNAseqTools/BenchMarks/Tophat2_counts.txt
    if [ $? -eq 0 ]
    then
        echo "Finished running htseq-count on sorted Tophat mapping "
    else
        echo "htseq-count could not be run on sorted tophat mapping"
        exit 1
    fi    
fi

