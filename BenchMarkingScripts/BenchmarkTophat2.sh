#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH
MAPPER=Tophat2
BASEDIR=~/Benchmark/$MAPPER
# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
if [ ! -d $BASEDIR ]
then
    echo "Making $MAPPER folder"
    mkdir $BASEDIR
    if [ $? -eq 0 ]
    then
        echo "Finished making $MAPPER folder"
    else
        echo "Could not create $MAPPER folder"
        exit 1
    fi    
fi
    
if [ ! -f ~/RNAseqTools/BenchMarks/${MAPPER}_time.txt ]   
then
    echo "Running $MAPPER mapping"
    /usr/bin/time -v -o ~/RNAseqTools/BenchMarks/${MAPPER}_time.txt tophat --keep-fasta-order --library-type fr-secondstrand --num-threads 8 \
      --transcriptome-index /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.inx \
      --output-dir $BASEDIR \
      /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Sequence/Bowtie2Index/genome \
      ~/Benchmark/FastQ/test_R1.fastq.gz ~/Benchmark/FastQ/test_R2.fastq.gz 
    if [ $? -eq 0 ]
    then
        echo "Finished $MAPPER mapping "
    else
        echo "$MAPPER mapping failed"
        exit 1
    fi    
fi

if [ ! -f ~/RNAseqTools/BenchMarks/${MAPPER}_flagstat.txt ]   
then
    echo "Running samtools flagstat $MAPPER mapping"
    samtools flagstat -n $BASEDIR/accepted_hits.bam \
      > ~/RNAseqTools/BenchMarks/${MAPPER}_flagstat.txt
    if [ $? -eq 0 ]
    then
        echo "Finished running samtools flagstat on $MAPPER mapping "
    else
        echo "samtools flagstat could not be run on $MAPPER mapping"
        exit 1
    fi    
fi


if [ ! -f $BASEDIR/accepted_hits_nsort.bam ]   
then
    echo "Sorting $MAPPER mapping"
    samtools sort -n $BASEDIR/accepted_hits.bam \
      $BASEDIR/accepted_hits_nsort 
    if [ $? -eq 0 ]
    then
        echo "Finished sorting $MAPPER mapping "
    else
        echo "Could not sort $MAPPER mapping"
        exit 1
    fi    
fi

if [ ! -f ~/RNAseqTools/BenchMarks/${MAPPER}_counts.txt ]   
then
    echo "Running htseq-count on sorted $MAPPER mapping"
    htseq-count -f bam -s reverse -t exon -i gene_id -m intersection-strict \
      $BASEDIR/accepted_hits_nsort.bam \
      /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.gtf \
      > ~/RNAseqTools/BenchMarks/${MAPPER}_counts.txt
    if [ $? -eq 0 ]
    then
        echo "Finished running htseq-count on sorted $MAPPER mapping "
    else
        echo "htseq-count could not be run on sorted $MAPPER mapping"
        exit 1
    fi    
fi

