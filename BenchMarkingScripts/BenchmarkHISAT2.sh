#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 8
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH
MAPPER=HISAT2
BASEDIR=~/Benchmark/$MAPPER
REF=/c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta
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

if [ ! -f $REF/genome.1.ht2 ] && [ ! -f $REF/genome.1.ht21 ]
then
    echo "Building index for $MAPPER"
    hisat2-build -p 8 $REF/genome.fa $REF/genome
    if [ $? -eq 0 ]
    then
        echo "Finished building index for $MAPPER"
    else
        echo "Building index for $MAPPER failed"
        exit 1
    fi    
fi

if [ ! -f /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/splicesites.txt ]
then
    echo "Making list of known splice sites"
    hisat2_extract_splice_sites.py \
      /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.gtf \
      > /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/splicesites.txt
    if [ $? -eq 0 ]
    then
        echo "Finished making list of known splice sites"
    else
        echo "Making list of known splice sites failed"
        exit 1
    fi    
fi

      

if [ ! -f ~/RNAseqTools/BenchMarks/${MAPPER}_time.txt ]   
then
    echo "Running $MAPPER mapping"
    /usr/bin/time -v -o ~/RNAseqTools/BenchMarks/${MAPPER}_time.txt hisat2 --fr --threads 8 \
      -x $REF/genome \
      -1 ~/Benchmark/FastQ/test_R1.fastq.gz -2 ~/Benchmark/FastQ/test_R2.fastq.gz \
      | samtools view -S -bo $BASEDIR/accepted_hits.bam -
    if [ $? -eq 0 ]
    then
        echo "Finished $MAPPER mapping "
    else
        echo "$MAPPER mapping failed"
        exit 1
    fi    
fi

if [ ! -f ~/RNAseqTools/BenchMarks/${MAPPER}_stats.txt ]   
then
    echo "Running RSeQC stats on $MAPPER mapping"
   bam_stat.py -i $BASEDIR/accepted_hits.bam > ~/RNAseqTools/BenchMarks/${MAPPER}_stats.txt
    if [ $? -eq 0 ]
    then
        echo "Finished running RSeQC stats on $MAPPER mapping "
    else
        echo "RSeQC stats could not be run on $MAPPER mapping"
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
