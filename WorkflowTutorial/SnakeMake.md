# Workflow Management

### Reproducible science
- whenever possible, automate *all* steps of analysis to keep a record of everything that is done and to avoid introduction of errors
- Avoid errors where different results are based on different versions of files
- Reduce effort required to correct errors or improve methods
- Adapt workflow to new datasets
- Promote data/knowledge sharing
- Make it easy to crowdsource error correction/pipeline improvements

### Bash:
- Hard coded file name:

```
bwa mem -t 8 GRCh38Decoy/Sequence/BWAIndex/genome.fa \ 
     FastQ/sample1_R1.fastq.gz FastQ/sample1_R2.fastq.gz \
     | samtools view -S -bo Mappings/sample1.bam -
samtools sort - o Mappings/sample1_sort.bam Mappings/sample1.bam
```
         
- Sample names as arguments:

```
for sample in $@
do
    bwa_mem -t 8 GRCh38Decoy/Sequence/BWAIndex/genome.fa \
        FastQ/${sample}_R1.fastq.gz FastQ/${sample}_R2.fastq.gz \
        | samtools view -S -bo Mappings/${sample}.bam
    samtools sort - o Mappings/${sample}_sort.bam Mappings/${sample}.bam
done
```
    - ```cat SampleList.txt | xargs -n 1 qsub MappingPipeline.sh```
    
- File Tests:

```
for sample in $@
do
    if [ ! -f $BASEDIR/Mappings/${sample}.bam ]
    then
        bwa_mem -t 8 GRCh38Decoy/Sequence/BWAIndex/genome.fa \
            FastQ/${sample}_R1.fastq.gz FastQ/${sample}_R2.fastq.gz \
            | samtools view -S -bo Mappings/${sample}.bam
    if [ ! -f Mappings/${sample}_sort.bam ]
    then
        samtools sort - o Mappings/${sample}_sort.bam Mappings/${sample}.bam
    fi
done
```

```
for sample in $@
do
    if test FastQ/${sample}_R1.fastq.gz -nt $BASEDIR/Mappings/${sample}.bam \
        || test FastQ/${sample}_R2.fastq.gz -nt $BASEDIR/Mappings/${sample}.bam
    then
        bwa_mem -t 8 GRCh38Decoy/Sequence/BWAIndex/genome.fa \
            FastQ/${sample}_R1.fastq.gz FastQ/${sample}_R2.fastq.gz \
            | samtools view -S -bo Mappings/${sample}.bam
    fi
    if test Mappings/${sample}.bam -nt Mappings/${sample}_sort.bam
    then
        samtools sort - o Mappings/${sample}_sort.bam Mappings/${sample}.bam
    fi
done
```

- Check exit status:

```
for sample in $@
do
   if test FastQ/${sample}_R1.fastq.gz -nt $BASEDIR/Mappings/${sample}.bam \
       || test FastQ/${sample}_R2.fastq.gz -nt $BASEDIR/Mappings/${sample}.bam
   then
       bwa_mem -t 8 GRCh38Decoy/Sequence/BWAIndex/genome.fa \
           FastQ/${sample}_R1.fastq.gz FastQ/${sample}_R2.fastq.gz \
           | samtools view -S -bo Mappings/${sample}.bam
       if [ $? -eq 0 ]
       then
           echo "Finished bwa_mem for $sample"
       else
           echo "bwa_mem failed for $sample"
           exit 1
       fi
   fi
   if test Mappings/{sample}.bam -nt Mappings/{sample}_sort.bam
   then
       samtools sort - o Mappings/{sample}_sort.bam Mappings/{sample}.bam
       if [ $? -eq 0 ]
       then
           echo "Finished samtools sort for $sample"
       else
           echo "samtools sort failed for $sample"
           exit 1
       fi
   fi
done
```

- What about:
    - inconsistently named inputs (including fastq files in different folders)
        - [bash paramter substitution](http://www.tldp.org/LDP/LG/issue18/bash.html)
        - ```find FastQ -name *.fastq.gz | sort | xargs -n 2 qsub MappingPipeline.sh```
    - requesting different resources for different steps of pipeline
        - [qsub -hold_jid argument](https://stackoverflow.com/questions/11525214/wait-for-set-of-qsub-jobs-to-complete)

- Make
- Snakemake
    - Keeps track of dependencies between files
    - Manages scheduling of job submission to cluster
    - Automatically re-runs necessary steps when part of the pipeline changes
    - Supports config files to abstract details of pipeline from inputs and outputs
    - Include python code within pipeline
    
- params.num_cores vs. [cluster.num_cores](http://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#cluster-configuration)
- includes
- [wrappers](https://snakemake-wrappers.readthedocs.io/en/stable)
- dag
- benchmarks
- conda env
- [R](http://snakemake.readthedocs.io/en/stable/snakefiles/utils.html#scripting-with-r)
- Scripts

- [Star source code](https://github.com/alexdobin/STAR/tree/master/source)

[Pipeline framework review](https://academic.oup.com/bib/article-lookup/doi/10.1093/bib/bbw020)