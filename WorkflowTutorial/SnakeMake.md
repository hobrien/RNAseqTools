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

- Pass sample names to script for paralle execution:
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

- Check exit status of each command:

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
- or use ```set -e``` to exit script after a command fails (see [here](http://www.davidpashley.com/articles/writing-robust-shell-scripts))

- What about:
    - inconsistently named inputs (including fastq files in different folders)
        - [bash paramter substitution](http://www.tldp.org/LDP/LG/issue18/bash.html)
        - ```find FastQ -name *.fastq.gz | sort | xargs -n 2 qsub MappingPipeline.sh```
    - requesting different resources for different steps of pipeline
        - [qsub -hold_jid argument](https://stackoverflow.com/questions/11525214/wait-for-set-of-qsub-jobs-to-complete)

### Make
- Used to compile [source code](https://github.com/alexdobin/STAR/tree/master/source) into binary
- Developed at a time when compilation was extremely resource-intensive
- Allows "nightly builds" where only modified code is re-compiled
- Builds a dependency tree from implicit wildcard rules
- Can be used to develop [bioinformatics pipelines](https://swcarpentry.github.io/make-novice)
- However: 
    - limited functionality and flexibility
    - perl-level syntax opacity
    - doesn't support parallelization
    
### Snakemake
![snakemake overview](https://snakemake.readthedocs.io/en/stable/_images/idea.png)
- Written in python
- Can be used to execute shell commands or python code blocks (in theory also [R code blocks](http://snakemake.readthedocs.io/en/stable/snakefiles/utils.html#scripting-with-r))
- Manages scheduling of job submission to [cluster](http://snakemake.readthedocs.io/en/stable/executable.html#cluster-execution) (or to the [cloud](http://snakemake.readthedocs.io/en/stable/executable.html#cloud-support))
    - pe smp and h_vmem can be specified as params in the Snakefile or in a [cluster config file](http://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#cluster-configuration)
    - cluster config files allow specification of default parameters
- Supports [config files](http://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html) to abstract details of pipeline from inputs and outputs
- Supports [benchmarking](http://snakemake.readthedocs.io/en/stable/tutorial/additional_features.html#benchmarking) to report CPU and memory usage
- Workflows can also be further abstracted by:
    - including [python code](http://snakemake.readthedocs.io/en/stable/project_info/faq.html#i-want-to-import-some-helper-functions-from-another-python-file-is-that-possible)
    - using the ```script``` command to [run code in a python script](http://snakemake.readthedocs.io/en/stable/tutorial/additional_features.html#using-custom-scripts), giving it access to variables defined in the Snakefile
    - including [additional Snakefiles](http://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#includes)
    - creating [sub-workflows](http://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#sub-workflows)
- [Conda environments](http://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management) can automatically be set up for each step of the analysis
- Many popular tools have [prewritten wrappers](https://snakemake-wrappers.readthedocs.io/en/stable) that automatically create the necessary environment and run the tools using the specified inputs, outputs, and paramaters

- Snakemake [documentation](https://snakemake.readthedocs.io/en/stable) and [tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html)
- Examples of:
    - a [Snakefile](https://github.com/hobrien/RNAseqTools/blob/master/Benchmarking/Snakefile)
    - including [additional Snakefiles](https://github.com/hobrien/RNAseqTools/blob/master/Benchmarking/bamQC)
    - a [config file](https://github.com/hobrien/RNAseqTools/blob/master/Benchmarking/config.yaml)
    - [cluster configuration](https://github.com/hobrien/RNAseqTools/blob/master/Benchmarking/cluster_config.yaml)
- use ```snakemake -n --dag | dot -Tsvg > dag.svg``` to produce a diagram of dependency tree:

![dag](https://github.com/hobrien/RNAseqTools/blob/master/Benchmarking/dag.png?raw=true)

### Alternatives to Snakemake
- Galaxy
![Galaxy BWA](http://galaxy.southgreen.fr/galaxy/static/style/cleaning_mapping_workflow_2.png)
- [Common Workflow Language (CWL)](https://github.com/common-workflow-language) / [Docker](https://www.docker.com/what-docker)
- See [here](https://academic.oup.com/bib/article-lookup/doi/10.1093/bib/bbw020) for more
