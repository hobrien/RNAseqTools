# Workflow Management

### Reproducible science
- Whenever possible, automate *all* steps of analysis to keep a record of everything that is done and to avoid introduction of errors
- This makes it easy to:
    - Avoid errors where different results are based on different versions of files
    - Reduce effort required to correct errors or improve methods
    - Adapt workflow to new datasets
    - Promote data/knowledge sharing
    - Crowdsource error correction/pipeline improvements

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
- Inputs, outputs, and [parameters](http://snakemake.readthedocs.io/en/stable/tutorial/advanced.html#step-4-rule-parameters) can be specified for each rule
- Need to include a "master rule" (usually called ```all```) which requires all of your desired outputs as input
- Intermediate files can be [automatically removed](http://snakemake.readthedocs.io/en/stable/tutorial/advanced.html#step-6-temporary-and-protected-files) once they are no longer needed
- Supports [benchmarking](http://snakemake.readthedocs.io/en/stable/tutorial/additional_features.html#benchmarking) to report CPU and memory usage and [Logging](http://snakemake.readthedocs.io/en/stable/tutorial/advanced.html#step-5-logging) of messages/errors
- Supports [config files](http://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html) to abstract details of pipeline from inputs and outputs
- Workflows can also be further abstracted by:
    - using ```include`` statements to import [python code](http://snakemake.readthedocs.io/en/stable/project_info/faq.html#i-want-to-import-some-helper-functions-from-another-python-file-is-that-possible)
    - using the ```script``` command to [execute a python script](http://snakemake.readthedocs.io/en/stable/tutorial/additional_features.html#using-custom-scripts), giving it access to variables defined in the Snakefile
    - using ```include``` statements to import rules from [other Snakefiles](http://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#includes)
    - creating [sub-workflows](http://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#sub-workflows)
- [Conda environments](http://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management) can automatically be set up for each step of the analysis
- Many popular tools have [prewritten wrappers](https://snakemake-wrappers.readthedocs.io/en/stable) that automatically create the necessary environment and run the tools using the specified inputs, outputs, and paramaters
- There is also a [repository](https://bitbucket.org/johanneskoester/snakemake-workflows) of example rules and workflows for NGS analyses

### Getting started with Snakemake
- Snakemake [documentation](https://snakemake.readthedocs.io/en/stable) and [tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html)
- Examples of:
    - a [Snakefile](https://github.com/hobrien/RNAseqTools/blob/master/Benchmarking/Snakefile)
    - including [additional Snakefiles](https://github.com/hobrien/RNAseqTools/blob/master/Benchmarking/bamQC)
    - a [config file](https://github.com/hobrien/RNAseqTools/blob/master/Benchmarking/config.yaml)
    - [cluster configuration](https://github.com/hobrien/RNAseqTools/blob/master/Benchmarking/cluster_config.yaml)

### Snakemake usage
- Do a dry run of workflow, printing commands to screen:
    - ```snakemake -np```
- Produce a diagram of dependency tree:
    - ```snakemake -n --dag | dot -Tsvg > dag.svg```

![dag](https://github.com/hobrien/RNAseqTools/blob/master/Benchmarking/dag.png?raw=true)

- Rerun rule (and all rules with it as a dependency):
    - ```snakemake -R RULENAME```
- Rerun on new input files:
    - ```snakemake -n -R `snakemake --list-input-changes` ```
- Rerun edited rules:
    - ```$ snakemake -n -R `snakemake --list-params-changes` ```
- Submit jobs to cluster:
    - ```snakemake --use-conda --cluster-config cluster_config.yaml --cluster "qsub -pe smp {cluster.num_cores} -l h_vmem={cluster.maxvmem}" -j 20```

- See [here](http://snakemake.readthedocs.io/en/stable/api_reference/snakemake.html) for additional command-line options

### Alternatives to Snakemake
- Galaxy
![Galaxy BWA](http://galaxy.southgreen.fr/galaxy/static/style/cleaning_mapping_workflow_2.png)
- [Common Workflow Language (CWL)](https://github.com/common-workflow-language) / [Docker](https://www.docker.com/what-docker)
- See [here](https://academic.oup.com/bib/article-lookup/doi/10.1093/bib/bbw020) for more

<br>
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br /><span xmlns:dct="http://purl.org/dc/terms/" href="http://purl.org/dc/dcmitype/Text" property="dct:title" rel="dct:type">Workflow Management</span> by <span xmlns:cc="http://creativecommons.org/ns#" property="cc:attributionName">Heath O'Brien</span> is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
