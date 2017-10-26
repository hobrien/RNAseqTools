# System Configuration
### Three ways to install software
- Build from source
- Precompiled binary
- Package manager
    - OSX: macPorts, fink
    - [Linux](https://www.digitalocean.com/community/tutorials/package-management-basics-apt-yum-dnf-pkg): yum (.rpm), apt-get (.deb)
    - pip, CRAN, CPAN, bioconductor
    - homebrew, conda
- Conda distributions:
    - [conda vs. anaconda vs. miniconda](https://bioconda.github.io/faqs.html#conda-anaconda-minconda)
        -"So: conda is a package manager, Miniconda is the conda installer, and Anaconda is a scientific Python distribution that also includes conda."
    - [Anaconda scientific python stack](https://docs.anaconda.com/anaconda/packages/pkg-docs): numpy, scipy, pandas, matplotlib, seaborn, etc
    - [r-essentials](https://conda.io/docs/user-guide/tasks/use-r-with-conda.html): "Anaconda for R"
    - [Bioconda](https://bioconda.github.io): bioinformatics-specific python packages (biopython, pysam, etc) + other bioinformatics tools (bwa, samtools, etc)
    - [condafordge](https://conda-forge.org): community contributed packages (*proceed with caution!*)

### Setting up your environment
- Copy executables to ~/bin
- Create symlinks to executables in ~/bin
    - ```ln -s ~/src/plink_mac/plink ~/bin/plink```
- Add locations of executables to your PATH variable
    - ```export PATH=~/src/plink_mac:$PATH```
- Use an environment management system

### Environmental Variables
- Variables that are available to subprocesses:
    - ```MYNAME=Heath```
    - ```echo $MYNAME```
        - "Heath"
    - ```echo "echo \$MYNAME" > test.sh```
    - ```bash test.sh```
        - ""
    - ```export MYNAME=Heath```
    - ```bash test.sh```
        - "Heath"
- Use environmental variables sparingly. It's usually better to pass variables to subprocesses as arguments:
    - ```echo "echo \$1" > test2.sh```
    - ```bash test2 $NAME```
        - "Heath"

- Environmental variables can be set automatically when a session is launched:
    - ```echo "MYNAME=Heath" >> ~/.bash_profile```
    - (start new session)
    - ```echo $NAME```
        - Heath
    - ```bash test.sh```
        - ""
    - ```echo "export MYNAME" >> ~/.bash_profile```
    - (start new session)
    - ```bash test.sh```
        - "Heath"
- If you want commands to execute every time you run a bash script, add the to .bashrc (but not on OSX)
- If you want commands to run in other shells, add them to .profile

### Making your scripts executable:
- In addition to being in your PATH, scripts need to be executable and to start with a "shebang" line:
    - ```./test.sh```
        - "-bash: ./test.sh: Permission denied"
    - ```chmod +x test.sh```
    - ```./test.sh```
        - "Heath"
    - ```echo print\(\"Heath\"\) > test.py```
    - ```python test.py```
        - "Heath"
    - ```chmod +x test.py```
    - ```./python.py```
        - "./test.py: line 1: syntax error near unexpected token `"Heath"'"  
    - ```which python```
        - "/Users/heo3/anaconda3/envs/python2/bin/python"
    - ```echo \#\!$(which python) > test.py```
    - ```echo print\(\"Heath\"\) >> test.py``` 
    - ```chmod +x test.py```
    - ```./test.py```
        - "Heath"
    - ```scp test.py mpmho@ravenlogin.arcca.cf.ac.uk:```
    - ```ssh mpmho@ravenlogin.arcca.cf.ac.uk "./test.py"```
        - "bash: ./test.py: /Users/heo3/anaconda3/envs/python2/bin/python: bad interpreter: No such file or directory"
    - ```echo \#\!/usr/bin/env\ python > test.py```
    - ```echo print\(\"Heath\"\) >> test.py``` 
    - ```chmod +x test.py```
    - ```./test.py```
       - "Heath"
    - ```scp test.py mpmho@ravenlogin.arcca.cf.ac.uk:```
    - ```ssh mpmho@ravenlogin.arcca.cf.ac.uk "./test.py"```
       - "Heath"

### Conda environments
- Not all version of software work together
- The most obvious example of this is python2 vs python3, but there are plenty of others
- Results aren't always the same depending on the version of software used
- This can be the source of very subtile bugs
- Its a therefore a good idea to isolate your environment before changing anything
- Managing python:
    - ```conda create -n py35 python=3.5```
    - ```source activate py35```
    - ```source deactivate py35```
- List environments:
    - ```conda info --envs``` (slow)
    - ```ls ~/anaconda/envs``` (fast)
- Backing up your environment
    - creating a clone:
        - ```conda create --name py35_backup --clone py35```
    - creating a spec file (can be used to recreate your environment on a different computer with the *same OS*):
        - ```conda list --explicit > spec-file.txt```
    - exporting to a yml file (can be used to recreate your environment on any computer):
        - ```conda env export > environment.yml```
- Roll back:
    - ```conda remove --name myenv --all```
    - from a clone:
        - ```conda create --name py35 --clone py35_backup```
    - from a spec file:
        - ```conda create --name py35_new --file spec-file.txt```
    - from a yml file:
        - ```conda env create -f environment.yml```
- See [here](https://conda.io/docs/user-guide/tasks/manage-environments.html) for more details

### Bash sessions
- If you are working on multiple projects, it can be useful to use different sessions
- This will allow you to quickly change environments, and also ensure that things like exported variables and your bash history are distinct
- You can open multiple terminal sessions, but this can get confusing, and it doesn't allow you to return to your sessions when you get disconnected
- [tmux](https://gist.github.com/MohamedAlaa/2961058) makes this easy:
    - ```tmux new -s GENEX-FB1```
    - (ctrl+b d)
    - ```tmux new -s GENEX-FB2```
    - (ctrl+b d)
    - ```tmux attach -t GENEX-FB1```
    - ```tmux switch -t GENEX-FB2```

### Avoiding passwords
- [ssh key pairs](http://www.linuxproblem.org/art_9.html) can be used to authenticate server connections:
    - ```ssh-keygen -t rsa```
    - ```ssh mpmho@ravenlogin.arcca.cf.ac.uk mkdir -p .ssh```
    - ```cat .ssh/id_rsa.pub | ssh mpmho@ravenlogin.arcca.cf.ac.uk 'cat >> .ssh/authorized_keys'```
    
### Further reading
- [conda myths and misconceptions](https://jakevdp.github.io/blog/2016/08/25/conda-myths-and-misconceptions)
- [why I promote conda](http://technicaldiscovery.blogspot.co.uk/2013/12/why-i-promote-conda.html)
- [a tmux crash course](https://robots.thoughtbot.com/a-tmux-crash-course)