# Workflow Management

### Reproducible science
- whenever possible, automate *all* steps of analysis to keep a record of everything that is done and to avoid introduction of errors
- Reduce opportunities for errors to enter workflow
- Avoid errors where different results are based on different versions of files
- Reduce effort required to correct errors or improve methods
- Adapt workflow to new datasets
- Promote data/knowledge sharing
- Make it easy to crowdsource error correction/pipeline improvements

### Workflow Solutions:
- Bash
- Make
- Snakemake
    - Keeps track of dependencies between files
    - Manages scheduling of job submission to cluster
    - Automatically re-runs necessary steps when part of the pipeline changes
    - Supports config files to abstract details of pipeline from inputs and outputs
    - Include python code within pipeline
    
