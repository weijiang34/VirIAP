# Virus Identification Pipeline

## Project directory

project-|--out/
        |--sample1/
        |--sample2/
        |--...
    |--jobs/
        |--job1.pbs
        |--...
    |--logs/
        |--job1.o
        |--...
    |--config.yaml
    |--completeness_status.csv
    


## - Batch Mode

### Modules 

-p/--project_dir:

1. must be specified, and will search/create a project directory according to the path (default: ./)

~~-c/--config~~:  
~~/# 1. a config file: .proj_dir/config.yaml~~

#### create:

-i/--input: 

1. a file contains a list of names of fasta files; OR
2. one or more fasta files

#### search:

\# A submodule which generate executable job files for a computing cluster  
\# generate jobs  
--generate:  
1. generate jobs for all samples
--batch_size:  
1. defines the number of samples a batch job contains (default=10)

\# submit jobs  
--submit:
1. submit all generated jobs


#### check:  
\# check completeness status of all samples
- input:  
    .proj_dir/
- output:
    .proj_dir/completeness_status.csv


#### extract:
\# extract putative viral contigs  
Built-in functions:  
\# a loop that extract putative contigs from all samples
- input:  
    .proj_dir/out/sample_name/
- output:  
    .proj_dir/out/sample_name/putative_contigs.fasta
    .proj_dir/out/sample_name/putative_summary.csv
- filtering strategy:  
    -l/--min_length:  
    1. minimum length for putative contigs. (default 3000)
    -c/confirmed_tools:  
    1. minimum number of tools to confirm. (default 2)


#### decontam:
\# Decontamination: filter out rRNAs from bac,euk,arc,mito
Built-in functions:  
\# a loop that extract putative contigs from all samples
- input:  
    .proj_dir/out/sample_name/putative_contigs.fasta
    .proj_dir/out/sample_name/putative_summary.csv
- output:  
    .proj_dir/out/sample_name/rRNAs.tsv

#### confirm:
\# output confirmed viral contigs  
Built-in functions:  
\# a loop that extract putative contigs from all samples
- input:  
    .proj_dir/out/sample_name/putative_contigs.fasta
    .proj_dir/out/sample_name/putative_summary.csv
    .proj_dir/out/sample_name/rRNAs.tsv
- output:  
    .proj_dir/out/sample_name/confirmed_contigs.fasta
    .proj_dir/out/sample_name/confirmed_summary.csv

#### merge:
\# merge all confirmed viral contigs into one fasta file  
- input:  
    .proj_dir/out
- output:  
    .proj_dir/merged_confirmed_contigs.fasta


## Resources

1. CAT - 192GB - 32cpus - 3k 48:00:00 - >10samples
2. VirSorter2 - 64GB - 32cpus - 2.3k 48:00:00 : 7-8samples
3. GeNomad - 64GB - 32cpus - 2.4k 48:00:00 : 7-8samples
4. ViraLM - 32GB - 2gpus - 16cpus - 12:00:00(7-8:00:00) : >3kUs 