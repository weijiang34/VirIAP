# Virome-Identification-and-Analysis-Pipeline (VirIAP)

Welcome to use this pipeline!  

## Installiaton  
Download this directory:
```
git clone https://github.com/weijiang34/VirIAP.git
```
Install necessary tools and databases by running: 
```
cd VirIAP/
bash setup.sh --tools
bash setup.sh --databases
```
NOTE: Both these two steps requires internet connection, please make sure your PC/server node has the internet connection. The installation will take some time because there are lots of softwares and databases required, please wait patiently.

*(Known issue)* When installing, the pipeline will automatically setup some environmental variables in *src/envs.py*. You can also open *envs.py* with text editor to see if all the entries are filled. If not, you can try:
```
bash setup.py --check_envs
```
to reset the environmental variables.

## Workflow:  
### 1. Vreate your project
#### 1.1 Create a project:
The pipeline takes projects as its working directory. To create a project, please use the following command:  
```
python path/to/viriap/src/main.py -p [your_project_folder_path] create -i [path_to_your_fasta(s)_file]
```
___-p/--project_dir___: search/create a project directory according to the path (default: ./)  
___-i/--input___: Must be specified. A file contains a list of names of fasta files; OR one or more fasta files.  
This will create a project folder under the paht you provided. By default, it will take the current folder ("./") as the project folder.  
*Important*: if you are not willing to use the current folder as the project folder, please specify a project path, otherwise it will make your current folder messy (with many additional project-related folders/files created).  

After creating a project folder, cd to it, so that you don't need to provide the folder path every time you run a command:
```
cd [your_project_folder_path]
```

#### 1.2 Config your project
After creating a project, you will see a *config.yaml* file under your project folder, it records some important information of this project. Before going into generateig jobs, you need to specify some parameters in the *config.yaml* file:  
**For general PBS users**, please specify (keep others unchanged):  
```  
job_manager: pbs  
pbs:  
    ncpus: 32  
    mail_addr: 'your email address for receiving the status of jobs'  
```
**For gadi users**, please specify (keep others unchanged):  
```
job_manager: pbs  
pbs:  
    ncpus: 32  
    mail_addr: 'your email address for receviing the status of jobs'  
    gadi:  
        -l storage: 'storage of your project'  
        -P project: 'your project code'  
```
**(Optional)** You can also specify how many files to be included in a job, by specifying:  
```
max_batch_size: 10
``` 

### 2. Identify viruses
#### 2.1 Generate virus identification jobs
After successfuly configured your project, you are ready to generate jobs for identifying viruses:  
```
python path/to/viriap/src/main.py -p ./ main.py search --generate
```

#### 2.2 Submit jobs (manually):
In case of using improper resoures, please double check the job headers to make sure the resources required are valid/proper, and then submit jobs manually:
**For all PBS users**:
```
qsub path/to/your_project/jobs/job.pbs
```

#### 2.3 Check complete status:
```
python path/to/viriap/src/main.py -p ./ check
```

#### 2.4 Extract putative contigs:
```
python path/to/viriap/src/main.py -p ./ extract
```

#### 2.5 Decontamination (remove rRNA):
```
python path/to/viriap/src/main.py -p ./ decontam
```

#### 2.6 Merge confirmed contigs:
```
python path/to/viriap/src/main.py -p ./ merge
```

#### 2.7 Deduplication:
```
python path/to/viriap/src/main.py -p ./ dedup
```

#### 2.8 Quality check:
```
python path/to/viriap/src/main.py -p ./ checkv_quality
```

#### 2.9 Make OVUs by clustering:
```
python path/to/viriap/src/main.py -p ./ cluster
```

### 3. Mapping & Abundance
#### 3.1 Building mapping index:
```
python path/to/viriap/src/main.py -p ./ mapping --manifest path/to/your/manifest.csv --indexing
```
***NOTE*** *: Please submit the jobs manually after this step was finished.*
#### 3.2 Mapping reads to representative contigs:
```
python path/to/viriap/src/main.py -p ./ mapping --manifest path/to/your/manifest.csv --mapping
```
***NOTE*** *: Please submit the jobs manually after this step was finished.*
#### 3.3 Calculating relative abundance metrices:
```
python path/to/viriap/src/main.py -p ./ mapping --manifest path/to/your/manifest.csv --count_matrixcount_matrix
```
Relative abundance of OVUs will be in *"all_counts.csv", "all_FPKM.csv", "all_TPM.csv"*.

### 4. Classification of OVUs
#### 4.1 Classification with vContact3
```
python path/to/viriap/src/main.py -p ./ classify --generate_job 
```
***NOTE*** *: Please submit the jobs manually after this step was finished.*
#### 4.2 Merge lineages from CAT, vContact3, and GeNomad
```
python path/to/viriap/src/main.py -p ./ classify --merge_lineage 
```
This will output a summary, with lineages, of the OVUs, named *"OVU_info.csv"* under folder *"OVU/"*.

<!-- ## Modules 

-p/--project_dir:

1. will search/create a project directory according to the path (default: ./)

### create:

-i/--input: 

1. a file contains a list of names of fasta files; OR
2. one or more fasta files' paths

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
    .proj_dir/merged_confirmed_contigs.fasta -->