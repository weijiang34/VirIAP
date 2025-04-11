
# VirIAP: Virome-Identification-and-Analysis-Pipeline

Welcome to use this pipeline!  

## Table of Contents
- [VirIAP: Virome-Identification-and-Analysis-Pipeline](#viriap-virome-identification-and-analysis-pipeline)
  - [Table of Contents](#table-of-contents)
  - [:rocket: New features](#rocket-new-features)
    - [**v0.2.0** - OVU (Operational Viral Units) updates!](#v020---ovu-operational-viral-units-updates)
  - [:gear: Installaton](#gear-installaton)
  - [:book: Workflow](#book-workflow)
    - [1. Create your project](#1-create-your-project)
      - [1.1 Create a project](#11-create-a-project)
      - [1.2 Configure your project](#12-configure-your-project)
    - [2. Identify viruses](#2-identify-viruses)
      - [2.1 Generate virus identification jobs](#21-generate-virus-identification-jobs)
      - [2.2 Submit jobs (manually)](#22-submit-jobs-manually)
      - [2.3 Check completion status](#23-check-completion-status)
      - [2.4 Extract putative contigs](#24-extract-putative-contigs)
      - [2.5 Decontamination (remove rRNA)](#25-decontamination-remove-rrna)
      - [2.6 Merge confirmed contigs](#26-merge-confirmed-contigs)
      - [2.7 Deduplication](#27-deduplication)
      - [2.8 Quality check](#28-quality-check)
    - [3. OVU construction, abundance classification, and classification](#3-ovu-construction-abundance-classification-and-classification)
      - [3.1 Construct OVUs](#31-construct-ovus)
      - [3.2 Esitmate abundance](#32-esitmate-abundance)
      - [3.3 Classification of OVUs](#33-classification-of-ovus)
  - [:scroll: Citation](#scroll-citation)
  - [:envelope: Contact](#envelope-contact)


## :rocket: New features  

### **v0.2.0** - OVU (Operational Viral Units) updates!  

**Now featuring:** Construction, classification, and abundance estimation of OVUs.
<details>
  <summary>More details</summary>  

**What’s new in v2?**  
This release focuses on **viral contig processing**, with enhanced:  
- ✔ **Clustering** → Group contigs into OVUs  
- ✔ **Mapping** → Assign reads for abundance estimation  
- ✔ **Classification** → Improved taxonomic annotation  

**Key benefits:**
- More accurate viral profiling
- Streamlined workflow for OVU-based analysis
- Better visualization of results
</details>

---  

For more updates information, please see [Updates](./updates.md).

## :gear: Installaton  
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
bash setup.sh --check_envs
```
to reset the environmental variables.

## :book: Workflow
### 1. Create your project
#### 1.1 Create a project
The pipeline takes projects as its working directory. To create a project, please use the following command:  
```
python path/to/viriap/src/main.py -p [your_project_folder_path] create -i [path_to_your_fasta(s)_file]
```
___-p/--project_dir___: search/create a project directory according to the path (default: ./)  
___-i/--input___: Must be specified. A file contains a list of names of fasta files; OR one or more fasta files.  
This will create a project folder under the path you provided. By default, it will take the current folder ("./") as the project folder.  
*Important*: if you are not willing to use the current folder as the project folder, please specify a project path, otherwise it will make your current folder messy (with many additional project-related folders/files created).  

After creating a project folder, cd to it, so that you don't need to provide the folder path every time you run a command:
```
cd [your_project_folder_path]
```

#### 1.2 Configure your project
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

#### 2.2 Submit jobs (manually)
In case of using improper resoures, please double check the job headers to make sure the resources required are valid/proper, and then submit jobs manually:
**For all PBS users**:
```
qsub path/to/your_project/jobs/job.pbs
```

#### 2.3 Check completion status
```
python path/to/viriap/src/main.py -p ./ check
```

#### 2.4 Extract putative contigs
```
python path/to/viriap/src/main.py -p ./ extract
```
Parameters for filtering:  
-l , defallt: 3000, minimum length for putative contigs.  
-c , default: 2, minimum number of tools to confirm a viral contig.  
--trusted , The trusted tool(s) to be used for classification, selected from {cat,vs2,gnm,vlm}, comma seperated. Default: 'cat'.  
--skip , List of tools to skip in the results, selected from {cat,vs2,gnm,vlm}, comma seperated. Default: None.

#### 2.5 Decontamination (remove rRNA)
```
python path/to/viriap/src/main.py -p ./ decontam
```

#### 2.6 Merge confirmed contigs
```
python path/to/viriap/src/main.py -p ./ merge
```

#### 2.7 Deduplication
```
python path/to/viriap/src/main.py -p ./ dedup
```

#### 2.8 Quality check
```
python path/to/viriap/src/main.py -p ./ checkv_quality
```

### 3. OVU construction, abundance classification, and classification 
#### 3.1 Construct OVUs 
```
python path/to/viriap/src/main.py -p ./ cluster
```
*If you find there are only small amount of viral contigs detected, you may use --no_checkv option to directly clustering un-checkved contigs.*

#### 3.2 Esitmate abundance 
- Building mapping index
    ```
    python path/to/viriap/src/main.py -p ./ mapping --manifest path/to/your/manifest.csv --indexing
    ```
    ***NOTE*** *: Please submit the jobs manually after this step was finished.*
- Mapping reads to representative contigs
    ```
    python path/to/viriap/src/main.py -p ./ mapping --manifest path/to/your/manifest.csv --mapping
    ```
***NOTE*** *: Please submit the jobs manually after this step was finished.*
- Calculating relative abundance metrices
    ```
    python path/to/viriap/src/main.py -p ./ mapping --manifest path/to/your/manifest.csv --count_matrixcount_matrix
    ```
    Relative abundance of OVUs will be in *"all_counts.csv", "all_FPKM.csv", "all_TPM.csv"*.

#### 3.3 Classification of OVUs
- Classification with vContact3
    ```
    python path/to/viriap/src/main.py -p ./ classify --generate_job 
    ```
***NOTE*** *: Please submit the jobs manually after this step was finished.*
- Merge lineages from CAT, vContact3, and GeNomad
    ```
    python path/to/viriap/src/main.py -p ./ classify --merge_lineage 
    ```
    Optional parameter:  
    *--include ; used together with --merge_lineage, select one or more from [CAT, VCT, GNM], split with ','. If not specified, will use: ```--include CAT,VCT,GNM``` by default.*  
    
    The order of tool names determines priority, the former the higher, meaning it will consider the first tools annotation as a base reference and expand lineages as detailed (to a lower rank) as possible according to the following tools classification.  
    If you only wish to use one of the tools classification lineage (e.g. VCT), please specify as this: ```--include VCT```. However, it's recommmended to use all the three annotations (especially CAT), and just put the preffered tool at the first place. (e.g. ```--include VCT,CAT,GNM```)

This will output a summary, with lineages, of the OVUs, named *"OVU_info.csv"* under folder *"OVU/"*.

## :scroll: Citation

If you use VirIAP in your research, please cite:
```
```

## :envelope: Contact

For any questions or support, please contact:  
- Email: wjiang34-c@my.cityu.edu.hk

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
