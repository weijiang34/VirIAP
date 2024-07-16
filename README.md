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
NOTE: Both these two steps requires internet connection, please make sure your PC/server node has the ineernet connection. The installation will take some time because there are lots of softwares and databases required, please wait atiently.

<!-- ***NOTE:***  
The installation requires the following steps, if the installation failed, you may aslo follow these steps.  
Prerequisites:  
1. Install CAT from: https://github.com/MGXlab/CAT_pack
2. Install Virsorter2 from: https://github.com/jiarong/VirSorter2
3. Install GeNomad from: https://portal.nersc.gov/genomad/installation.html
4. Install ViraLM from: https://github.com/ChengPENG-wolf/ViraLM

Create environment for VIP:
```
conda create -n vip -c bioconda -c conda-forge seqkit checkv barrnap pandas
conda activate vip 
```
Download necessary databases:
```
# CAT DB (in this pipeline, we used NR database):
mkdir ./dependencies/CAT_pack_nr_db 
cd ./dependencies/CAT_pack_nr_db
wget tbb.bio.uu.nl/tina/CAT_pack_prepare/20240422_CAT_nr.tar.gz.
tar -xvzf 20240422_CAT_nr.tar.gz
cd ../..
# checkv DB:
checkv download_database ./dependencies/checkvdb
``` -->


## Usage:  
### Create a project:
The pipeline takes projects as its working directory. To create a project, please use the following command:  
```
python main.py -p [your_project_folder_path] create -i [path_to_your_fasta(s)_file]
```
___-p/--project_dir___: search/create a project directory according to the path (default: ./)  
___-i/--input___: Must be specified. A file contains a list of names of fasta files; OR one or more fasta files.  
This will create a project folder under the paht you provided. By default, it will take the current folder ("./") as the project folder.  
__Important__: if you are not willing to use the current folder as the project folder, please specify a project path, otherwise it will make your current folder messy (with many additional project-related folders/files created).  
After creating a project folder, cd to it, so that you don't need to provide the folder path every time you run a command:
```
cd [your_project_folder_path]
```
### Generate jobs:

```
python -p ./ main.py search --generate
```

### Submit jobs (optional):

```
python -p ./ main.py search --submit
```

### Check complete status:

```
python -p ./ main.py check
```

### Extract putative contigs:

```
python -p ./ main.py extract
```

### Decontamination (remove rRNA):

```
python -p ./ main.py decontam
```

### Extract confirmed contigs:

```
python -p ./ main.py confirm
```

### Merge confirmed contigs:

```
python -p ./ main.py merge
```



## Modules 

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
    .proj_dir/merged_confirmed_contigs.fasta