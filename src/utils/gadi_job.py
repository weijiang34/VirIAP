import pandas as pd
import numpy as np
import os
import warnings
import envs

# job configuration class 
class GadiJobConfig():
    def __init__(self, job_name: str=None, walltime: str=None, mem: str=None, ncpus: str=None, ngpus: str=None, 
                node_type: str=None, jobfs: str=None, project: str=None, storage: str=None, 
                log_dir: str=None, mail_addr: str=None):
        
        # specific parameters
        self.job_name = job_name
        self.node_type = node_type
        self.walltime = walltime
        self.mem = mem
        self.ncpus = ncpus
        self.ngpus = ngpus
        
        # default parameters
        self.jobfs = jobfs
        self.project = project
        self.storage = storage
        self.log_dir = log_dir
        self.log_o = os.path.join(log_dir, f"{self.job_name}.o")
        self.log_e = os.path.join(log_dir, f"{self.job_name}.e")
        self.mail_addr = mail_addr

class GadiJob():
    pbs_head = ["#!/bin/bash\n"]
    # initialize a job with job configuration and command lines
    def __init__(self, index_s: str=None, index_e: str=None, tool_name: str=None, job_name: str=None, config_template=None, commands_template=None, 
                 walltime: str=None, mem: str=None, ncpus: str=None, ngpus: str=None,  
                 node_type: str=None, jobfs: str=None, project: str=None, storage: str=None, 
                 log_dir: str=None, mail_addr: str=None, file_list: list=[], config: dict={},
                 out_dir: str=os.path.join(os.getcwd(),"out"), prj_dir: str=os.getcwd()):
        
        # quick-run config templates to select
        config_templates = {
            "CAT_10_nr":GadiJobConfig(job_name=f"CAT_{index_s}_{index_e}", walltime="48:00:00", mem="192GB", ncpus="32", ngpus="0", 
                                        node_type="normalsl", jobfs="2GB", project=config["gadi"]["-P project"], storage=config["gadi"]["-l storage"], 
                                        log_dir=os.path.join(prj_dir,"logs"), mail_addr=config["gadi"]["-M mail_addr"]), 
            "VirSorter2_10":GadiJobConfig(job_name=f"VS2_{index_s}_{index_e}", walltime="48:00:00", mem="192GB", ncpus="32", ngpus="0", 
                                        node_type="normalsl", jobfs="2GB", project=config["gadi"]["-P project"], storage=config["gadi"]["-l storage"], 
                                        log_dir=os.path.join(prj_dir,"logs"), mail_addr=config["gadi"]["-M mail_addr"]),
            "geNomad_10":GadiJobConfig(job_name=f"GNM_{index_s}_{index_e}", walltime="48:00:00", mem="192GB", ncpus="32", ngpus="0", 
                                        node_type="normalsl", jobfs="2GB", project=config["gadi"]["-P project"], storage=config["gadi"]["-l storage"], 
                                        log_dir=os.path.join(prj_dir,"logs"), mail_addr=config["gadi"]["-M mail_addr"]),
            "ViraLM_10_GPU":GadiJobConfig(job_name=f"VLM_{index_s}_{index_e}", walltime="12:00:00", mem="64GB", ncpus="32", ngpus="2", 
                                        node_type="dgxa100", jobfs="2GB", project=config["gadi"]["-P project"], storage=config["gadi"]["-l storage"], 
                                        log_dir=os.path.join(prj_dir,"logs"), mail_addr=config["gadi"]["-M mail_addr"]),
        }
        
        # valid options, some parameters' value can only be selected from these valid options
        valid = {"tools":["CAT", "VS2", "GNM", "VLM"]}
        # Three ways to name a job: 
        #     1. use a config template; 
        #     2. sepecify tool_name, index_s, and index_e to auto-generate a job_name; 
        #     3. specify job_name directly; 
        if config_template==None:
            self.job_configs = GadiJobConfig(
                job_name=None, walltime=walltime, mem=mem, ncpus=ncpus, ngpus=ngpus, 
                node_type=node_type, jobfs=jobfs, project=project, storage=storage, log_dir=log_dir, mail_addr=mail_addr
            ) # GadiJobConfig Class
            if job_name==None:
                if tool_name in valid["tools"]:
                    self.job_configs.job_name = f"{tool_name}_{index_s}_{index_e}"
                elif tool_name not in valid["tools"]:
                    raise ValueError(f"{tool_name} must be one of: {valid['tools']}.")
            else: self.job_configs.job_name = job_name
        elif config_template!=None:
            if config_template in list(config_templates.keys()): self.job_configs = config_templates[config_template]
            else: raise ValueError(f"config_template '{config_template}' must be one of: {list(config_templates.keys())}.")
        
        # quick-run commands templates to select
        commands_templates = {
            "CAT":[
                f'threads={self.job_configs.ncpus}',
                f'CAT_dbPath={envs.CAT_DB_PATH}', # parameter:
                '',
                "file_list=(",
                "{}".format('\n'.join(f'"{item}"' for item in file_list)),
                ")",
                'for file in ${file_list[@]};do',
                '\techo "$file"',
                '',
                "\tfileHeader=$(echo $file | awk -F/ '{print $(NF-1)}')", # fileHeader
                f'\tout_dir={out_dir}/$fileHeader/CAT_results',
                '',
                '\tif [ -f "$out_dir/$fileHeader.nr.contig2classification.with_names.txt" ];then',
                '\t\techo "$fileHeader CAT has already finished. Continue to the next one."',
                '\t\tcontinue',
                '\tfi',
                '',
                '\tif [ ! -d "$out_dir/" ]; then',
                '\t\tmkdir $out_dir/',
                '\tfi',
                '',
                f'    {envs.CAT_PATH}/CAT_pack contigs -c $file \\',
                '                -d $CAT_dbPath/db \\',
                '                -t $CAT_dbPath/tax \\',
                '                -n $threads \\',
                '                --force \\',
                '                -o $out_dir/$fileHeader.nr',
                '                ',
                f'    {envs.CAT_PATH}/CAT_pack add_names -i $out_dir/$fileHeader.nr.contig2classification.txt \\',
                '                  -o $out_dir/$fileHeader.nr.contig2classification.with_names.txt \\',
                '                  -t $CAT_dbPath/tax \\',
                '                  --force \\',
                '                  --only_official ',
                '                  ',
                f'    {envs.CAT_PATH}/CAT_pack summarise -c $file \\',
                '                  -i $out_dir/$fileHeader.nr.contig2classification.with_names.txt \\',
                '                  -o $out_dir/$fileHeader.nr.summary.txt',
                '',
                '    echo CAT on "$file" finished!',
                '    echo "############################################################" ',
                'done'
            ],
            "VS2":[
                f'threads={self.job_configs.ncpus}',
                '',
                # 'source /g/data1b/oo46/wj6768/miniconda3/bin/activate /g/data1b/oo46/wj6768/miniconda3/envs/vs2',
                ' ',
                "file_list=(",
                "{}".format('\n'.join(f'"{item}"' for item in file_list)),
                ")",
                'for file in ${file_list[@]};do',
                '    echo "$file"',
                "    fileHeader=$(echo $file | cut -d '/' -f 9)", # fileHeader
                f'    out_dir={out_dir}/$fileHeader/VirSorter2_results', 
                '',
                '\tif [ -f "$out_dir/$fileHeader-final-viral-score.tsv" ];then',
                '\t\techo "$fileHeader VirSorter2 has already finished. Continue to the next one."',
                '\t\tcontinue',
                '\tfi',
                '',
                '    virsorter run -w $out_dir -i $file -l $fileHeader --include-groups dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae -j $threads all --rm-tmpdir',
                '    echo VirSorter2 on "$file" finished!',
                '    echo "############################################################"',
                'done'
            ],
            "GNM":[
                f'threads={self.job_configs.ncpus}',
                f'GeNomad_dbPath={envs.GENOMAD_DB_PATH}', # parameter:
                '',
                # 'source /g/data1b/oo46/wj6768/miniconda3/bin/activate /g/data1b/oo46/wj6768/miniconda3/envs/genomad',
                ' ',
                "file_list=(",
                "{}".format('\n'.join(f'"{item}"' for item in file_list)),
                ")",
                'for file in ${file_list[@]};do',
                '    echo "$file"',
                "    fileHeader=$(echo $file | awk -F/ '{print $(NF-1)}')", # fileHeader
                f'    out_dir={out_dir}/$fileHeader/GeNomad_results', 
                '',
                '\tif [ -f "$out_dir/final.contigs_summary/final.contigs_virus_summary.tsv" ];then',
                '\t\tif [ -f "$out_dir/final.contigs_summary/final.contigs_plasmid_summary.tsv" ] ;then',
                '\t\t\techo "$fileHeader geNomad has already finished. Continue to the next one."',
                '\t\t\tcontinue',
                '\t\tfi',
                '\tfi',
                '',
                '    genomad end-to-end $file $out_dir/ $GeNomad_dbPath --restart --threads $threads --cleanup',
                '    echo GeNomad on "$file" finished!',
                '    echo "############################################################" ',
                'done'
            ],
            "VLM":[
                'threads=32',
                f'ViraLMPath={envs.VIRALM_PATH}',
                '',
                # 'source /g/data1b/oo46/wj6768/miniconda3/bin/activate /g/data1b/oo46/wj6768/miniconda3/envs/viralm',
                ' ',
                "file_list=(",
                "{}".format('\n'.join(f'"{item}"' for item in file_list)),
                ")",
                'for file in ${file_list[@]};do',
                '    echo "$file"',
                "    fileHeader=$(echo $file | awk -F/ '{print $(NF-1)}')", # fileHeader
                f'    out_dir={out_dir}/$fileHeader/ViraLM_results',
                '',
                '\tif [ -f "$out_dir/result_final.csv" ];then',
                '\t\techo "$fileHeader ViraLM has already finished. Continue to the next one."',
                '\t\tcontinue',
                '\tfi',
                '',
                '    cd $ViraLMPath',
                '    python viralm.py --input $file --output $out_dir/',
                '    echo ViraLM on "$file" finished!',
                '    echo "############################################################" ',
                'done'
            ]
        }
        # exception check of commands templates input
        if commands_template==None:
            warnings.warn(f"Commands template are empty, please specify a template to use: {list(commands_templates.keys())}.", stacklevel=2)
            self.command_lines = [] # currently only support using template commands
        elif commands_template!= None:
            if commands_template in list(commands_templates.keys()): self.command_lines = commands_templates[commands_template]
            else: raise ValueError(f"commands_template '{commands_template}' must be one of: {list(commands_templates.keys())}.")

        # create lines of texts for pbs jobs
        self.job_texts = GadiJob.pbs_head + [
            f"# Job Name:",
            f"#PBS -N {self.job_configs.job_name}",
            f"",
            f"# Project Info:",
            f"#PBS -P {self.job_configs.project}",
            f"#PBS -l storage={self.job_configs.storage}",
            f"",
            f"# Log Output:",
            f"#PBS -o {self.job_configs.log_o}",
            f"#PBS -e {self.job_configs.log_e}",
            f"#PBS -j oe",
            f"",
            f"# Mailing:",
            f"#PBS -m abe",
            f"#PBS -M {self.job_configs.mail_addr}",
            f"",
            f"# Resources Allocation:",
            f"#PBS -q {self.job_configs.node_type}",
            f"#PBS -l walltime={self.job_configs.walltime}",
            f"#PBS -l mem={self.job_configs.mem}",
            f"#PBS -l ncpus={self.job_configs.ncpus}",
            "" if self.job_configs.ngpus=="0" else f"#PBS -l ngpus={self.job_configs.ngpus}",
            f"#PBS -l jobfs={self.job_configs.jobfs}",
        ] + self.command_lines
            
    def preview(self):
        for line in self.job_texts:
            print(line)
    def generate_pbs_file(self, job_dir=os.getcwd()):
        # the default job dir is the current working dir
        # if job dir does not exist, create it
        if not os.path.exists(job_dir):
            os.mkdir(job_dir)
        self.job_texts = [line+"\n" for line in self.job_texts]
        with open(os.path.join(job_dir, f"gadi_job_{self.job_configs.job_name}.pbs"), 'w') as f:
            f.writelines(self.job_texts)

    # Experimental functions
    # def update_config(self, config):
    #     self.job_configs = config
    #     return self
    # def add_command(self, cmd):
    #     self.job_texts.append(cmd)
    #     print(f"Command added on line {len(self.job_texts)}.")
    #     return self
    # def del_command(self, line_index):
    #     self.job_texts.append(cmd)
    #     del self.job_texts[line_index]
    #     print(f"Line {len(self.job_texts)} removed.")
    #     return self

def chunk_dataframe(df: pd.DataFrame, size: int=10):
    num_of_chunks = len(df) // size
    remaining_data = len(df) % size
    chunks = []
    for i in range(num_of_chunks):
        chunk = df.iloc[i*size:(i+1)*size, :]
        chunks.append(chunk)
    if remaining_data!= 0:
        chunk = df.iloc[num_of_chunks*size:, :]
        chunks.append(chunk)
    return chunks

def generate_jobs(project_dir: str=os.getcwd(), batch_size=10, config: dict={}):
    df = pd.read_csv(os.path.join(project_dir,"completeness_status.csv"), header=0, index_col=None, sep=',')
    data_chunks = chunk_dataframe(df, size=batch_size)
    if config["server_type"]=="gadi":
        for chunk in data_chunks:
            # CAT
            gadi_job = GadiJob(
                prj_dir=project_dir, out_dir=os.path.join(project_dir,"out"), tool_name="CAT", index_s=chunk.index.to_list()[0]+1, index_e=chunk.index.to_list()[-1]+1,
                config_template="CAT_10_nr", commands_template="CAT", file_list=chunk["path"].to_list(), config=config
            )
            gadi_job.generate_pbs_file(job_dir=os.path.join(f"{project_dir}","jobs"))
            # VS2
            gadi_job = GadiJob(
                prj_dir=project_dir, out_dir=os.path.join(project_dir,"out"), tool_name="VS2", index_s=chunk.index.to_list()[0]+1, index_e=chunk.index.to_list()[-1]+1,
                config_template="VirSorter2_10", commands_template="VS2", file_list=chunk["path"].to_list(), config=config
            )
            gadi_job.generate_pbs_file(job_dir=os.path.join(f"{project_dir}","jobs"))
            # GNM
            gadi_job = GadiJob(
                prj_dir=project_dir, out_dir=os.path.join(project_dir,"out"), tool_name="GNM", index_s=chunk.index.to_list()[0]+1, index_e=chunk.index.to_list()[-1]+1,
                config_template="geNomad_10", commands_template="GNM", file_list=chunk["path"].to_list(), config=config
            )
            gadi_job.generate_pbs_file(job_dir=os.path.join(f"{project_dir}","jobs"))
            # VLM
            gadi_job = GadiJob(
                prj_dir=project_dir, out_dir=os.path.join(project_dir,"out"), tool_name="VLM", index_s=chunk.index.to_list()[0]+1, index_e=chunk.index.to_list()[-1]+1,
                config_template="ViraLM_10_GPU", commands_template="VLM", file_list=chunk["path"].to_list(), config=config
            )
            gadi_job.generate_pbs_file(job_dir=os.path.join(f"{project_dir}","jobs"))

def submit_jobs(project_dir, file_index: str=""):
    job_list = os.listdir(os.path.join(project_dir,"jobs"))
    bash_commands = [
        "job_list=(\n",
        "{}\n".format('\n'.join(f'"{job}"' for job in job_list)),
        ")\n",
        'for job in ${job_list[@]};do\n',
        '   qsub job\n',
        "\n",
    ]
    with open(os.path.join(project_dir,"submit_jobs_tmp.sh"), 'w') as f:
        f.writelines("#!/bin/bash\n")
        f.writelines(bash_commands)
    os.system("chmod +x ./submit_jobs_tmp.sh")
    os.system("./submit_jobs_tmp.sh")
    os.remove("./submit_jobs_tmp.sh")
    