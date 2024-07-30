import os
import envs
import pandas as pd

job_managers = ["pbs","gadi"]

class BashHeader():
    def __init__(self, job_name, ncpus) -> None:
        self.job_name = job_name
        self.ncpus = ncpus

    def get_header(self):
        header = []
        header.append(f"#!/bin/bash")

        header = [x+"\n" for x in header]
        return header

class PBSHeader(BashHeader):
    def __init__(self, job_name, ncpus, ngpus, mem, walltime, mail_addr, log_o, log_e) -> None:
        super().__init__(job_name, ncpus)
        self.ngpus = ngpus
        self.mem = mem
        self.walltime = walltime
        self.mail_addr = mail_addr
        self.log_o = log_o
        self.log_e = log_e
    
    # return a list of header lines
    def get_header(self):
        header = []
        header.append(f"#!/bin/bash")
        header.append(f"#PBS -N {self.job_name}")
        header.append(f"#PBS -o {self.log_o}")
        header.append(f"#PBS -e {self.log_e}")
        header.append(f"#PBS -j oe")
        header.append(f"#PBS -m abe")
        header.append(f"#PBS -M {self.mail_addr}")
        header.append(f"#PBS -l ncpus={self.ncpus}")
        if self.ngpus>0:
            header.append(f"#PBS -l ngpus={self.ngpus}")
        header.append(f"#PBS -l mem={self.mem}")
        header.append(f"#PBS -l walltime={self.walltime}")

        header = [x+"\n" for x in header]
        return header
        
class GadiHeader(PBSHeader):
    def __init__(self, job_name, ncpus, ngpus, mem, walltime, mail_addr, log_o, log_e, project, storage, node_type, jobfs) -> None:
        super().__init__(job_name, ncpus, ngpus, mem, walltime, mail_addr, log_o, log_e)
        self.project = project
        self.storage = storage
        self.node_type = node_type
        self.jobfs = jobfs

    # return a list of header lines
    def get_header(self):
        header = []
        header.append(f"#!/bin/bash")
        header.append(f"#PBS -N {self.job_name}")
        header.append(f"#PBS -o {self.log_o}")
        header.append(f"#PBS -e {self.log_e}")
        header.append(f"#PBS -j oe")
        header.append(f"#PBS -m abe")
        header.append(f"#PBS -M {self.mail_addr}")
        header.append(f"#PBS -l ncpus={self.ncpus}")
        if self.ngpus>0:
            header.append(f"#PBS -l ngpus={self.ngpus}")
        header.append(f"#PBS -l mem={self.mem}")
        header.append(f"#PBS -l walltime={self.walltime}")

        header.append(f"#PBS -P {self.project}")
        header.append(f"#PBS -l storage={self.storage}")
        if self.ngpus>0:
            header.append(f"#PBS -q dgxa100")
        else:
            header.append(f"#PBS -q {self.node_type}")
        header.append(f"#PBS -l jobfs={self.jobfs}")

        header = [x+"\n" for x in header]
        return header

class Job():
    def __init__(self, job_manager, job_header, commands) -> None:
        self.job_manager = job_manager
        self.job_header = job_header
        self.commands = commands
        self.scripts = self.job_header.get_header() + self.commands
        
    
    def save_job(self, job_dir):
        suffix = {"pbs":"pbs", "gadi":"pbs", "bash":"sh"}
        with open(os.path.join(job_dir, f"{self.job_header.job_name}.{suffix[self.job_manager]}"), 'w') as f:
            f.writelines(self.scripts)

    def preview(self):
        for line in self.scripts:
            print(line, end="")


def generate_CAT_commands(job_header, out_dir, file_list):
    commands = [
        f'threads={job_header.ncpus}',
        f'CAT_dbPath={envs.CAT_PACK_DB_PATH}', # parameter:
        '',
        "file_list=(",
        "{}".format('\n'.join(f'"{item}"' for item in file_list)),
        ")",
        'for file in ${file_list[@]};do',
        '\techo "$file"',
        '',
        "\tfileHeader=$(echo $file | awk -F/ '{print $(NF)}' | cut -d '.' -f1)", # fileHeader
        f'\tout_dir={out_dir}/$fileHeader/CAT_results',
        '',
        '\tif [ ! -d "$out_dir/" ]; then',
        '\t\tmkdir -p $out_dir/',
        '\tfi',
        '',
        '\tif [ -f "$out_dir/$fileHeader.nr.contig2classification.with_names.txt" ];then',
        '\t\techo "$fileHeader CAT has already finished. Continue to the next one."',
        '\t\tcontinue',
        '\tfi',
        '',
        f'    {envs.CAT_PACK_PATH} contigs -c $file \\',
        '                -d $CAT_dbPath/db \\',
        '                -t $CAT_dbPath/tax \\',
        '                -n $threads \\',
        '                --force \\',
        '                -o $out_dir/$fileHeader.nr',
        '                ',
        f'    {envs.CAT_PACK_PATH} add_names -i $out_dir/$fileHeader.nr.contig2classification.txt \\',
        '                  -o $out_dir/$fileHeader.nr.contig2classification.with_names.txt \\',
        '                  -t $CAT_dbPath/tax \\',
        '                  --force \\',
        '                  --only_official ',
        '                  ',
        f'    {envs.CAT_PACK_PATH} summarise -c $file \\',
        '                  -i $out_dir/$fileHeader.nr.contig2classification.with_names.txt \\',
        '                  -o $out_dir/$fileHeader.nr.summary.txt',
        '',
        '    echo CAT on "$file" finished!',
        '    echo "############################################################" ',
        'done'
    ]
    commands = [x+"\n" for x in commands]
    return commands

def generate_VS2_commands(job_header, out_dir, file_list):
    commands = [
        f"source {envs.CONDA_PATH}/bin/activate vs2\n",
        f'threads={job_header.ncpus}',
        '',
        ' ',
        "file_list=(",
        "{}".format('\n'.join(f'"{item}"' for item in file_list)),
        ")",
        'for file in ${file_list[@]};do',
        '    echo "$file"',
        "    fileHeader=$(echo $file | awk -F/ '{print $(NF)}' | cut -d '.' -f1)", # fileHeader
        f'    out_dir={out_dir}/$fileHeader/VirSorter2_results', 
        '',
        '\tif [ ! -d "$out_dir/" ]; then',
        '\t\tmkdir -p $out_dir/',
        '\tfi',
        '',
        '\tif [ -f "$out_dir/$fileHeader-final-viral-score.tsv" ];then',
        '\t\techo "$fileHeader VirSorter2 has already finished. Continue to the next one."',
        '\t\tcontinue',
        '\tfi',
        '',
        f'    {envs.VIRSORTER2_PATH} run -w $out_dir -i $file -l $fileHeader --include-groups \"dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae\" -j $threads --rm-tmpdir all',
        '    echo VirSorter2 on "$file" finished!',
        '    echo "############################################################"',
        'done'
    ]
    commands = [x+"\n" for x in commands]
    return commands

def generate_GNM_commands(job_header, out_dir, file_list):
    commands = [
        f"source {envs.CONDA_PATH}/bin/activate genomad\n",
        f'threads={job_header.ncpus}',
        f'GeNomad_dbPath={envs.GENOMAD_DB_PATH}', # parameter:
        '',
        # 'source /g/data1b/oo46/wj6768/miniconda3/bin/activate /g/data1b/oo46/wj6768/miniconda3/envs/genomad',
        ' ',
        "file_list=(",
        "{}".format('\n'.join(f'"{item}"' for item in file_list)),
        ")",
        'for file in ${file_list[@]};do',
        '    echo "$file"',
        "    fileHeader=$(echo $file | awk -F/ '{print $(NF)}' | cut -d '.' -f1)", # fileHeader
        f'    out_dir={out_dir}/$fileHeader/GeNomad_results', 
        '',
        '\tif [ ! -d "$out_dir/" ]; then',
        '\t\tmkdir -p $out_dir/',
        '\tfi',
        "\theader=$(echo $file | awk -F \'/\' \'{print $NF}\')",
        '\tif [ -f "$out_dir/${header%.*}_summary/${header%.*}_virus_summary.tsv" ];then',
        '\t\tif [ -f "$out_dir/${header%.*}_summary/${header%.*}_plasmid_summary.tsv" ] ;then',
        '\t\t\techo "$fileHeader geNomad has already finished. Continue to the next one."',
        '\t\t\tcontinue',
        '\t\tfi',
        '\tfi',
        '',
        f'    {envs.GENOMAD_PATH} end-to-end $file $out_dir/ $GeNomad_dbPath --restart --threads $threads --cleanup',
        '    echo GeNomad on "$file" finished!',
        '    echo "############################################################" ',
        'done'
    ]
    commands = [x+"\n" for x in commands]
    return commands

def generate_VLM_commands(job_header, out_dir, file_list):
    commands = [
        f"source {envs.CONDA_PATH}/bin/activate viralm\n",
        f'threads={job_header.ncpus}',
        f'ViraLMPath={os.path.dirname(envs.VIRALM_PATH)}',
        '',
        # 'source /g/data1b/oo46/wj6768/miniconda3/bin/activate /g/data1b/oo46/wj6768/miniconda3/envs/viralm',
        ' ',
        "file_list=(",
        "{}".format('\n'.join(f'"{item}"' for item in file_list)),
        ")",
        'for file in ${file_list[@]};do',
        '    echo "$file"',
        "    fileHeader=$(echo $file | awk -F/ '{print $(NF)}' | cut -d '.' -f1)", # fileHeader
        f'    out_dir={out_dir}/$fileHeader/ViraLM_results',
        '',
        '\tif [ ! -d "$out_dir/" ]; then',
        '\t\tmkdir -p $out_dir/',
        '\tfi',
        '',
        '\tif [ -f "$out_dir/result_final.csv" ];then',
        '\t\techo "$fileHeader ViraLM has already finished. Continue to the next one."',
        '\t\tcontinue',
        '\telse',
        '\t\trm -r $out_dir',
        '\t\techo "$fileHeader ViraLM rerunning."',
        '\tfi',
        '',
        '    cd $ViraLMPath',
        f'    python {envs.VIRALM_PATH} --input $file --output $out_dir/',
        '    echo ViraLM on "$file" finished!',
        '    echo "############################################################" ',
        'done'
    ]
    commands = [x+"\n" for x in commands]
    return commands

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
    if config["job_manager"]=="bash":
        for chunk in data_chunks:
            # CAT
            cat_job_header = BashHeader(
                job_name=f"CAT_{chunk.index.to_list()[0]+1}_{chunk.index.to_list()[-1]+1}",
                ncpus=config["pbs"]["ncpus"],
            )
            cat_job = Job(
                job_manager=config["job_manager"],
                job_header=cat_job_header,
                commands=generate_CAT_commands(
                    job_header=cat_job_header,
                    out_dir=os.path.join(project_dir,"out"),
                    file_list=chunk["path"].to_list()
                ),
            )
            cat_job.save_job(job_dir=os.path.join(project_dir,"jobs"))
            # VS2
            vs2_job_header = BashHeader(
                job_name=f"VS2_{chunk.index.to_list()[0]+1}_{chunk.index.to_list()[-1]+1}",
                ncpus=config["pbs"]["ncpus"],
            )
            vs2_job = Job(
                job_manager=config["job_manager"],
                job_header=vs2_job_header,
                commands=generate_VS2_commands(
                    job_header=vs2_job_header,
                    out_dir=os.path.join(project_dir,"out"),
                    file_list=chunk["path"].to_list()
                ),
            )
            vs2_job.save_job(job_dir=os.path.join(project_dir,"jobs"))
            # GNM
            gnm_job_header = BashHeader(
                job_name=f"GNM_{chunk.index.to_list()[0]+1}_{chunk.index.to_list()[-1]+1}",
                ncpus=config["pbs"]["ncpus"],
            )
            gnm_job = Job(
                job_manager=config["job_manager"],
                job_header=gnm_job_header,
                commands=generate_GNM_commands(
                    job_header=gnm_job_header,
                    out_dir=os.path.join(project_dir,"out"),
                    file_list=chunk["path"].to_list()
                ),
            )
            gnm_job.save_job(job_dir=os.path.join(project_dir,"jobs"))
            # VLM
            vlm_job_header = PBSHeader(
                job_name=f"VLM_{chunk.index.to_list()[0]+1}_{chunk.index.to_list()[-1]+1}",
                ncpus=16,
            )
            vlm_job = Job(
                job_manager=config["job_manager"],
                job_header=vlm_job_header,
                commands=generate_VLM_commands(
                    job_header=vlm_job_header,
                    out_dir=os.path.join(project_dir,"out"),
                    file_list=chunk["path"].to_list()
                ),
            )
            vlm_job.save_job(job_dir=os.path.join(project_dir,"jobs"))
    if config["job_manager"]=="pbs":
        for chunk in data_chunks:
            # CAT
            cat_job_header = PBSHeader(
                job_name=f"CAT_{chunk.index.to_list()[0]+1}_{chunk.index.to_list()[-1]+1}",
                ncpus=config["pbs"]["ncpus"],
                ngpus=0,
                mem="192GB",
                walltime="48:00:00",
                mail_addr=config["pbs"]["mail_addr"],
                log_o=os.path.join(project_dir,"logs",f"CAT_{chunk.index.to_list()[0]+1}_{chunk.index.to_list()[-1]+1}.o"),
                log_e=os.path.join(project_dir,"logs",f"CAT_{chunk.index.to_list()[0]+1}_{chunk.index.to_list()[-1]+1}.e"),
            )
            cat_job = Job(
                job_manager=config["job_manager"],
                job_header=cat_job_header,
                commands=generate_CAT_commands(
                    job_header=cat_job_header,
                    out_dir=os.path.join(project_dir,"out"),
                    file_list=chunk["path"].to_list()
                ),
            )
            cat_job.save_job(job_dir=os.path.join(project_dir,"jobs"))
            # VS2
            vs2_job_header = PBSHeader(
                job_name=f"VS2_{chunk.index.to_list()[0]+1}_{chunk.index.to_list()[-1]+1}",
                ncpus=config["pbs"]["ncpus"],
                ngpus=0,
                mem=f"192GB",
                walltime="48:00:00",
                mail_addr=config["pbs"]["mail_addr"],
                log_o=os.path.join(project_dir,"logs",f"VS2_{chunk.index.to_list()[0]+1}_{chunk.index.to_list()[-1]+1}.o"),
                log_e=os.path.join(project_dir,"logs",f"VS2_{chunk.index.to_list()[0]+1}_{chunk.index.to_list()[-1]+1}.e"),
            )
            vs2_job = Job(
                job_manager=config["job_manager"],
                job_header=vs2_job_header,
                commands=generate_VS2_commands(
                    job_header=vs2_job_header,
                    out_dir=os.path.join(project_dir,"out"),
                    file_list=chunk["path"].to_list()
                ),
            )
            vs2_job.save_job(job_dir=os.path.join(project_dir,"jobs"))
            # GNM
            gnm_job_header = PBSHeader(
                job_name=f"GNM_{chunk.index.to_list()[0]+1}_{chunk.index.to_list()[-1]+1}",
                ncpus=config["pbs"]["ncpus"],
                ngpus=0,
                mem="192GB",
                walltime="48:00:00",
                mail_addr=config["pbs"]["mail_addr"],
                log_o=os.path.join(project_dir,"logs",f"GNM_{chunk.index.to_list()[0]+1}_{chunk.index.to_list()[-1]+1}.o"),
                log_e=os.path.join(project_dir,"logs",f"GNM_{chunk.index.to_list()[0]+1}_{chunk.index.to_list()[-1]+1}.e"),
            )
            gnm_job = Job(
                job_manager=config["job_manager"],
                job_header=gnm_job_header,
                commands=generate_GNM_commands(
                    job_header=gnm_job_header,
                    out_dir=os.path.join(project_dir,"out"),
                    file_list=chunk["path"].to_list()
                ),
            )
            gnm_job.save_job(job_dir=os.path.join(project_dir,"jobs"))
            # VLM
            vlm_job_header = PBSHeader(
                job_name=f"VLM_{chunk.index.to_list()[0]+1}_{chunk.index.to_list()[-1]+1}",
                ncpus=16,
                ngpus=1,
                mem="64GB",
                walltime="12:00:00",
                mail_addr=config["pbs"]["mail_addr"],
                log_o=os.path.join(project_dir,"logs",f"VLM_{chunk.index.to_list()[0]+1}_{chunk.index.to_list()[-1]+1}.o"),
                log_e=os.path.join(project_dir,"logs",f"VLM_{chunk.index.to_list()[0]+1}_{chunk.index.to_list()[-1]+1}.e"),
            )
            vlm_job = Job(
                job_manager=config["job_manager"],
                job_header=vlm_job_header,
                commands=generate_VLM_commands(
                    job_header=vlm_job_header,
                    out_dir=os.path.join(project_dir,"out"),
                    file_list=chunk["path"].to_list()
                ),
            )
            vlm_job.save_job(job_dir=os.path.join(project_dir,"jobs"))
    if config["job_manager"]=="gadi":
        for chunk in data_chunks:
            # CAT
            cat_job_header = GadiHeader(
                job_name=f"CAT_{chunk.index.to_list()[0]+1}_{chunk.index.to_list()[-1]+1}",
                ncpus=config["pbs"]["ncpus"],
                ngpus=0,
                mem="192GB",
                walltime="48:00:00",
                mail_addr=config["pbs"]["mail_addr"],
                log_o=os.path.join(project_dir,"logs",f"CAT_{chunk.index.to_list()[0]+1}_{chunk.index.to_list()[-1]+1}.o"),
                log_e=os.path.join(project_dir,"logs",f"CAT_{chunk.index.to_list()[0]+1}_{chunk.index.to_list()[-1]+1}.e"),
                project=config["pbs"]["gadi"]["-P project"],
                storage=config["pbs"]["gadi"]["-l storage"],
                node_type="normalsl",
                jobfs="2GB",
            )
            cat_job = Job(
                job_manager=config["job_manager"],
                job_header=cat_job_header,
                commands=generate_CAT_commands(
                    job_header=cat_job_header,
                    out_dir=os.path.join(project_dir,"out"),
                    file_list=chunk["path"].to_list()
                ),
            )
            cat_job.save_job(job_dir=os.path.join(project_dir,"jobs"))
            # VS2
            vs2_job_header = GadiHeader(
                job_name=f"VS2_{chunk.index.to_list()[0]+1}_{chunk.index.to_list()[-1]+1}",
                ncpus=config["pbs"]["ncpus"],
                ngpus=0,
                mem=f"192GB",
                walltime="48:00:00",
                mail_addr=config["pbs"]["mail_addr"],
                log_o=os.path.join(project_dir,"logs",f"VS2_{chunk.index.to_list()[0]+1}_{chunk.index.to_list()[-1]+1}.o"),
                log_e=os.path.join(project_dir,"logs",f"VS2_{chunk.index.to_list()[0]+1}_{chunk.index.to_list()[-1]+1}.e"),
                project=config["pbs"]["gadi"]["-P project"],
                storage=config["pbs"]["gadi"]["-l storage"],
                node_type="normalsl",
                jobfs="2GB",
            )
            vs2_job = Job(
                job_manager=config["job_manager"],
                job_header=vs2_job_header,
                commands=generate_VS2_commands(
                    job_header=vs2_job_header,
                    out_dir=os.path.join(project_dir,"out"),
                    file_list=chunk["path"].to_list()
                ),
            )
            vs2_job.save_job(job_dir=os.path.join(project_dir,"jobs"))
            # GNM
            gnm_job_header = GadiHeader(
                job_name=f"GNM_{chunk.index.to_list()[0]+1}_{chunk.index.to_list()[-1]+1}",
                ncpus=config["pbs"]["ncpus"],
                ngpus=0,
                mem="192GB",
                walltime="48:00:00",
                mail_addr=config["pbs"]["mail_addr"],
                log_o=os.path.join(project_dir,"logs",f"GNM_{chunk.index.to_list()[0]+1}_{chunk.index.to_list()[-1]+1}.o"),
                log_e=os.path.join(project_dir,"logs",f"GNM_{chunk.index.to_list()[0]+1}_{chunk.index.to_list()[-1]+1}.e"),
                project=config["pbs"]["gadi"]["-P project"],
                storage=config["pbs"]["gadi"]["-l storage"],
                node_type="normalsl",
                jobfs="2GB",
            )
            gnm_job = Job(
                job_manager=config["job_manager"],
                job_header=gnm_job_header,
                commands=generate_GNM_commands(
                    job_header=gnm_job_header,
                    out_dir=os.path.join(project_dir,"out"),
                    file_list=chunk["path"].to_list()
                ),
            )
            gnm_job.save_job(job_dir=os.path.join(project_dir,"jobs"))
            # VLM
            vlm_job_header = GadiHeader(
                job_name=f"VLM_{chunk.index.to_list()[0]+1}_{chunk.index.to_list()[-1]+1}",
                ncpus=16,
                ngpus=1,
                mem="64GB",
                walltime="12:00:00",
                mail_addr=config["pbs"]["mail_addr"],
                log_o=os.path.join(project_dir,"logs",f"VLM_{chunk.index.to_list()[0]+1}_{chunk.index.to_list()[-1]+1}.o"),
                log_e=os.path.join(project_dir,"logs",f"VLM_{chunk.index.to_list()[0]+1}_{chunk.index.to_list()[-1]+1}.e"),
                project=config["pbs"]["gadi"]["-P project"],
                storage=config["pbs"]["gadi"]["-l storage"],
                node_type="dgxa100",
                jobfs="2GB",
            )
            vlm_job = Job(
                job_manager=config["job_manager"],
                job_header=vlm_job_header,
                commands=generate_VLM_commands(
                    job_header=vlm_job_header,
                    out_dir=os.path.join(project_dir,"out"),
                    file_list=chunk["path"].to_list()
                ),
            )
            vlm_job.save_job(job_dir=os.path.join(project_dir,"jobs"))

if __name__=="__main__":
    pass
    # PBS_header = PBSHeader(
    #     job_name="test",
    #     ncpus=32,
    #     ngpus=0,
    #     mem="64GB",
    #     walltime="1:00:00",
    #     mail_addr="test@test.com",
    #     log_o="path",
    #     log_e="path",
    # )
    # Gadi_header = GadiHeader(
    #     job_name="test",
    #     ncpus=32,
    #     ngpus=0,
    #     mem="64GB",
    #     walltime="1:00:00",
    #     mail_addr="test@test.com",
    #     log_o="path",
    #     log_e="path",
    #     project="mp96",
    #     storage="gdata/oo46+gdata/mp96",
    #     node_type="normalsl",
    #     jobfs="2GB",
    # )

    # commands = [x+"\n" for x in ["cd ./", "echo \"hello world!\""]]
    # test_job = Job(job_header=BashHeader(job_name="test"), commands=commands)
    # test_job.save_job(path=os.path.join(os.getcwd(), f"{test_job.job_header.job_name}.pbs"))
 