import pandas as pd
import numpy as np
import os
from Bio import SeqIO
import logging
logging.basicConfig(
    level=logging.DEBUG, 
    format="%(asctime)s [%(levelname)s]: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
import time

import envs
from utils import job_management
from functools import reduce

def extract_putative_contigs_single_sample(prj_dir, fileHeader, fasta_path, min_len=3000, num_tools=2, trusted=None, skip=None):
    '''
    extract putative contigs from the results of:
    {
        cat: CAT,
        vs2: VirSorter2, 
        gnm: GeNoma, 
        vlm: ViraLM
    }
    Prameters:
        Positional:
            prj_dir: str, the directory of the project.
            fileHeader: str, the header of the input file.
            fasta_path: str, the path of the input file.
        Optional:
            min_len: int, the minimum length of the contig to be considered as putative. Default: 3000.
            num_tools: int, a contig must be classified as viral by at least this many tools to be considered putative. Default: 2.
            trusted: str, the trusted tool(s) to be used for classification, selected from {cat,vs2,gnm,vlm}, comma seperated. Default: 'cat'.
            skip: str, a comma-separated list of tools to skip in the results, selected from {cat,vs2,gnm,vlm}, comma seperated. Default: None.
    Returns:
        
    '''
    # define functions to process results from different tools
    def process_cat():
        file_path = os.path.join(prj_dir,"out",f"{fileHeader}","CAT_results",f"{fileHeader}.nr.contig2classification.with_names.txt")
        if not os.path.isfile(file_path):
            logging.warning(f"{file_path} does not exist. Skipped CAT_pack.")
            return None
        cat = pd.read_table(file_path, sep='\t', header=0).rename({"# contig":"contig"},axis=1)
        cat = cat[cat["classification"]=="taxid assigned"]
        cat.columns = ["cat_" + x for x in cat.columns.values.tolist()]
        cat = cat.rename({"cat_contig":"seq_name", "cat_superkingdom":"cat_category"}, axis=1)
        cat["cat_category"] = cat["cat_category"].str.split(":", expand=True)[0]
        cat = cat.loc[:,["seq_name","cat_category"]].reset_index(drop=True)
        return cat
    def process_vs2():
        file_path = os.path.join(prj_dir,"out",f"{fileHeader}","VirSorter2_results", f"{fileHeader}-final-viral-score.tsv")
        if not os.path.isfile(file_path):
            logging.warning(f"{file_path} does not exist. Skipped VirSorter2.")
            return None
        vs2 = pd.read_table(file_path, sep='\t', header=0)
        if vs2.shape[0]>0:
            vs2["completeness"] = vs2["seqname"].str.split("\|\|", expand=True)[1]
            vs2["seqname"] = vs2["seqname"].str.split("\|\|", expand=True)[0]
            vs2["category"] = "Viruses"
            vs2.columns = ["vs2_" + x for x in vs2.columns.values.tolist()]
            
            vs2 = vs2.rename({"vs2_seqname":"seq_name"}, axis=1)
            vs2 = vs2[vs2["vs2_max_score_group"]!="RNA"]
            vs2 = vs2.loc[:, ["seq_name", "vs2_category"]].reset_index(drop=True)
        else:
            vs2["completeness"] = None
            vs2["seqname"] = None
            vs2["category"] = None
            vs2.columns = ["vs2_" + x for x in vs2.columns.values.tolist()]
            vs2 = vs2.rename({"vs2_seqname":"seq_name"}, axis=1)
            vs2 = vs2[vs2["vs2_max_score_group"]!="RNA"]
            vs2 = vs2.loc[:, ["seq_name", "vs2_category"]].reset_index(drop=True)
        return vs2
    def process_gnm():
        gnm_v_path = os.path.join(
            prj_dir,"out",f"{fileHeader}","GeNomad_results",
            f"{'.'.join(fasta_path.split('/')[-1].split('.')[:-1])}_summary",
            f"{'.'.join(fasta_path.split('/')[-1].split('.')[:-1])}_virus_summary.tsv"
        )
        gnm_p_path = os.path.join(
            prj_dir,"out",f"{fileHeader}","GeNomad_results",
            f"{'.'.join(fasta_path.split('/')[-1].split('.')[:-1])}_summary",
            f"{'.'.join(fasta_path.split('/')[-1].split('.')[:-1])}_plasmid_summary.tsv"
        )
        if not os.path.isfile(gnm_v_path) or not os.path.isfile(gnm_p_path):
            logging.warning(f"{gnm_v_path} or {gnm_p_path} does not exist. Skipped GeNomad.")
            return None
        gnm_v = pd.read_table(gnm_v_path, sep='\t', header=0)
        gnm_p = pd.read_table( gnm_p_path, sep='\t', header=0)
        gnm_v["category"] = "Viruses"
        gnm_p["category"] = "Plasmids"
        gnm = pd.merge(gnm_v, gnm_p, 'outer')
        gnm.columns = ["gnm_" + x for x in gnm.columns.values.tolist()]
        gnm = gnm.rename({"gnm_seq_name":"seq_name"}, axis=1)
        gnm = gnm.loc[:,["seq_name", "gnm_category"]].reset_index(drop=True)
        return gnm
    def process_vlm():
        file_path = os.path.join(prj_dir,"out",f"{fileHeader}","ViraLM_results",f"result_{fileHeader}.csv")
        if not os.path.isfile(file_path):
            logging.warning(f"{file_path} does not exist. Skipped ViraLM.")
            return None
        vlm = pd.read_table(file_path, sep=',', header=0)
        vlm = vlm[vlm["virus_score"]>=0.8]
        vlm["category"] = "Viruses"
        vlm.columns = ["vlm_" + x for x in vlm.columns.values.tolist()]
        vlm = vlm.rename({"vlm_seq_name":"seq_name"}, axis=1)
        vlm = vlm.loc[:, ["seq_name", "vlm_category"]].reset_index(drop=True)
        return vlm
    # define a function to acquire length info from fasta file
    def get_sequence_lengths(fasta_file):
        """
        Read FASTA file, return a dictionary with sequence names and their lengths.
        
        Parameters:
            fasta_file (str): FASTA filepath.
        
        return:
            dict: {seq_name: seq_length}
        """
        seq_lengths = {}
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq_lengths[record.id] = len(record.seq)
        df = pd.DataFrame.from_dict(seq_lengths, orient='index', columns=['length']).reset_index().rename({'index':'seq_name'}, axis=1)
        return df
    # define a function to filter putative
    def filter_putative(all, min_len, num_tools, trusted):
        '''
        Filter putative contigs based on length, number of tools, and trusted tools.
        '''
        # select entries with length >= min_len
        putative_min_len = all[all["length"]>=min_len].reset_index(drop=True).set_index(["seq_name","length"])
        # select entries with v_count >= num_tools
        putative_min_len["v_count"] = putative_min_len.apply(lambda x: x[x=='Viruses'].count(), axis=1)
        putative_min_len = putative_min_len[(putative_min_len["v_count"]>=num_tools) & (putative_min_len["gnm_category"]!="Plasmids")].reset_index().astype({"length":int}).astype(str)
        # select entries with trusted tools
        cols = all.columns[all.columns.str.contains("_category")].tolist()
        trusted_cols = [col for col in cols if col.startswith(tuple(trusted.split(',')))] if trusted is not None else []
        putative_trusted = all[all[trusted_cols].eq("Viruses").any(axis=1)].reset_index(drop=True)
        # combine filtered and trusted as putative
        putative = all.set_index('seq_name').loc[list(set(putative_min_len['seq_name'].tolist()) | set(putative_trusted['seq_name'].tolist())),:].reset_index()
        # save putative summary
        putative.to_csv(os.path.join(prj_dir,"out",fileHeader,"putative_summary.csv"), index=None)
        return putative
    # define a function to extract sequences from fasta file
    def extract_sequences(fasta_file, seq_names, output_file=None):
        sequences = {}
        temp_records = []  # store records to write, temporarily
        with open(fasta_file, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                if record.id in seq_names:
                    sequences[record.id] = str(record.seq)
                    temp_records.append(record)  # collect records
        if output_file and temp_records:
            with open(output_file, "w") as out:  # overwrite the output file
                SeqIO.write(temp_records, out, "fasta")
        return sequences
    # list of results to check (subset of ["cat", "vs2", "gnm", "vlm"])
    default_set = {"cat", "vs2", "gnm", "vlm"}
    skip_set = set(skip.split(',')) if skip is not None else set()
    check_set = default_set - skip_set
    # iterate through the tools and process the results
    df_list = []
    for tool in check_set:
        if tool == "cat":
            df_list.append(process_cat())
        elif tool == "vs2":
            df_list.append(process_vs2())
        elif tool == "gnm":
            df_list.append(process_gnm())
        elif tool == "vlm":
            df_list.append(process_vlm())
    df_list = [df for df in df_list if df is not None]
    if df_list == []:
        logging.error("No results found. Please check the input files.")
        return
    for df in df_list:
        if 'seq_name' not in df.columns:
            raise ValueError("All DataFrame must contain 'seq_name' column.")
    # merge results
    all = reduce(lambda left, right: pd.merge(left, right, on='seq_name', how='outer'), df_list)
    # read completeness status
    completeness_status = pd.read_csv(os.path.join(prj_dir, "completeness_status.csv"), sep=',', header=0, index_col=None)
    # acquire length info from fasta file
    fasta_path = completeness_status[completeness_status['fileHeader']==fileHeader].iloc[0]['path']

    seq_length_info = get_sequence_lengths(fasta_path)
    # merge length info into summary
    all = pd.merge(all, seq_length_info, on="seq_name", how='left')
    # filter putative
    putative = filter_putative(all=all, min_len=min_len, num_tools=num_tools, trusted=trusted)
    # extract fasta
    putative_seq_names = putative["seq_name"].tolist()
    putative_fasta_path = os.path.join(prj_dir, 'out', fileHeader,'putative_contigs.fasta')

    start_time = time.time()
    putative_sequences = extract_sequences(fasta_file=fasta_path, seq_names=putative_seq_names, output_file=putative_fasta_path)
    elapsed_time = time.time() - start_time
    logging.info(f"Extracted {len(putative_sequences)} putative sequences to {putative_fasta_path} in {elapsed_time:.2f} seconds.")

    return

def find_rRNAs_single_file(prj_dir, fileHeader, threads=32):
    
    # pre-run check
    files_to_check = [
        os.path.join(prj_dir, 'out', fileHeader, 'putative_contigs.fasta')
    ]
    for file in files_to_check:
        if not os.path.exists(file):
            print(f"{file} does not exist, skipped.")
            return
    
    if os.path.exists(os.path.join(prj_dir, "out", fileHeader, "rRNAs.tsv")):
        os.remove(os.path.join(prj_dir, "out", fileHeader, "rRNAs.tsv"))
    bash_commands = [
        f"source /g/data1b/oo46/wj6768/miniconda3/bin/activate /g/data1b/oo46/wj6768/miniconda3/envs/mpa\n",
        f"barrnap --kingdom bac --threads {threads} {os.path.join(prj_dir, 'out', fileHeader, 'putative_contigs.fasta')} | sed '1d' >> {os.path.join(prj_dir, 'out', fileHeader, 'rRNAs.tsv')}\n",
        f"barrnap --kingdom arc --threads {threads} {os.path.join(prj_dir, 'out', fileHeader, 'putative_contigs.fasta')} | sed '1d' >> {os.path.join(prj_dir, 'out', fileHeader, 'rRNAs.tsv')}\n",
        f"barrnap --kingdom euk --threads {threads} {os.path.join(prj_dir, 'out', fileHeader, 'putative_contigs.fasta')} | sed '1d' >> {os.path.join(prj_dir, 'out', fileHeader, 'rRNAs.tsv')}\n",
        f"barrnap --kingdom mito --threads {threads} {os.path.join(prj_dir, 'out', fileHeader, 'putative_contigs.fasta')} | sed '1d' >> {os.path.join(prj_dir, 'out', fileHeader, 'rRNAs.tsv')}\n"
    ]
    with open(os.path.join(prj_dir, "find_rRNAs_tmp.sh"), 'w') as f:
        f.writelines("#!/bin/bash\n")
        f.writelines(bash_commands)
    os.system(f"chmod +x {os.path.join(prj_dir, 'find_rRNAs_tmp.sh')}")
    os.system(os.path.join(prj_dir, "find_rRNAs_tmp.sh"))
    os.remove(os.path.join(prj_dir, "find_rRNAs_tmp.sh"))
    if os.path.exists(f"{os.path.join(prj_dir, 'out', fileHeader, 'putative_contigs.fasta.fai')}"):
        os.remove(f"{os.path.join(prj_dir, 'out', fileHeader, 'putative_contigs.fasta.fai')}")
    
def extract_decontaminated_contigs_single_file(prj_dir, fileHeader):
    
    # pre-run check
    files_to_check = [
        os.path.join(prj_dir, 'out', fileHeader, 'putative_summary.csv'),
        os.path.join(prj_dir, 'out', fileHeader, 'rRNAs.tsv'),
        os.path.join(prj_dir, "completeness_status.csv")
    ]
    for file in files_to_check:
        if not os.path.exists(file):
            return
    
    putative_summary = pd.read_table(os.path.join(prj_dir, 'out', fileHeader, 'putative_summary.csv'), sep=',', header=0).astype({"length":int})
    rRNAs_summary = pd.read_table(os.path.join(prj_dir, 'out', fileHeader, 'rRNAs.tsv'), header=None, names=["seq_name","source","type","start","end","score","strand","phase","attributes"])
    completeness_status = pd.read_csv(os.path.join(prj_dir, "completeness_status.csv"), sep=',', header=0, index_col=None)
    
    confirmed_summary = putative_summary.copy()
    confirmed_summary = confirmed_summary[~confirmed_summary["seq_name"].isin(rRNAs_summary["seq_name"])]
    confirmed_summary.to_csv(os.path.join(prj_dir, 'out', fileHeader, 'decontaminated_summary.csv'), sep=',', index=None)

    # extract confirmed fasta
    bash_commands = [
        f"seqkit grep -f <(sed '1d' {os.path.join(prj_dir, 'out', fileHeader, 'decontaminated_summary.csv')} | cut -f1 -d',') {completeness_status[completeness_status['fileHeader']==fileHeader].iloc[0]['path']} > {os.path.join(prj_dir, 'out', fileHeader, 'decontaminated_contigs.fasta')}\n"
    ]
    with open(os.path.join(prj_dir, "extract_decontaminated_contigs_tmp.sh"), 'w') as f:
        f.writelines("#!/bin/bash\n")
        f.writelines(bash_commands)
    os.system(f"chmod +x {os.path.join(prj_dir, 'extract_decontaminated_contigs_tmp.sh')}")
    os.system(os.path.join(prj_dir, "extract_decontaminated_contigs_tmp.sh"))
    os.remove(os.path.join(prj_dir, "extract_decontaminated_contigs_tmp.sh"))

def extract_putative_contigs_multi_samples(prj_dir, min_len=3000, num_tools=2, trusted=None, skip=None):
    '''
    multi-sample implementation of extract_putative_contigs_single_sample(), which:
        extracts putative contigs from the results of:
            {
                cat: CAT,
                vs2: VirSorter2, 
                gnm: GeNoma, 
                vlm: ViraLM
            }
    Prameters:
        Positional:
            prj_dir: str, the directory of the project.
        Optional:
            min_len: int, the minimum length of the contig to be considered as putative. Default: 3000.
            num_tools: int, a contig must be classified as viral by at least this many tools to be considered putative. Default: 2.
            trusted: str, the trusted tool(s) to be used for classification, selected from {cat,vs2,gnm,vlm}, comma seperated. Default: 'cat'.
            skip: str, a list of tools to skip in the results, selected from {cat,vs2,gnm,vlm}, comma seperated. Default: None.
    Returns:
        
    '''
    fileHeader_list = pd.read_csv(os.path.join(prj_dir,"completeness_status.csv"),sep=',',header=0,index_col=None).loc[:,"fileHeader"].tolist()
    status = pd.read_csv(os.path.join(prj_dir,"completeness_status.csv"),sep=',',header=0,index_col=None)
    for fileHeader in fileHeader_list:
        # if os.path.isfile(os.path.join(prj_dir,"out",fileHeader,"putative_contigs.fasta")):
        #     print(f"{fileHeader} has finished, skip.")
        #     continue
        # else:
        extract_putative_contigs_single_sample(
            prj_dir=prj_dir, fileHeader=fileHeader, fasta_path=status.set_index("fileHeader").loc[fileHeader,"path"], 
            min_len=min_len, num_tools=num_tools,
            trusted=trusted, skip=skip
        )
        # break
def find_rRNAs_multi_files(prj_dir, threads=32):
    fileHeader_list = pd.read_csv(os.path.join(prj_dir,"completeness_status.csv"),sep=',',header=0,index_col=None).loc[:,"fileHeader"].tolist()
    for fileHeader in fileHeader_list:
        # if os.path.isfile(os.path.join(prj_dir,"out",fileHeader,"rRNAs.tsv")):
        #     print(f"{fileHeader} has finished, skip.")
        #     continue
        # else:
        find_rRNAs_single_file(prj_dir=prj_dir, fileHeader=fileHeader, threads=threads)
def extract_decontaminated_contigs_multi_files(prj_dir):
    fileHeader_list = pd.read_csv(os.path.join(prj_dir,"completeness_status.csv"),sep=',',header=0,index_col=None).loc[:,"fileHeader"].tolist()
    for fileHeader in fileHeader_list:
        # if os.path.isfile(os.path.join(prj_dir,"out",fileHeader,"decontaminated_contigs.fasta")):
        #     print(f"{fileHeader} has finished, skip.")
        #     continue
        # else:
        extract_decontaminated_contigs_single_file(prj_dir=prj_dir, fileHeader=fileHeader)

def merge_confirmed_contigs(prj_dir):
    fileHeader_list = pd.read_csv(os.path.join(prj_dir,"completeness_status.csv"),sep=',',header=0,index_col=None).loc[:,"fileHeader"].tolist()
    if not os.path.isdir(os.path.join(prj_dir,"OVU")):
        os.makedirs(os.path.join(prj_dir,"OVU"), exist_ok=True)
    with open(os.path.join(prj_dir,"OVU","merged_decontaminated_contigs.fasta"), 'w') as merged_confirmed_contigs:
        for fileHeader in fileHeader_list:
            if not os.path.exists(os.path.join(prj_dir, 'out', fileHeader, 'decontaminated_contigs.fasta')):
                print(f"{os.path.join(prj_dir, 'out', fileHeader, 'decontaminated_contigs.fasta')} not exist. Exiting.")
                return
        for fileHeader in fileHeader_list:
            with open(os.path.join(prj_dir, 'out', fileHeader, 'decontaminated_contigs.fasta'), 'r') as fasta:
                sequence = []
                for line in fasta:
                    if line.startswith('>'):
                        sequence.append(f">{fileHeader}_{line[1:]}")
                    else:
                        sequence.append(line)
                merged_confirmed_contigs.writelines(sequence)

def dedup(prj_dir):
    merged_fasta = os.path.join(prj_dir,"OVU","merged_decontaminated_contigs.fasta")
    dedup_details = os.path.join(prj_dir,"OVU","merged_decontaminated_contigs_dedup_detail.txt")
    dup = os.path.join(prj_dir,"OVU","merged_decontaminated_contigs_dup.fasta")
    dedup = os.path.join(prj_dir,"OVU","merged_decontaminated_contigs_dedup.fasta")
    if not os.path.isdir(os.path.join(prj_dir,"OVU")):
        os.makedirs(os.path.join(prj_dir,"OVU"), exist_ok=True)
    bash_commands = [
        f"source {envs.CONDA_PATH}/bin/activate {envs.MAIN_ENV_NAME}\n",
        f"cat {merged_fasta} | seqkit rmdup -s -D {dedup_details} -d {dup} -o {dedup}\n",
    ]
    with open(os.path.join(prj_dir, "dedup_tmp.sh"), 'w') as f:
        f.writelines("#!/bin/bash\n")
        f.writelines(bash_commands)
    os.system(f"chmod +x {os.path.join(prj_dir, 'dedup_tmp.sh')}")
    os.system(os.path.join(prj_dir, "dedup_tmp.sh"))
    os.remove(os.path.join(prj_dir, "dedup_tmp.sh"))

def check_quality(prj_dir, config):
    dedup = os.path.join(prj_dir,"OVU","merged_decontaminated_contigs_dedup.fasta")
    quality_check_dir = os.path.join(prj_dir,'OVU','quality_check')
    quality_filtered_fasta = os.path.join(prj_dir,'OVU','quality_filtered_viral_contigs.fasta')
    log_dir = os.path.join(prj_dir,'OVU')
    job_dir = os.path.join(prj_dir,'OVU')
    if config['job_manager'] in ['pbs', 'gadi', 'bash']:
        threads = config['ncpus']
    if not os.path.isdir(os.path.join(prj_dir,"OVU")):
        os.makedirs(os.path.join(prj_dir,"OVU"), exist_ok=True)
    bash_commands = [
        f"source {envs.CONDA_PATH}/bin/activate {envs.MAIN_ENV_NAME}",
        f"checkv end_to_end {dedup} {quality_check_dir} -d {envs.CHECKV_DB_PATH} -t {threads}",
        f"seqkit grep -f <(awk -F \'\\t\' \'{{if ($6 == 0 && $8 == \"Not-determined\") {{next;}} print $1}}\' {quality_check_dir}/quality_summary.tsv | sed \'1d\') {dedup} > {quality_filtered_fasta}",
    ]
    bash_commands = [x+"\n" for x in bash_commands]
    if config['job_manager']=='pbs':
        check_quality_job_header = job_management.PBSHeader(
            job_name="check_quality",
            ncpus=threads,
            ngpus=0,
            mem=f"{int(threads*4)}GB",
            walltime="10:00:00",
            mail_addr=config['pbs']['mail_addr'],
            log_o=f"{log_dir}/check_quality.o",
            log_e=f"{log_dir}/check_quality.e",
        )
        check_quality_job = job_management.Job(
            job_manager='pbs',job_header=check_quality_job_header, commands=bash_commands,
        )
        check_quality_job.save_job(job_dir=job_dir)
    elif config['job_manager']=='gadi':
        check_quality_job_header = job_management.GadiHeader(
            job_name="check_quality",
            ncpus=threads,
            ngpus=0,
            mem=f"{int(threads*4)}GB",
            walltime="10:00:00",
            mail_addr=config['pbs']['mail_addr'],
            log_o=f"{log_dir}/check_quality.o",
            log_e=f"{log_dir}/check_quality.e",
            project=config['pbs']['gadi']['-P project'],
            storage=config['pbs']['gadi']['-l storage'],
            node_type="normalsl",
            jobfs="2GB",
        )
        check_quality_job = job_management.Job(
            job_manager='gadi',job_header=check_quality_job_header, commands=bash_commands,
        )
        check_quality_job.save_job(job_dir=job_dir)
    elif config['job_manager']=='bash':
        check_quality_job_header = job_management.BashHeader(
            job_name="check_quality",
            ncpus=threads
        )
        check_quality_job = job_management.Job(
            job_manager='bash',job_header=check_quality_job_header, commands=bash_commands,
        )
        check_quality_job.save_job(job_dir=job_dir)


def cluster(prj_dir, config):
    job_dir = os.path.join(prj_dir,"OVU")
    log_dir = os.path.join(prj_dir,"OVU")
    options = {
        "checked": os.path.join(prj_dir,'OVU','quality_filtered_viral_contigs.fasta'),
        "unchecked": os.path.join(prj_dir,"OVU","merged_decontaminated_contigs_dedup.fasta"),
    }
    fasta_file = os.path.join(prj_dir,'OVU','quality_filtered_viral_contigs.fasta')
    script_path = os.path.join(envs.INSTALLATION_PATH,"src")
    if config['job_manager'] in ['pbs', 'gadi', 'bash']:
        threads = config['ncpus']
    bash_commands = [
        f"source {envs.CONDA_PATH}/bin/activate {envs.MAIN_ENV_NAME}",
        f"echo \"make blast db ...\"",
        f"makeblastdb -in {fasta_file} -out {prj_dir}/OVU/blastdb_for_anicluster/blastdb_for_anicluster -dbtype nucl",
        f"echo \"blasting ...\"",
        f"blastn -query {fasta_file} -db {prj_dir}/OVU/blastdb_for_anicluster/blastdb_for_anicluster -out {prj_dir}/OVU/filtered_blast.tsv -outfmt '6 std qlen slen' -max_target_seqs 25000 -perc_identity 90",
        f"echo \"blast finished\"",
        f"python {os.path.join(script_path, 'blastani.py')} -i {prj_dir}/OVU/filtered_blast.tsv -o {prj_dir}/OVU/filtered_ani.tsv",
        f"echo \"compute ANI finished\"",
        f"python {os.path.join(script_path, 'cluster.py')} --fna {fasta_file} --ani {prj_dir}/OVU/filtered_ani.tsv --out {prj_dir}/OVU/filtered_clusters.tsv --min_ani 95 --min_qcov 0 --min_tcov 85",
        f"echo \"cluster finished\"",
        f"echo \"extract representatives ...\"",
        f"seqkit grep -f <(cat {prj_dir}/OVU/filtered_clusters.tsv | cut -f1) {fasta_file} > {prj_dir}/OVU/rep_contigs.fasta",
        f"echo \"finished\"",
    ]
    bash_commands = [x+"\n" for x in bash_commands]
    if config['job_manager']=='pbs':
        cluster_job_header = job_management.PBSHeader(
            job_name="cluster",
            ncpus=threads,
            ngpus=0,
            mem=f"{int(threads*4)}GB",
            walltime="10:00:00",
            mail_addr=config['pbs']['mail_addr'],
            log_o=f"{log_dir}/cluster.o",
            log_e=f"{log_dir}/cluster.e",
        )
        cluster_job = job_management.Job(
            job_manager='pbs',job_header=cluster_job_header, commands=bash_commands,
        )
        cluster_job.save_job(job_dir=job_dir)
    elif config['job_manager']=='gadi':
        cluster_job_header = job_management.GadiHeader(
            job_name="cluster",
            ncpus=threads,
            ngpus=0,
            mem=f"{int(threads*4)}GB",
            walltime="10:00:00",
            mail_addr=config['pbs']['mail_addr'],
            log_o=f"{log_dir}/cluster.o",
            log_e=f"{log_dir}/cluster.e",
            project=config['pbs']['gadi']['-P project'],
            storage=config['pbs']['gadi']['-l storage'],
            node_type="normalsl",
            jobfs="2GB",
        )
        cluster_job = job_management.Job(
            job_manager='gadi',job_header=cluster_job_header, commands=bash_commands,
        )
        cluster_job.save_job(job_dir=job_dir)
    elif config['job_manager']=='bash':
        cluster_job_header = job_management.BashHeader(
            job_name="cluster",
            ncpus=threads
        )
        cluster_job = job_management.Job(
            job_manager='bash',job_header=cluster_job_header, commands=bash_commands,
        )
        cluster_job.save_job(job_dir=job_dir)
    
    return

# def find_genomad_header():
#     genomad_out_header = os.path.join("/g/data/oo46/wj6768/Healthy_Virome_HOAM_STOOL","out","HOAM22501", "GeNomad_results", "*_summary", "*_virus_summary.tsv")
#     # genomad_out_header = os.path.basename(glob.glob("/g/data/oo46/wj6768/Healthy_Virome_HOAM_STOOL/out/HOAM22501/GeNomad_results/*_summary")[0]).split()
#     print(genomad_out_header)
#     return 

if __name__=="__main__":
    pass
    # find_genomad_header()