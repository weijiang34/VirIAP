import pandas as pd
import numpy as np
import os

def extract_putative_contigs_single_sample(prj_dir, fileHeader, min_len=3000):
    # pre-run check
    files_to_check = [
        os.path.join(prj_dir,"out",f"{fileHeader}","CAT_results", f"{fileHeader}.nr.contig2classification.with_names.txt"),
        os.path.join(prj_dir,"out",f"{fileHeader}","VirSorter2_results", f"{fileHeader}-final-viral-score.tsv"),
        os.path.join(prj_dir,"out",f"{fileHeader}", "GeNomad_results", "final.contigs_summary", "final.contigs_virus_summary.tsv"),
        os.path.join(prj_dir,"out",f"{fileHeader}","GeNomad_results","final.contigs_summary", "final.contigs_plasmid_summary.tsv"),
        os.path.join(prj_dir,"out",f"{fileHeader}","ViraLM_results","result_final.csv")
    ]
    for file in files_to_check:
        if not os.path.exists(file):
            return
    
    cat = pd.read_table(os.path.join(prj_dir,"out",f"{fileHeader}","CAT_results",
                                     f"{fileHeader}.nr.contig2classification.with_names.txt"), sep='\t', header=0).rename({"# contig":"contig"},axis=1)
    vs2 = pd.read_table(os.path.join(prj_dir,"out",f"{fileHeader}","VirSorter2_results",
                                     f"{fileHeader}-final-viral-score.tsv"), sep='\t', header=0)
    gnm_v = pd.read_table(os.path.join(prj_dir,"out",f"{fileHeader}","GeNomad_results","final.contigs_summary",
                                     "final.contigs_virus_summary.tsv"), sep='\t', header=0)
    gnm_p = pd.read_table(os.path.join(prj_dir,"out",f"{fileHeader}","GeNomad_results","final.contigs_summary",
                                     "final.contigs_plasmid_summary.tsv"), sep='\t', header=0)
    vlm = pd.read_table(os.path.join(prj_dir,"out",f"{fileHeader}","ViraLM_results","result_final.csv"), sep=',', header=0)
    completeness_status = pd.read_csv(os.path.join(prj_dir, "completeness_status.csv"), sep=',', header=0, index_col=None)
    
    # processing vs2 and gnm reports; adding names
    cat = cat[cat["classification"]=="taxid assigned"]
    cat.columns = ["cat_" + x for x in cat.columns.values.tolist()]
    cat = cat.rename({"cat_contig":"seq_name", "cat_superkingdom":"cat_category"}, axis=1)
    cat = cat.loc[:,["seq_name","cat_category"]].reset_index(drop=True)
    
    vs2["completeness"] = vs2["seqname"].str.split("\|\|", expand=True)[1]
    vs2["seqname"] = vs2["seqname"].str.split("\|\|", expand=True)[0]
    vs2["category"] = "Viruses"
    vs2.columns = ["vs2_" + x for x in vs2.columns.values.tolist()]
    vs2 = vs2.rename({"vs2_seqname":"seq_name"}, axis=1)
    vs2 = vs2[vs2["vs2_max_score_group"]!="RNA"]
    vs2 = vs2.loc[:, ["seq_name", "vs2_category"]].reset_index(drop=True)
    
    gnm_v["category"] = "Viruses"
    gnm_p["category"] = "Plasmids"
    gnm = pd.merge(gnm_v, gnm_p, 'outer')
    gnm.columns = ["gnm_" + x for x in gnm.columns.values.tolist()]
    gnm = gnm.rename({"gnm_seq_name":"seq_name"}, axis=1)
    gnm = gnm.loc[:,["seq_name", "gnm_category"]].reset_index(drop=True)
    
    vlm = vlm[vlm["virus_score"]>=0.8]
    vlm["category"] = "Viruses"
    vlm.columns = ["vlm_" + x for x in vlm.columns.values.tolist()]
    vlm = vlm.rename({"vlm_seq_name":"seq_name"}, axis=1)
    vlm = vlm.loc[:, ["seq_name", "vlm_category"]].reset_index(drop=True)
    
    all = pd.merge(cat, vs2, on="seq_name", how='outer')
    all = pd.merge(all, gnm, on="seq_name", how='outer')
    all = pd.merge(all, vlm, on="seq_name", how='outer')
    all.to_csv(os.path.join(prj_dir, "out", fileHeader, "putative_summary_tmp.csv"), sep='\t', index=None)
    bash_commands_1 = [
        f"grep -F -f <(sed \'1d\' {os.path.join(prj_dir, 'out', fileHeader, 'putative_summary_tmp.csv')} | cut -f1) {completeness_status[completeness_status['fileHeader']==fileHeader].iloc[0]['path']} -w | cut -d\' \' -f1,4 | sed \'s/^>//\' | sed \'s/ /\\t/\' | sed \'s/len=//\' > {os.path.join(prj_dir, 'out', fileHeader, 'putative_contigs_length_info_tmp.csv')}\n"
    ]
    with open(os.path.join(prj_dir, "find_len_tmp.sh"), 'w') as f:
        f.writelines("#!/bin/bash\n")
        f.writelines(bash_commands_1)
    os.system(f"chmod +x {os.path.join(prj_dir, 'find_len_tmp.sh')}")
    os.system(os.path.join(prj_dir, "find_len_tmp.sh"))
    os.remove(os.path.join(prj_dir, "find_len_tmp.sh"))
    
    length = pd.read_table(os.path.join(prj_dir,"out",fileHeader,"putative_contigs_length_info_tmp.csv"), header=None, names=["seq_name", "length"])
    # merge length info into summary
    all = pd.merge(all, length, on="seq_name", how='left')
    # filter putative
    putative = all.copy()
    putative_min_len = putative[putative["length"]>=min_len].reset_index(drop=True).set_index(["seq_name","length"])
    def CountVirus(row):
        count = row[row=="Viruses"].count()
        return count
    putative_min_len["v_count"] = putative_min_len.apply(CountVirus, axis=1)
    putative_min_len = putative_min_len[((putative_min_len["v_count"]>=2) | (putative_min_len["cat_category"]=="Viruses")) & (putative_min_len["gnm_category"]!="Plasmids")].reset_index().astype({"length":int}).astype(str)
    putative_min_len.to_csv(os.path.join(prj_dir,"out",fileHeader,"putative_summary.csv"), index=None)
    os.remove(os.path.join(prj_dir,"out",fileHeader,"putative_summary_tmp.csv"))
    os.remove(os.path.join(prj_dir, "out", fileHeader, "putative_contigs_length_info_tmp.csv"))

    # extract fasta
    bash_commands_2 = [
        f"seqkit grep -f <(sed '1d' {os.path.join(prj_dir, 'out', fileHeader, 'putative_summary.csv')} | cut -f1 -d',') {completeness_status[completeness_status['fileHeader']==fileHeader].iloc[0]['path']} > {os.path.join(prj_dir, 'out', fileHeader,'putative_contigs.fasta')}\n"
    ]
    with open(os.path.join(prj_dir, "extract_putative_contigs_tmp.sh"), 'w') as f:
        f.writelines("#!/bin/bash\n")
        f.writelines(bash_commands_2)
    os.system(f"chmod +x {os.path.join(prj_dir, 'extract_putative_contigs_tmp.sh')}")
    os.system(os.path.join(prj_dir, "extract_putative_contigs_tmp.sh"))
    os.remove(os.path.join(prj_dir, "extract_putative_contigs_tmp.sh"))
    

def extract_putative_contigs_multi_samples(prj_dir, fileHeader_list, min_len=3000):
    for fileHeader in fileHeader_list:
        extract_putative_contigs_single_sample(prj_dir=prj_dir, fileHeader=fileHeader, min_len=min_len)
    
def find_rRNAs_single_file(prj_dir, fileHeader, threads=64):
    
    # pre-run check
    files_to_check = [
        os.path.join(prj_dir, 'out', fileHeader, 'putative_contigs.fasta')
    ]
    for file in files_to_check:
        if not os.path.exists(file):
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
    if os.path.exists(f"{os.path.join(prj_dir, 'out', fileHeader, 'putative_contigs.fasta')}.fai"):
        os.remove(f"{os.path.join(prj_dir, 'out', fileHeader, 'putative_contigs.fasta.fai')}")

def find_rRNAs_multi_files(prj_dir, fileHeader_list, threads=64):
    for fileHeader in fileHeader_list:
        find_rRNAs_single_file(prj_dir=prj_dir, fileHeader=fileHeader, threads=threads)
    
def extract_confirmed_contigs_single_file(prj_dir, fileHeader):
    
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
    confirmed_summary.to_csv(os.path.join(prj_dir, 'out', fileHeader, 'confirmed_summary.csv'), sep=',', index=None)

    # extract confirmed fasta
    bash_commands = [
        f"seqkit grep -f <(sed '1d' {os.path.join(prj_dir, 'out', fileHeader, 'confirmed_summary.csv')} | cut -f1 -d',') {completeness_status[completeness_status['fileHeader']==fileHeader].iloc[0]['path']} > {os.path.join(prj_dir, 'out', fileHeader, 'confirmed_contigs.fasta')}\n"
    ]
    with open(os.path.join(prj_dir, "extract_confirmed_contigs_tmp.sh"), 'w') as f:
        f.writelines("#!/bin/bash\n")
        f.writelines(bash_commands)
    os.system(f"chmod +x {os.path.join(prj_dir, 'extract_confirmed_contigs_tmp.sh')}")
    os.system(os.path.join(prj_dir, "extract_confirmed_contigs_tmp.sh"))
    os.remove(os.path.join(prj_dir, "extract_confirmed_contigs_tmp.sh"))
    
def extract_confirmed_contigs_multi_files(prj_dir, fileHeader_list):
    for fileHeader in fileHeader_list:
        extract_confirmed_contigs_single_file(prj_dir=prj_dir, fileHeader=fileHeader)
        
def merge_confirmed_contigs(prj_dir, fileHeader_list):
    with open(os.path.join(prj_dir,"merged_confirmed_contigs.fasta"), 'w') as merged_confirmed_contigs:
        for fileHeader in fileHeader_list:
            if not os.path.exists(os.path.join(prj_dir, 'out', fileHeader, 'confirmed_contigs.fasta')):
                print(f"{os.path.join(prj_dir, 'out', fileHeader, 'confirmed_contigs.fasta')} not exist. Exiting.")
                return
            with open(os.path.join(prj_dir, 'out', fileHeader, 'confirmed_contigs.fasta'), 'r') as fasta:
                sequence = []
                for line in fasta:
                    if line.startswith('>'):
                        sequence.append(f">{fileHeader}_{line[1:]}")
                    else:
                        sequence.append(line)
                merged_confirmed_contigs.writelines(sequence)

if __name__=="__main__":
    pass