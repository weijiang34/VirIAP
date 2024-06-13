import pandas as pd
import numpy as np
import os

def check_complete_singlefile(series, prj_dir: str="."):
    # if all final files exist, then complete
    s = series
    columns_to_check = ["CAT", "VirSorter2", "GeNomad", "ViraLM", "putative", "decontamination", "confirmed"]
    # Verification of CAT: "./out/{}/CAT_results/{}.nr.contig2classification.with_names.txt"
    if os.path.exists(os.path.join(prj_dir,"out",s["fileHeader"],"CAT_results",f"{s['fileHeader']}.nr.contig2classification.with_names.txt")):
        s["CAT"] = True
    # Verification of VirSorter2: "./out/{}/VirSorter2_results/{}-final-viral-score.tsv"
    if os.path.exists(os.path.join(prj_dir,"out",s["fileHeader"],"VirSorter2_results",f"{s['fileHeader']}-final-viral-score.tsv")):
        s["VirSorter2"] = True
    # Verification of GeNomad: "./out/{}/GeNomad_results/final.contigs_summary/final.contigs_plasmid_summary.tsv and final.contigs_virus_summary.tsv"
    if os.path.exists(os.path.join(prj_dir,"out",s["fileHeader"],"GeNomad_results","final.contigs_summary","final.contigs_virus_summary.tsv")) and os.path.exists(os.path.join(prj_dir,"out",s["fileHeader"],"GeNomad_results","final.contigs_summary","final.contigs_plasmid_summary.tsv")):
        s["GeNomad"] = True
    # Verification of ViraLM: "./out/{}/ViraLM_results/result_final.csv"
    if os.path.exists(os.path.join(prj_dir,"out",s["fileHeader"],"ViraLM_results","result_final.csv")): 
        s["ViraLM"] = True
    # check putative
    if os.path.exists(os.path.join(prj_dir,"out",s["fileHeader"],"putative_contigs.fasta")): 
        s["putative"] = True
    # check decontamination
    if os.path.exists(os.path.join(prj_dir,"out",s["fileHeader"],"rRNAs.tsv")): 
        s["decontamination"] = True
    # check confirmed
    if os.path.exists(os.path.join(prj_dir,"out",s["fileHeader"],"confirmed_contigs.fasta")): 
        s["confirmed"] = True
    
    if s[columns_to_check].all():
        s["completed"] = True
    
    return s

def check_complete_multifile(prj_dir: str):
    status = pd.read_csv(os.path.join(prj_dir,"completeness_status.csv"),sep=',',header=0,index_col=None)
    status = status.apply(lambda x: check_complete_singlefile(x, prj_dir=prj_dir), axis=1)
    status = status[["path","fileHeader","completed","CAT","VirSorter2","GeNomad","ViraLM","putative","decontamination","confirmed"]]
    status.to_csv(os.path.join(prj_dir, "completeness_status.csv"), sep=',', index=None)

