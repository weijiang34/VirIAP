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
    else: s["CAT"] = False
    # Verification of VirSorter2: "./out/{}/VirSorter2_results/{}-final-viral-score.tsv"
    if os.path.exists(os.path.join(prj_dir,"out",s["fileHeader"],"VirSorter2_results",f"{s['fileHeader']}-final-viral-score.tsv")):
        s["VirSorter2"] = True
    else: s["VirSorter2"] = False
    # Verification of GeNomad: "./out/{}/GeNomad_results/final.contigs_summary/final.contigs_plasmid_summary.tsv and final.contigs_virus_summary.tsv"
    # print(f"{'.'.join(s['path'].split('/')[-1].split('.')[:-1])}")
    if os.path.exists(os.path.join(prj_dir,"out",s["fileHeader"],"GeNomad_results",f"{'.'.join(s['path'].split('/')[-1].split('.')[:-1])}_summary",f"{'.'.join(s['path'].split('/')[-1].split('.')[:-1])}_virus_summary.tsv")) and os.path.exists(os.path.join(prj_dir,"out",s["fileHeader"],"GeNomad_results",f"{'.'.join(s['path'].split('/')[-1].split('.')[:-1])}_summary",f"{'.'.join(s['path'].split('/')[-1].split('.')[:-1])}_plasmid_summary.tsv")):
        s["GeNomad"] = True
    else: s["GeNomad"] = False
    # Verification of ViraLM: "./out/{}/ViraLM_results/result_final.csv"
    if os.path.exists(os.path.join(prj_dir,"out",s["fileHeader"],"ViraLM_results",f"result_{s['fileHeader']}.csv")): 
        s["ViraLM"] = True
    else: s["ViraLM"] = False
    # check putative
    if os.path.exists(os.path.join(prj_dir,"out",s["fileHeader"],"putative_contigs.fasta")): 
        s["putative"] = True
    else: s["putative"] = False
    # check decontamination
    if os.path.exists(os.path.join(prj_dir,"out",s["fileHeader"],"rRNAs.tsv")): 
        s["decontamination"] = True
    else: s["decontamination"] = False
    # check confirmed
    if os.path.exists(os.path.join(prj_dir,"out",s["fileHeader"],"decontaminated_contigs.fasta")): 
        s["confirmed"] = True
    else: s["confirmed"] = False
    # all check
    if s[columns_to_check].all():
        s["completed"] = True
    else: s["completed"] = False
    
    return s

def check_complete_multifile(prj_dir: str):
    status = pd.read_csv(os.path.join(prj_dir,"completeness_status.csv"),sep=',',header=0,index_col=None)
    status = status.apply(lambda x: check_complete_singlefile(x, prj_dir=prj_dir), axis=1)
    status = status[["path","fileHeader","completed","CAT","VirSorter2","GeNomad","ViraLM","putative","decontamination","confirmed"]]
    status.to_csv(os.path.join(prj_dir, "completeness_status.csv"), sep=',', index=None)
    
    out = status[["completed","CAT","VirSorter2","GeNomad","ViraLM","putative","decontamination","confirmed"]]
    out = out.apply(lambda x: str(x[x==True].count())+f"/{status.shape[0]}", axis=0)
    summary = pd.DataFrame(out).T
    print(summary.to_string(index=False, justify='center', ))

