import os
import envs
import numpy as np
import pandas as pd
from Bio import SeqIO
from utils import job_management

def anno_vContact3(prj_dir, config):
    os.makedirs(os.path.join(prj_dir,"Classification"), exist_ok=True)
    job_dir = f"{prj_dir}/Classification"
    log_dir = f"{prj_dir}/Classification"
    threads = config['pbs']['ncpus']

    bash_commands = [
        f"source {envs.CONDA_PATH}/bin/activate vcontact3",
        f"{envs.VCONTACT3_PATH} run --nucleotide {prj_dir}/OVU/rep_contigs.fasta --output {prj_dir}/Classification/anno_vcontact3 --db-domain \"prokaryotes\" --db-version {envs.VCONTACT3_DB_VERSION} --db-path {envs.VCONTACT3_DB_PATH} -t {threads}", # -e {{cytoscape,d3js,tree,pyupset,profiles}}
    ]
    bash_commands = [x+"\n" for x in bash_commands]
    if config['job_manager']=='pbs':
        cluster_job_header = job_management.PBSHeader(
            job_name="classify_vContact3",
            ncpus=threads,
            ngpus=0,
            mem="32GB",
            walltime="12:00:00",
            mail_addr=config['pbs']['mail_addr'],
            log_o=os.path.join(log_dir, f"{log_dir}/classify_vContact3.o"),
            log_e=os.path.join(log_dir, f"{log_dir}/classify_vContact3.e"),
        )
        cluster_job = job_management.Job(
            job_manager='pbs',job_header=cluster_job_header, commands=bash_commands,
        )
        cluster_job.save_job(job_dir=job_dir)
    elif config['job_manager']=='gadi':
        cluster_job_header = job_management.GadiHeader(
            job_name=f"classify_vContact3",
            ncpus=threads,
            ngpus=0,
            mem="32GB",
            walltime="12:00:00",
            mail_addr=config['pbs']['mail_addr'],
            log_o=os.path.join(log_dir, f"{log_dir}/classify_vContact3.o"),
            log_e=os.path.join(log_dir, f"{log_dir}/classify_vContact3.e"),
            project=config['pbs']['gadi']['-P project'],
            storage=config['pbs']['gadi']['-l storage'],
            node_type="normalsl",
            jobfs="2GB",
        )
        cluster_job = job_management.Job(
            job_manager='gadi',job_header=cluster_job_header, commands=bash_commands,
        )
        cluster_job.save_job(job_dir=job_dir)

    return

def summarise_OVUs(prj_dir, include=["CAT", "VCT", "GNM"]):
    # include CAT and genomad annotations
    def include_CAT_genomad(filtered_clusters_path, OVU_info_tmp_path):
        ovu = pd.read_table(filtered_clusters_path, sep='\t', header=None, names=["representative_contig", "contigs_in_cluster"])
        ovu["OVU"] = [f"OVU_{x}" for x in range(ovu.shape[0])]
        ovu["fileHeader"] = ovu["representative_contig"].str.split("_").apply(lambda x: '_'.join(x[:-2]))
        ovu["contig"] = ovu["representative_contig"].str.split("_").apply(lambda x: '_'.join([x[-2], x[-1]]))
        ovu["cluster_size"] = ovu["contigs_in_cluster"].str.split(",").apply(len)
        ovu = ovu.sort_values(by="representative_contig", ascending=True)
        status = pd.read_csv(os.path.join(prj_dir,"completeness_status.csv"), header=0)
        ovu_fileHeaders = []
        for idx, fileHeader in enumerate(ovu["fileHeader"].unique().tolist()):
            print(idx, fileHeader)
            cat_out = pd.read_table(
                os.path.join(prj_dir,"out",fileHeader,"CAT_results",f"{fileHeader}.nr.contig2classification.with_names.txt"),
                # f"../out/{fileHeader}/CAT_results/{fileHeader}.nr.contig2classification.with_names.txt",
                header=0, sep='\t', index_col=None
            ).rename({"# contig":"contig"}, axis=1)
            gnm_fileHeader = '.'.join(status[status["fileHeader"]==fileHeader]['path'].iloc[0].split('/')[-1].split('.')[:-1])
            gnm_out = pd.read_table(
                os.path.join(prj_dir,"out",fileHeader,"GeNomad_results",f"{gnm_fileHeader}_annotate",f"{gnm_fileHeader}_taxonomy.tsv"),
                # f"../out/{fileHeader}/GeNomad_results/final.contigs_annotate/final.contigs_taxonomy.tsv",
                header=0, sep='\t', index_col=None
            ).rename({"seq_name":"contig", "lineage":"genomad_lineage"}, axis=1)
            ovu_current_fileHeader = ovu[ovu["fileHeader"]==fileHeader]
            ovu_current_fileHeader_annotation = pd.merge(ovu_current_fileHeader, cat_out, on="contig", how='left')
            ovu_current_fileHeader_annotation = pd.merge(ovu_current_fileHeader_annotation, gnm_out, on="contig", how='left')
            ovu_fileHeaders.append(ovu_current_fileHeader_annotation)
        ovu_annotations = pd.concat(ovu_fileHeaders, axis=0).reset_index(drop=True)
        ovu_annotations = ovu_annotations.set_index(["OVU", "representative_contig", "cluster_size"]).reset_index()
        ovu_annotations = ovu_annotations.sort_values(by="cluster_size", ascending=False).reset_index(drop=True)
        ovu_annotations = ovu_annotations.astype({"superkingdom":str,"phylum":str,"class":str,"order":str,"family":str,"genus":str,"species":str})
        ovu_annotations["cat_lineage"] = ovu_annotations.apply(lambda row: ';'.join([row["superkingdom"].split(':')[0],row["phylum"].split(':')[0],row["class"].split(':')[0],row["order"].split(':')[0],row["family"].split(':')[0],row["genus"].split(':')[0],row["species"].split(':')[0]]), axis=1)
        ovu_annotations = ovu_annotations.drop(["classification","superkingdom","phylum","class","order","family","genus","species"], axis=1)
        tmp_df = ovu_annotations["contigs_in_cluster"]
        ovu_annotations = ovu_annotations.drop("contigs_in_cluster", axis=1)
        ovu_annotations["contigs_in_cluster"] = tmp_df
        ovu_annotations.to_csv(OVU_info_tmp_path, index=None)
        return ovu_annotations
    # include vcontact3 annotations
    def incldue_vContact3(OVU_info_tmp_path, reps_lineage_path):
        ovu_annotations = pd.read_csv(OVU_info_tmp_path, header=0)
        # vcontact3_summary = pd.read_csv(os.path.join(prj_dir,"Classification","anno_vcontact3","final_assignments.csv"), header=0).fillna("").astype({"family (prediction)":"str"})
        # vcontact3 = vcontact3_summary[vcontact3_summary["Reference"]==False].drop("index",axis=1).reset_index(drop=True).rename({
        #     "GenomeName":"representative_contig",
        #     "realm (prediction)": "realm",
        #     "phylum (prediction)": "phylum",
        #     "class (prediction)": "class",
        #     "order (prediction)": "order",
        #     "family (prediction)": "family",
        #     "genus (prediction)": "genus",
        # },axis=1)
        vcontact3 = pd.DataFrame({
            "index":[],
            "representative_contig":[],
            "RefSeqID":[],
            "Proteins":[],
            "Reference":[],
            "Size (Kb)":[],
            "realm (Reference)":[],
            "realm":[],
            "phylum (Reference)":[],
            "phylum":[],
            "class (Reference)":[],
            "class":[],
            "order (Reference)":[],
            "order":[],
            "family (Reference)":[],
            "family":[],
            "subfamily (Reference)":[],
            "subfamily":[],
            "genus (Reference)":[],
            "genus":[],
            "network":[],
        })
        
        # vcontact3["vContact3_lineage"] = vcontact3.apply(lambda row: ";".join([row["realm"], row["phylum"], row["class"], row["order"], row["family"], row["genus"]]), axis=1)
        vcontact3["vContact3_lineage"] = vcontact3.apply(lambda row: "", axis=1)

        reps_lineage = pd.DataFrame({
            "OVU": ovu_annotations["OVU"],
            "representative_contig": ovu_annotations["representative_contig"],
            "cluster_size": ovu_annotations["cluster_size"],
            "CAT_lineage": ovu_annotations["cat_lineage"],#.str.replace("no support",""),
            # "vContact3_lineage": vcontact3.apply(lambda x: ";".join([str(x["realm"]),str(x["phylum"]),str(x["class"]),str(x["order"]),str(x["family"]),str(x["genus"])]), axis=1),
            "GeNomad_lineage": ovu_annotations["genomad_lineage"],
            "contigs_in_cluster": ovu_annotations["contigs_in_cluster"]
        })
        reps_lineage = pd.merge(reps_lineage, vcontact3.loc[:, ["representative_contig", "vContact3_lineage"]], how='left', on='representative_contig')
        reps_lineage = reps_lineage.loc[:,["OVU", "representative_contig","cluster_size","CAT_lineage","vContact3_lineage","GeNomad_lineage","contigs_in_cluster"]]
        
        reps_lineage["CAT_lineage"] = reps_lineage["CAT_lineage"].apply(lambda x: ";".join(["" if value in [None, "nan", "no support"] else value for value in str(x).split(";")]))
        # print(reps_lineage["CAT_lineage"])
        reps_lineage["CAT_lineage"] = reps_lineage["CAT_lineage"].apply(lambda x: ";".join([value.replace("*","") for value in str(x).split(";")]))
        # print(reps_lineage["CAT_lineage"])
        reps_lineage["CAT_lineage"] = reps_lineage["CAT_lineage"].apply(lambda x: ";".join([""]*7 if "Viruses" not in str(x).split(";") else str(x).split(";")))
        # print(reps_lineage["CAT_lineage"])
        # print(reps_lineage["vContact3_lineage"])
        reps_lineage["vContact3_lineage"] = reps_lineage["vContact3_lineage"].apply(lambda x: ";".join(["" if (value in [None, "nan", "no support", "No Realm", "No prediction"]) or ("novel" in value) else value for value in str(x).split(";")]))
        # print(reps_lineage["vContact3_lineage"])
        reps_lineage["vContact3_lineage"] = reps_lineage["vContact3_lineage"].apply(lambda x: ";".join([value.split("|")[0] for value in str(x).split(";")]))
        reps_lineage["GeNomad_lineage"] = reps_lineage["GeNomad_lineage"].apply(lambda x: ";".join(["" if value in [None, "nan", "no support"] else value for value in str(x).split(";")]))
        # print(reps_lineage["vContact3_lineage"])
        reps_lineage.to_csv(reps_lineage_path, index=None)
        return reps_lineage

    def merge_lineage(lineage_list):
        ranks = ["phylum","class","order","family","genus","species"]
        current_lineage = {"phylum":"","class":"","order":"","family":"","genus":"","species":""}

        def find_lowest_valid_rank(lineage):
            for idx, key in enumerate(reversed(list(lineage.keys()))):
                if lineage[key] in ["", None, "nan", "no support"]:
                    continue
                else:
                    return key
            return key

        def fill_lineage(lineage, template):
            for rank in ranks[ranks.index(find_lowest_valid_rank(lineage)):]:
                if rank in template.keys():
                    lineage[rank] = template[rank]
                else: lineage[rank] = ""
            return lineage

        for lineage in lineage_list:
            if find_lowest_valid_rank(current_lineage)=="phylum" or lineage[find_lowest_valid_rank(current_lineage)]==current_lineage[find_lowest_valid_rank(current_lineage)]:
                current_lineage = fill_lineage(current_lineage, lineage)
            else:
                continue
        return current_lineage

    def to_lineage(str, tool):
        lineage = {"phylum":"","class":"","order":"","family":"","genus":"","species":""}
        if tool=="CAT":
            CAT_lineage = str.split(";")
            if len(CAT_lineage)>=2:
                for idx,value in enumerate(CAT_lineage[1:]):
                    lineage[list(lineage.keys())[idx]] = value
        if tool=="vContact3":
            vcontact3_lineage = str.split(";")
            if len(vcontact3_lineage)>=2:
                for idx,value in enumerate(vcontact3_lineage[1:]):
                    lineage[list(lineage.keys())[idx]] = value
        if tool=="GeNomad":
            genoamd_lineage = str.split(";")
            if len(genoamd_lineage)>=4:
                for idx,value in enumerate(genoamd_lineage[3:]):
                    lineage[list(lineage.keys())[idx]] = value
        return lineage

    def length_info(OVU_info, quality_filtered_contigs_path):
        contig_id_list = []
        contig_length_list = []
        with open(quality_filtered_contigs_path, 'r') as f:
            for record in SeqIO.parse(f, "fasta"):
                contig_id_list.append(record.id)
                contig_length_list.append(len(record.seq))
        df_contig_length = pd.DataFrame({
            "contig_id": contig_id_list,
            "contig_length": contig_length_list,
        })
        OVU_info["cluster_length"] = OVU_info.apply(lambda s: df_contig_length[df_contig_length["contig_id"].isin(s["contigs_in_cluster"].split(","))]["contig_length"].sum(), axis=1)
        OVU_info["cluster_median_length"] = OVU_info.apply(lambda s: np.median(df_contig_length[df_contig_length["contig_id"].isin(s["contigs_in_cluster"].split(","))]["contig_length"].tolist()), axis=1)
        OVU_info["cluster_mean_length"] = OVU_info["cluster_length"]/OVU_info["cluster_size"]
        OVU_info = OVU_info.loc[:,['OVU', 'representative_contig', 'cluster_size', 'cluster_length', 'cluster_median_length', 'cluster_mean_length', 'lineage', 'contigs_in_cluster']]
        
        return OVU_info

    filtered_clusters_path = os.path.join(prj_dir,"OVU","filtered_clusters.tsv")
    OVU_info_tmp_path = os.path.join(prj_dir,"OVU","OVU_info_tmp.csv")
    reps_lineage_path = os.path.join(prj_dir,"OVU","reps_lineage.csv")
    OVU_info_path = os.path.join(prj_dir,"OVU","OVU_info.csv")

    # include_CAT_genomad(filtered_clusters_path, OVU_info_tmp_path=OVU_info_tmp_path)
    incldue_vContact3(OVU_info_tmp_path=OVU_info_tmp_path, reps_lineage_path=reps_lineage_path)
    
    reps_lineage = pd.read_csv(reps_lineage_path, header=0)
    tool_name_map = {"CAT":"CAT", "VCT":"vContact3", "GNM":"GeNomad"}
    # lineage_list = [to_lineage(str(x[tool_name_map[tool]+"_lineage"])) for tool in include]
    # reps_lineage["lineage"] = reps_lineage.apply(lambda x: ";".join(merge_lineage([to_lineage(str(x["CAT_lineage"]),"CAT"),to_lineage(str(x["vContact3_lineage"]),"vContact3"),to_lineage(str(x["GeNomad_lineage"]),"GeNomad")]).values()), axis=1)
    reps_lineage["lineage"] = reps_lineage.apply(lambda x: ";".join(merge_lineage([to_lineage(str(x[tool_name_map[tool]+"_lineage"]), tool_name_map[tool]) for tool in include]).values()), axis=1)
    reps_lineage.to_csv(reps_lineage_path, index=None)

    OVU_info = reps_lineage.drop(["CAT_lineage","vContact3_lineage","GeNomad_lineage"],axis=1).loc[:,["OVU","representative_contig","cluster_size","lineage","contigs_in_cluster"]]
    
    OVU_info = length_info(
        OVU_info=OVU_info, 
        quality_filtered_contigs_path=os.path.join(prj_dir,"OVU","quality_filtered_viral_contigs.fasta")
    )
    
    OVU_info.to_csv(OVU_info_path,index=None)

    # os.remove(OVU_info_tmp_path)
    
    return

def save_rank_level_relative_abundance_TPM(OVU_info_path, all_TPM_path, relative_abundance_path):
        # save rank level TPMs for OVUs to excel
        OVU_info = pd.read_csv(OVU_info_path, header=0)
        all_TPM = pd.read_csv(all_TPM_path, header=0).drop(["Chr","Start","End","Strand","Length"], axis=1).rename({"contig": "representative_contig"},axis=1)
        relative_abundance = pd.DataFrame({
            "OVU": OVU_info["OVU"],
            "representative_contig": OVU_info["representative_contig"],
            "cluster_size": OVU_info["cluster_size"],
            "phylum": OVU_info["lineage"].str.split(";",expand=True)[0].replace("","unclassified/unassigned"),
            "class": OVU_info["lineage"].str.split(";",expand=True)[1].replace("","unclassified/unassigned"),
            "order": OVU_info["lineage"].str.split(";",expand=True)[2].replace("","unclassified/unassigned"),
            "family": OVU_info["lineage"].str.split(";",expand=True)[3].replace("","unclassified/unassigned"),
            "genus": OVU_info["lineage"].str.split(";",expand=True)[4].replace("","unclassified/unassigned"),
            "species": OVU_info["lineage"].str.split(";",expand=True)[5].replace("","unclassified/unassigned"),
        })
        relative_abundance = pd.merge(left=relative_abundance, right=all_TPM, how='left', on='representative_contig').drop(["OVU","representative_contig","cluster_size"], axis=1)
        with pd.ExcelWriter(relative_abundance_path, mode='w') as writer:
            ranks = relative_abundance.columns[:6]
            for rank in ranks.tolist():
                rank_level_relative_abundance = relative_abundance.drop(ranks.drop(rank), axis=1).groupby(rank).agg('sum')
                rank_level_relative_abundance.to_excel(writer, sheet_name=rank, engine='openpyxl')

        return

if __name__=='__main__':
    pass
    # anno_vContact3(prj_dir=PRJ_DIR, threads=THREADS)