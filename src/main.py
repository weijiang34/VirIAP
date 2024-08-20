import os
import argparse
from utils import create, config, job_management, gadi_job, check_completeness, post_process, mapping, classify
import pandas as pd

def main():
    pwd = os.getcwd()
    parser = argparse.ArgumentParser(prog="VirIAP", description="This is a commandline software for the pipeline.")
    parser.add_argument("-p","--prj_dir", default=pwd, help="Specify your project directory. (e.g. './test')")
    parser.add_argument("--max_batch_size", type=int, default=10, help="Defines the maximum number of samples a batch job contains.")
    parser.add_argument("--dryrun", action="store_true", help="Run the commands but not generate anything. (Used for testing.)")
    # subparser_search.add_argument("-c", "--config", help="The config file of this project in yaml format.")

    subparsers  = parser.add_subparsers(
        title="Modules",
        dest="modules",
        description="Modules that proceed different functions.", 
        help="Please specify one and only one of these options.", 
        required=True
    )

    subparser_create = subparsers.add_parser("create", help="Create a project under [-p/--porject_dir]")
    subparser_create.add_argument("-i","--input", nargs="+", action="extend", help="A file contains a list of absolute paths to fasta files; OR one or more fasta files")
    
    subparser_search = subparsers.add_parser("search", help="A submodule which generate/submit executable job files;")
    search_option = subparser_search.add_mutually_exclusive_group(required=True)
    search_option.add_argument("--generate", action="store_true", help="Generate jobs to ./jobs")
    search_option.add_argument("--submit", action="store_true", help="Submit jobs under ./jobs")
    
    subparser_check = subparsers.add_parser("check", help ="Check completeness status of all samples.")

    subparser_extract = subparsers.add_parser("extract", help="Extract putative contigs.")
    subparser_extract.add_argument("-l", "--min_length", type=int, default=3000, help="Minimum length for putative contigs.")
    subparser_extract.add_argument("-c", "--num_tools", type=int, default=2, help="Minimum number of tools to confirm.")

    subparser_filter = subparsers.add_parser("decontam", help="Decontamination: filter out rRNAs from bac,euk,arc,mito.")
    # subparser_confirm = subparsers.add_parser("confirm", help="Output confirmed viral contigs.")
    subparser_merge = subparsers.add_parser("merge", help="Merge all confirmed viral contigs into one fasta file.")
    subparser_dedup = subparsers.add_parser("dedup", help="Remove exactly the same contigs.")
    subparser_quality_check = subparsers.add_parser("check_quality", help="Use CheckV to check the quality of merged viral contigs.")
    subparser_cluster = subparsers.add_parser("cluster", help="Use ANI and AF results from blast all against all to cluster viral contigs.")
    subparser_cluster.add_argument("--no_checkv", action="store_true", help="If specified, will use un-checkv-checked contigs for clustering instead of using checkv-checked contigs.")
    
    subparser_mapping = subparsers.add_parser("mapping", help="Map clean paired-end reads to representative contigs using strobealign, and calculate relative abundance.")
    subparser_mapping.add_argument("--manifest", required=True, type=str, help="a three column csv file, with columns: fileHeader,fq1,fq2")
    mapping_option = subparser_mapping.add_mutually_exclusive_group(required=True)
    mapping_option.add_argument("--indexing", action="store_true", help="Generate jobs for building strobealign index with representative viral contigs.")
    mapping_option.add_argument("--mapping", action="store_true", help="Generate jobs for mapping batches of samples to the representatives.")
    mapping_option.add_argument("--count_matrix", action="store_true", help="Count number of reads and calculate FPKM and TPM.")

    subparser_classify = subparsers.add_parser("classify", help="Classify OVU representatives and merge lineage results from CAT, GeNomad and vContact3.")
    calssification_option = subparser_classify.add_mutually_exclusive_group(required=True)
    calssification_option.add_argument("--generate_job", action="store_true", help="generate classiification job to run vContact3.")
    calssification_option.add_argument("--merge_lineage", action="store_true", help="Merge lineage info from CAT, GeNomad, and vContact3.")

    args = parser.parse_args()

    project_dir = os.path.abspath(args.prj_dir)

    # Whether to dryrun
    if args.dryrun==True:
        print("Dryrun for testing.")
        return
    # Begin execution
    if args.modules=="create":
        create.create_project(project_dir)
        create.parse_input(project_dir, args.input)
        create.make_sample_dirs(project_dir)
        check_completeness.check_complete_multifile(prj_dir=project_dir)

    if args.modules=="search":
        proj_config = config.read_project_config(os.path.join(project_dir,"config.yaml"))
        if args.generate==True:
            if args.max_batch_size==None: batch_size=proj_config["max_batch_size"]
            else: batch_size=args.max_batch_size
            job_management.generate_jobs(project_dir=project_dir, config=proj_config, batch_size=batch_size)
        elif args.submit==True:
            print("Please submit the jobs manually, just to make sure all the resoures requiered are valid and proper.")
            # gadi_job.submit_jobs(project_dir)

    if args.modules=="check":
        check_completeness.check_complete_multifile(prj_dir=project_dir)
    if args.modules in ["extract","decontam","confirm"]:
        fileHeader_list = pd.read_csv(os.path.join(project_dir,"completeness_status.csv"),sep=',',header=0,index_col=None)
        if args.modules=="extract":
            post_process.extract_putative_contigs_multi_samples(prj_dir=project_dir, fileHeader_list=fileHeader_list.loc[:,"fileHeader"].tolist(), min_len=args.min_length, num_tools=args.num_tools)
        if args.modules=="decontam":
            proj_config = config.read_project_config(os.path.join(project_dir,"config.yaml"))
            if proj_config["job_manager"] in ["pbs","gadi"]:
                post_process.find_rRNAs_multi_files(prj_dir=project_dir, fileHeader_list=fileHeader_list.loc[:,"fileHeader"].tolist(), threads=proj_config['pbs']["ncpus"]) # [fileHeader_list["completed"]==False]
                post_process.extract_decontaminated_contigs_multi_files(prj_dir=project_dir, fileHeader_list=fileHeader_list.loc[:,"fileHeader"].tolist())
    if args.modules=="merge":
        fileHeader_list = pd.read_csv(os.path.join(project_dir,"completeness_status.csv"),sep=',',header=0,index_col=None)
        post_process.merge_confirmed_contigs(prj_dir=project_dir, fileHeader_list=fileHeader_list.loc[:,"fileHeader"].tolist())
    if args.modules=="dedup":
        post_process.dedup(prj_dir=project_dir)
    if args.modules=="check_quality":
        proj_config = config.read_project_config(os.path.join(project_dir,"config.yaml"))
        post_process.check_quality(prj_dir=project_dir, config=proj_config)
    if args.modules=="cluster":
        proj_config = config.read_project_config(os.path.join(project_dir,"config.yaml"))
        if args.no_checkv:
            post_process.cluster(prj_dir=project_dir, config=proj_config, checked="unchecked")
        else:
            post_process.cluster(prj_dir=project_dir, config=proj_config, checked="checked")
    if args.modules=="mapping":
        proj_config = config.read_project_config(os.path.join(project_dir,"config.yaml"))
        if args.indexing==True:
            mapping.indexing(prj_dir=project_dir, config=proj_config)
        if args.mapping==True:
            mapping.mapping(prj_dir=project_dir, manifest=args.manifest, config=proj_config)
        if args.count_matrix==True:
            mapping.count_matrix(prj_dir=project_dir, manifest=args.manifest)
    if args.modules=="classify":
        os.makedirs(os.path.join(project_dir,"Classification"), exist_ok=True)
        proj_config = config.read_project_config(os.path.join(project_dir,"config.yaml"))
        if args.generate_job==True:
            classify.anno_vContact3(prj_dir=project_dir, config=proj_config)
        if args.merge_lineage==True:
            classify.summarise_OVUs(prj_dir=project_dir)

    

if __name__ == "__main__":
    main()
    