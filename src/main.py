import os
import argparse
from utils import create, config, gadi_job, check_completeness, post_process
import pandas as pd

def main():
    pwd = os.getcwd()
    parser = argparse.ArgumentParser(prog="VIP", description="This is a commandline software for the pipeline.")
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
    subparser_extract.add_argument("-c", "--num_confirmed_tools", type=int, default=2, help="Minimum number of tools to confirm.")

    subparser_filter = subparsers.add_parser("decontam", help="Decontamination: filter out rRNAs from bac,euk,arc,mito.")
    subparser_confirm = subparsers.add_parser("confirm", help="Output confirmed viral contigs.")
    subparser_merge = subparsers.add_parser("merge", help="Merge all confirmed viral contigs into one fasta file.")

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
            gadi_job.generate_jobs(project_dir=project_dir, config=proj_config, batch_size=batch_size)
        elif subparser_search.submit==True:
            gadi_job.submit_jobs(project_dir)

    if args.modules=="check":
        check_completeness.check_complete_multifile(prj_dir=project_dir)
    if args.modules in ["extract","decontam","confirm","merge"]:
        fileHeader_list = pd.read_csv(os.path.join(project_dir,"completeness_status.csv"),sep=',',header=0,index_col=None)
        if args.modules=="extract":
            post_process.extract_putative_contigs_multi_samples(prj_dir=project_dir, fileHeader_list=fileHeader_list[fileHeader_list["completed"]==False].loc[:,"fileHeader"].tolist(), min_len=args.min_length)
        if args.modules=="decontam":
            post_process.find_rRNAs_multi_files(prj_dir=project_dir, fileHeader_list=fileHeader_list[fileHeader_list["completed"]==False].loc[:,"fileHeader"].tolist(), threads=64)
        if args.modules=="confirm":
            post_process.extract_confirmed_contigs_multi_files(prj_dir=project_dir, fileHeader_list=fileHeader_list[fileHeader_list["completed"]==False].loc[:,"fileHeader"].tolist())
    if args.modules=="merge":
        post_process.merge_confirmed_contigs(prj_dir=project_dir, fileHeader_list=os.listdir(os.path.join(project_dir, "out")))

if __name__ == "__main__":
    main()
    