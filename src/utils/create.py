import os
from utils import config
import pandas as pd

def create_project(proj_dir):
    # initialize
    if not os.path.exists(proj_dir):
        print(f"Creating project: {proj_dir}")
        os.mkdir(proj_dir)
    else:
        print(f"{proj_dir} already exists. Switching to project: {os.path.basename(proj_dir)}")

    project_completeness = {
        "out": False,
        "jobs": False,
        "logs": False,
        "config.yaml": False
    }
    if os.path.exists(os.path.join(proj_dir, "out")):
        project_completeness["out"] = True
    else:
        os.makedirs(os.path.join(proj_dir, "out"))
        
    if os.path.exists(os.path.join(proj_dir, "jobs")):
        project_completeness["jobs"] = True
    else:
        os.makedirs(os.path.join(proj_dir, "jobs"))
        
    if os.path.exists(os.path.join(proj_dir, "logs")):
        project_completeness["logs"] = True
    else:
        os.makedirs(os.path.join(proj_dir, "logs"))
        
    if os.path.exists(os.path.join(proj_dir, "config.yaml")):
        project_completeness["config.yaml"] = True
    else:
        config.init_project_config(os.path.join(proj_dir, "config.yaml"))

def parse_input(proj_dir, input):
    if input==None:
        file_list = []
    elif len(input)==1:
        with open(input[0],'r') as f:
            first_line = f.readline()
            f.seek(0)
            if first_line.startswith(">"):
                file_list = input
            else:
                file_list = [s.rstrip('\n') for s in f.readlines()]
    else:
        file_list = input
    print(f"Total {len(file_list)} files input.")
    df = pd.DataFrame({
        "path": file_list,
        # to be refined
        "fileHeader": [file_path.split("/")[-1].split(".")[0] for file_path in file_list],
        "completed": False,
        "CAT": False,
        "VirSorter2": False,
        "GeNomad": False,
        "ViraLM": False,
        "putative": False,
        "decontamination": False,
        "confirmed": False,
    })
    df.to_csv(os.path.join(proj_dir,"completeness_status.csv"), sep=',', index=None)
    
def make_sample_dirs(proj_dir):
    df = pd.read_csv(os.path.join(proj_dir,"completeness_status.csv"), header=0, index_col=0, sep=',')
    fileHeader_list = df["fileHeader"].to_list()
    os.makedirs(os.path.join(proj_dir,"OVU"), exist_ok=True)
    os.makedirs(os.path.join(proj_dir,"Abundance"), exist_ok=True)
    os.makedirs(os.path.join(proj_dir,"Classification"), exist_ok=True)
    for fileHeader in fileHeader_list:
        os.makedirs(os.path.join(proj_dir,"out",fileHeader,"CAT_results"), exist_ok=True)
        os.makedirs(os.path.join(proj_dir,"out",fileHeader,"VirSorter2_results"), exist_ok=True)
        os.makedirs(os.path.join(proj_dir,"out",fileHeader,"GeNomad_results"), exist_ok=True)
        os.makedirs(os.path.join(proj_dir,"out",fileHeader,"ViraLM_results"), exist_ok=True)
         
        
if __name__=="__main__":
    pass