import os
from ruamel.yaml import YAML

def init_project_config(path):
    project_dir = os.path.abspath(os.path.split(path)[0])
    project_name = os.path.basename(project_dir)
    default_config = {
        "project_dir": project_dir,
        "name": project_name,
        "max_batch_size": 10,
        "job_manager": "pbs",
        "pbs": {
            "ncpus": 32,
            "mail_addr": "",
            "gadi": {
                "-l storage": "",
                "-P project": "",
            }
        },
        "bash": {
            "ncpus": 32
        }
    }
    yaml = YAML()
    with open(path, 'w', encoding="utf-8") as f:
        yaml.dump(default_config, f)

def read_project_config(path):
    yaml = YAML()
    with open(path, 'r', encoding="utf-8") as f:
        config = yaml.load(f)
    return config

if __name__=="__main__":
    pass