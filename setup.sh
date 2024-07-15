#!/bin/bash

WORKING_DIR=$(pwd)
CONDA_ENVS_PATH=$(conda info | grep 'envs directories' | cut -d':' -f2 | sed 's/^ //')
CONDA_PATH=$(conda info | grep 'envs directories' | cut -d':' -f2 | sed 's/^ //' | sed 's/\/envs$//')

install_tools() {
    echo "INFO: Creating conda env: vip ..."
    conda create -n viriap -c bioconda -c conda-forge -f $WORKING_DIR/environment.yml
    mkdir -p $WORKING_DIR/dependencies
    # CAT
    echo "INFO: [1/4] Installing CAT_pack ..."
    if [ -d $WORKING_DIR/dependencies/CAT_pack ]; then
        rm -r $WORKING_DIR/dependencies/CAT_pack
    fi 
    cd $WORKING_DIR/dependencies
    git clone https://github.com/MGXlab/CAT_pack.git
    # VirSorter2
    echo "INFO: [2/4] Installing Virsorter2 ..."
    if conda info --envs | grep -q -w "vs2"; then
        echo -e "\tEnv: 'vs2' already exists."
    else
        echo -e "\tEnv: 'vs2' not existed, creating..."
        # create vs2 environment
        conda create -n vs2 -c conda-forge -c bioconda virsorter=2 --yes
        if [ $? -eq 0 ]; then
            echo -e "\tEnv: 'vs2' created."
        fi
    fi
    # GeNomad
    echo "INFO: [3/4] Installing GeNomad ..."
    if conda info --envs | grep -q -w "genomad"; then
        echo -e "\tEnv: 'genomad' already exists."
    else
        echo -e "\tEnv: 'genomad' not existed, creating..."
        # create genomad environment
        conda create -n genomad -c conda-forge -c bioconda genomad --yes
        if [ $? -eq 0 ]; then
            echo -e "\tEnv: 'genomad' created."
        fi
    fi
    # ViraLM
    echo "INFO: [4/4] Installing ViraLM ..."
    ENV_NAME="viralm"
    if conda info --envs | grep -q -w "viralm"; then
        echo -e "\tEnv: 'viralm' already exists."
    else
        echo -e "\tEnv: 'viralm' not existed, creating..."
        cd $WORKING_DIR/dependencies
        git clone https://github.com/ChengPENG-wolf/ViraLM.git
        cd $WORKING_DIR/dependencies/ViraLM
        # create viralm environment
        conda env create -f viralm.yaml -n viralm
        if [ $? -eq 0 ]; then
            echo -e "\tEnv: 'viralm' created."
        fi
    fi
    # tool path check and env variables settings
    if [ -f $WORKING_DIR/src/envs.py ];then
        rm $WORKING_DIR/src/envs.py
    fi
    if [ -d $CONDA_PATH ];then
        echo "CONDA_PATH = \"$CONDA_PATH\"" >> $WORKING_DIR/src/envs.py
    fi
    # CAT_pack  
    if [ -f $WORKING_DIR/dependencies/CAT_pack/CAT_pack/CAT_pack ]; then
        echo -e "CAT_PACK_PATH = \"${WORKING_DIR}/dependencies/CAT_pack/CAT_pack/CAT_pack\"" >> $WORKING_DIR/src/envs.py
    else
        echo "CAT_PACK_PATH = " >> $WORKING_DIR/src/envs.py
    fi
    # Virsorter2
    if [ -f $CONDA_ENVS_PATH/vs2/bin/virsorter ]; then
        echo -e "VIRSORTER2_PATH = \"${CONDA_ENVS_PATH}/vs2/bin/virsorter\"" >> $WORKING_DIR/src/envs.py
    else 
        echo "VIRSORTER2_PATH = " >> $WORKING_DIR/src/envs.py
    fi
    # genomad
    if [ -f $CONDA_ENVS_PATH/genomad/bin/genomad ]; then
        echo -e "GENOMAD_PATH = \"${CONDA_ENVS_PATH}/genomad/bin/genomad\"" >> $WORKING_DIR/src/envs.py
    else 
        echo "GENOMAD_PATH = " >> $WORKING_DIR/src/envs.py
    fi
    # viralm
    if [ -f $WORKING_DIR/dependencies/ViraLM/viralm.py ]; then
        echo -e "VIRALM_PATH = \"${WORKING_DIR}/dependencies/ViraLM/viralm.py\"" >> $WORKING_DIR/src/envs.py
    else 
        echo "VIRALM_PATH = " >> $WORKING_DIR/src/envs.py
    fi
    # checkv
    if [ -f $CONDA_ENVS_PATH/vip/bin/checkv ]; then
        echo -e "CHECKV_PATH = \"$CONDA_ENVS_PATH/vip/bin/checkv\"" >> $WORKING_DIR/src/envs.py
    else 
        echo "CHECKV_PATH = " >> $WORKING_DIR/src/envs.py
    fi
    # strobealign
    if [ -f $CONDA_ENVS_PATH/vip/bin/strobealign ]; then
        echo -e "STROBEALIGN_PATH = \"$CONDA_ENVS_PATH/vip/bin/strobealign\"" >> $WORKING_DIR/src/envs.py
    else 
        echo "STROBEALIGN_PATH = " >> $WORKING_DIR/src/envs.py
    fi
    # samtools
    if [ -f $CONDA_ENVS_PATH/vip/bin/samtools ]; then
        echo -e "SAMTOOLS_PATH = \"$CONDA_ENVS_PATH/vip/bin/samtools\"" >> $WORKING_DIR/src/envs.py
    else 
        echo "SAMTOOLS_PATH = " >> $WORKING_DIR/src/envs.py
    fi
    # featureCounts
    if [ -f $CONDA_ENVS_PATH/vip/bin/featureCounts ]; then
        echo -e "FEATURECOUNTS_PATH = \"$CONDA_ENVS_PATH/vip/bin/featureCounts\"" >> $WORKING_DIR/src/envs.py
    else 
        echo "FEATURECOUNTS_PATH = " >> $WORKING_DIR/src/envs.py
    fi
}

prepare_databases() {
    # checkv
    echo -e "\tPreparing checkv db ..."
    if [ ! -d $WORKING_DIR/dependencies/checkvdb ]; then
        mkdir $WORKING_DIR/dependencies/checkvdb
    fi
    $CONDA_ENVS_PATH/vip/bin/checkv download_database $WORKING_DIR/dependencies/checkvdb/
    # CAT_pack nr
    echo -e "\tPreparing CAT_pack_nr_db ..."
    if [ -d $WORKING_DIR/dependencies/CAT_pack_nr_db ]; then
        rm -rf $WORKING_DIR/dependencies/CAT_pack_nr_db
    fi
    cd $WORKING_DIR/dependencies
    mkdir $WORKING_DIR/dependencies/CAT_pack_nr_db
    cd $WORKING_DIR/dependencies/CAT_pack_nr_db
    wget -c tbb.bio.uu.nl/tina/CAT_pack_prepare/20240422_CAT_nr.tar.gz
    tar -xvzf 20240422_CAT_nr.tar.gz
    # vs2 db
    echo -e "\tPreparing vs2_db ..."
    if [ -d $WORKING_DIR/dependencies/vs2_db ]; then
        rm -rf $WORKING_DIR/dependencies/vs2_db
    fi
    if conda info --envs | grep -q -w "vs2"; then
        cd $WORKING_DIR/dependencies
        source $CONDA_PATH/bin/activate vs2 
        virsorter setup -d vs2_db -j 4
        conda deactivate
    else
        echo -e "[WARNING]: Conda environment: 'vs2' not existed, vs2_db installation failed."
    fi
    # genomad db
    echo -e "\tPreparing genomad db ..."
    if [ -d $WORKING_DIR/dependencies/genomad_db ]; then
        rm -rf $WORKING_DIR/dependencies/genomad_db
    fi
    if conda info --envs | grep -q -w "genomad"; then
        source $CONDA_PATH/bin/activate genomad
        genomad download-database ./dependencies/
        conda deactivate
    else
        echo -e "[WARNING]: Conda environment: 'genomad' not existed, genomad db installation failed."
    fi
    # viralm model
    echo -e "\tPreparing viralm model ..."
    if [ -d $WORKING_DIR/dependencies/ViraLM/model ]; then
        rm -rf $WORKING_DIR/dependencies/ViraLM/model
    fi 
    if conda info --envs | grep -q -w "viralm"; then
        cd $WORKING_DIR/dependencies/ViraLM
        source $CONDA_PATH/bin/activate viralm
        # download and setup the model
        gdown --id 1EQVPmFbpLGrBLU0xCtZBpwvXrtrRxic1
        tar -xzvf model.tar.gz -C .
        rm model.tar.gz
        conda deactivate
    else
        echo -e "[WARNING]: Conda environment: 'viralm' not existed, viralm model installation failed."
    fi

    # db path check and env variables settings
    # CAT_pack_nr_db
    if [ -d $WORKING_DIR/dependencies/CAT_pack_nr_db ]; then
        echo -e "CAT_PACK_DB_PATH = \"${WORKING_DIR}/dependencies/CAT_pack_nr_db\"" >> $WORKING_DIR/src/envs.py
    else
        echo "CAT_PACK_DB_PATH = " >> $WORKING_DIR/src/envs.py
    fi
    # genomad_db
    if [ -d $WORKING_DIR/dependencies/genomad_db ]; then
        echo -e "GENOMAD_DB_PATH = \"${WORKING_DIR}/dependencies/genomad_db\"" >> $WORKING_DIR/src/envs.py
    else 
        echo "GENOMAD_DB_PATH = " >> $WORKING_DIR/src/envs.py
    fi
    # checkvdb
    if [ -d $WORKING_DIR/dependencies/checkvdb ] && compgen -d $WORKING_DIR/dependencies/checkvdb > /dev/null; then
        echo -e "CEHCKV_DB_PATH = \"${$compgen -d $WORKING_DIR/dependencies/checkvdb}\"" >> $WORKING_DIR/src/envs.py
    else 
        echo "CEHCKV_DB_PATH = " >> $WORKING_DIR/src/envs.py
    fi

}

case $1 in
    "--all")
        install_tools
        prepare_databases
        ;;
    "--tools")
        install_tools
        ;;
    "--databases")
        prepare_databases
        ;;
    *)
        echo -e "options: '--tools' '--databases' '--all'"
        ;;
esac