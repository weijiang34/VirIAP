#!/bin/bash

WORKING_DIR=$(pwd)
CONDA_ENVS_PATH=$(conda info | grep 'envs directories' | cut -d':' -f2 | sed 's/^ //')
CONDA_PATH=$(conda info | grep 'envs directories' | cut -d':' -f2 | sed 's/^ //' | sed 's/\/envs$//')
MAIN_ENV_NAME="viriap"

install_tools() {
    
    if conda info --envs | grep -q -w "$MAIN_ENV_NAME"; then
        echo -e "Env: '$MAIN_ENV_NAME' already exists."
    else
        echo "INFO: Creating conda env: $MAIN_ENV_NAME ..."
        conda create -n $MAIN_ENV_NAME pandas ruamel.yaml strobealign samtools --yes
        source $CONDA_PATH/bin/activate $MAIN_ENV_NAME
        conda activate $MAIN_ENV_NAME
        conda install bioconda::barrnap bioconda::seqkit bioconda::subread bioconda::checkv --yes
        if [ $? -eq 0 ]; then
            echo -e "\tEnv: '$MAIN_ENV_NAME' created."
        fi
    fi
    # conda env create -n viriap -f environment.yml
    mkdir -p $WORKING_DIR/dependencies
    # CAT
    echo "INFO: [1/5] Installing CAT_pack ..."
    if [ -d $WORKING_DIR/dependencies/CAT_pack ]; then
        rm -rf $WORKING_DIR/dependencies/CAT_pack
    fi 
    cd $WORKING_DIR/dependencies
    git clone https://github.com/MGXlab/CAT_pack.git
    # VirSorter2
    echo "INFO: [2/5] Installing Virsorter2 ..."
    if conda info --envs | grep -q -w "vs2"; then
        echo -e "\tEnv: 'vs2' already exists."
    else
        echo -e "\tEnv: 'vs2' not existed, creating..."
        # create vs2 environment
        conda create -n vs2 -c conda-forge -c bioconda virsorter=2 --yes
        # conda create -n vs2 bioconda::virsorter --yes
        if [ $? -eq 0 ]; then
            echo -e "\tEnv: 'vs2' created."
        fi
    fi
    # GeNomad
    echo "INFO: [3/5] Installing GeNomad ..."
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
    echo "INFO: [4/5] Installing ViraLM ..."
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
    # vcontact3
    echo "INFO: [5/5] Installing vContact3 ..."
    if conda info --envs | grep -q -w "vcontact3"; then
        echo -e "\tEnv: 'vcontact3' already exists."
    else
        echo -e "\tEnv: 'vcontact3' not existed, creating..."
        conda create -n vcontact3 bioconda::vcontact3 --yes
        if [ $? -eq 0 ]; then
            echo -e "\tEnv: 'vcontact3' created."
        fi
    fi
    # tool path check and env variables settings
    if [ -f $WORKING_DIR/src/envs.py ];then
        rm $WORKING_DIR/src/envs.py
    fi
    # conda path check
    if [ -d $CONDA_PATH ];then
        echo "CONDA_PATH = \"$CONDA_PATH\"" >> $WORKING_DIR/src/envs.py
    fi
    echo "INSTALLATION_PATH = \"$WORKING_DIR\"" >> $WORKING_DIR/src/envs.py
    # main env check
    if [ -d $CONDA_ENVS_PATH/$MAIN_ENV_NAME ]; then
        echo -e "MAIN_ENV_NAME = \"${MAIN_ENV_NAME}\"" >> $WORKING_DIR/src/envs.py
    else
        echo "MAIN_ENV_NAME = " >> $WORKING_DIR/src/envs.py
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
    if [ -f $CONDA_ENVS_PATH/$MAIN_ENV_NAME/bin/checkv ]; then
        echo -e "CHECKV_PATH = \"$CONDA_ENVS_PATH/$MAIN_ENV_NAME/bin/checkv\"" >> $WORKING_DIR/src/envs.py
    else 
        echo "CHECKV_PATH = " >> $WORKING_DIR/src/envs.py
    fi
    # strobealign
    if [ -f $CONDA_ENVS_PATH/$MAIN_ENV_NAME/bin/strobealign ]; then
        echo -e "STROBEALIGN_PATH = \"$CONDA_ENVS_PATH/$MAIN_ENV_NAME/bin/strobealign\"" >> $WORKING_DIR/src/envs.py
    else 
        echo "STROBEALIGN_PATH = " >> $WORKING_DIR/src/envs.py
    fi
    # samtools
    if [ -f $CONDA_ENVS_PATH/$MAIN_ENV_NAME/bin/samtools ]; then
        echo -e "SAMTOOLS_PATH = \"$CONDA_ENVS_PATH/$MAIN_ENV_NAME/bin/samtools\"" >> $WORKING_DIR/src/envs.py
    else 
        echo "SAMTOOLS_PATH = " >> $WORKING_DIR/src/envs.py
    fi
    # featureCounts
    if [ -f $CONDA_ENVS_PATH/$MAIN_ENV_NAME/bin/featureCounts ]; then
        echo -e "FEATURECOUNTS_PATH = \"$CONDA_ENVS_PATH/$MAIN_ENV_NAME/bin/featureCounts\"" >> $WORKING_DIR/src/envs.py
    else 
        echo "FEATURECOUNTS_PATH = " >> $WORKING_DIR/src/envs.py
    fi
    # vContact3 
    if [ -f $CONDA_ENVS_PATH/$MAIN_ENV_NAME/bin/featureCounts ]; then
        echo -e "VCONTACT3_PATH = \"$CONDA_ENVS_PATH/vcontact3/bin/vcontact3\"" >> $WORKING_DIR/src/envs.py
    else 
        echo "VCONTACT3_PATH = " >> $WORKING_DIR/src/envs.py
    fi
}

prepare_databases() {
    # checkv
    echo -e "Preparing checkv db ..."
    if [ ! -d $WORKING_DIR/dependencies/checkvdb ]; then
        mkdir $WORKING_DIR/dependencies/checkvdb
        $CONDA_ENVS_PATH/$MAIN_ENV_NAME/bin/checkv download_database $WORKING_DIR/dependencies/checkvdb/
    else
        echo -e "\tcheckv db exists."
    fi
    # CAT_pack nr
    echo -e "Preparing CAT_pack_nr_db ..."
    if [ ! -d $WORKING_DIR/dependencies/CAT_pack_nr_db ]; then
        mkdir $WORKING_DIR/dependencies/CAT_pack_nr_db
        cd $WORKING_DIR/dependencies/CAT_pack_nr_db
        wget -c tbb.bio.uu.nl/tina/CAT_pack_prepare/20240422_CAT_nr.tar.gz
        tar -xvzf 20240422_CAT_nr.tar.gz
    else
        echo -e "\tCAT_pack_nr_db exists."
    fi
    # vs2 db
    echo -e "Preparing vs2_db ..."
    if [ ! -d $WORKING_DIR/dependencies/vs2_db ]; then
        if conda info --envs | grep -q -w "vs2"; then
            cd $WORKING_DIR/dependencies
            source $CONDA_PATH/bin/activate vs2 
            virsorter setup -d vs2_db -j 4
            conda deactivate
        else
            echo -e "[WARNING]: Conda environment: 'vs2' not existed, vs2_db installation failed."
        fi
    else
        echo -e "\tvs2_db exists."
    fi
    # genomad db
    echo -e "Preparing genomad db ..."
    if [ ! -d $WORKING_DIR/dependencies/genomad_db ]; then
        if conda info --envs | grep -q -w "genomad"; then
            source $CONDA_PATH/bin/activate genomad
            genomad download-database $WORKING_DIR/dependencies/
            conda deactivate
        else
            echo -e "[WARNING]: Conda environment: 'genomad' not existed, genomad db installation failed."
        fi
    else
        echo -e "\tgenomad db exists."
    fi
    # viralm model
    echo -e "Preparing viralm model ..."
    if [ ! -d $WORKING_DIR/dependencies/ViraLM/model ]; then
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
    else
        echo -e "\tviralm model exists."
    fi
    # vcontact3 db
    echo -e "Preparing vContact3_db ..."
    if [ ! -d $WORKING_DIR/dependencies/vcontact3_db ]; then
        if conda info --envs | grep -q -w "vcontact3"; then
            cd $WORKING_DIR/dependencies
            mkdir $WORKING_DIR/dependencies/vcontact3_db
            source $CONDA_PATH/bin/activate vcontact3
            # download and setup the model
            vcontact3 prepare_databases --get-version "220" --set-location $WORKING_DIR/dependencies/vcontact3_db
            conda deactivate
        else
            echo -e "[WARNING]: Conda environment: 'vcontact3' not existed, vcontact3_db installation failed."
        fi
    else
        echo -e "\tvContact3_db exists."
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
    if [ -d $WORKING_DIR/dependencies/checkvdb ] && compgen -d $WORKING_DIR/dependencies/checkvdb/ > /dev/null; then
        echo -e "CHECKV_DB_PATH = \"$(compgen -d $WORKING_DIR/dependencies/checkvdb/)\"" >> $WORKING_DIR/src/envs.py
    else 
        echo "CHECKV_DB_PATH = " >> $WORKING_DIR/src/envs.py
    fi
    # vcontact3_db
    if [ -d $WORKING_DIR/dependencies/vcontact3_db ] && compgen -d $WORKING_DIR/dependencies/vcontact3_db/ > /dev/null; then
        echo -e "VCONTACT3_DB_PATH = \"$WORKING_DIR/dependencies/vcontact3_db\"" >> $WORKING_DIR/src/envs.py
        echo -e "VCONTACT3_DB_VERSION = \"220\"" >> $WORKING_DIR/src/envs.py
    else 
        echo "VCONTACT3_DB_PATH = " >> $WORKING_DIR/src/envs.py
        echo "VCONTACT3_DB_VERSION = " >> $WORKING_DIR/src/envs.py
    fi

}

check_envs() {
    # tool path check and env variables settings
    if [ -f $WORKING_DIR/src/envs.py ];then
        rm $WORKING_DIR/src/envs.py
    fi
    # conda path check
    if [ -d $CONDA_PATH ];then
        echo "CONDA_PATH = \"$CONDA_PATH\"" >> $WORKING_DIR/src/envs.py
    fi
    echo "INSTALLATION_PATH = \"$WORKING_DIR\"" >> $WORKING_DIR/src/envs.py
    # main env check
    if [ -d $CONDA_ENVS_PATH/$MAIN_ENV_NAME ]; then
        echo -e "MAIN_ENV_NAME = \"${MAIN_ENV_NAME}\"" >> $WORKING_DIR/src/envs.py
    else
        echo "MAIN_ENV_NAME = " >> $WORKING_DIR/src/envs.py
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
    if [ -f $CONDA_ENVS_PATH/$MAIN_ENV_NAME/bin/checkv ]; then
        echo -e "CHECKV_PATH = \"$CONDA_ENVS_PATH/$MAIN_ENV_NAME/bin/checkv\"" >> $WORKING_DIR/src/envs.py
    else 
        echo "CHECKV_PATH = " >> $WORKING_DIR/src/envs.py
    fi
    # strobealign
    if [ -f $CONDA_ENVS_PATH/$MAIN_ENV_NAME/bin/strobealign ]; then
        echo -e "STROBEALIGN_PATH = \"$CONDA_ENVS_PATH/$MAIN_ENV_NAME/bin/strobealign\"" >> $WORKING_DIR/src/envs.py
    else 
        echo "STROBEALIGN_PATH = " >> $WORKING_DIR/src/envs.py
    fi
    # samtools
    if [ -f $CONDA_ENVS_PATH/$MAIN_ENV_NAME/bin/samtools ]; then
        echo -e "SAMTOOLS_PATH = \"$CONDA_ENVS_PATH/$MAIN_ENV_NAME/bin/samtools\"" >> $WORKING_DIR/src/envs.py
    else 
        echo "SAMTOOLS_PATH = " >> $WORKING_DIR/src/envs.py
    fi
    # featureCounts
    if [ -f $CONDA_ENVS_PATH/$MAIN_ENV_NAME/bin/featureCounts ]; then
        echo -e "FEATURECOUNTS_PATH = \"$CONDA_ENVS_PATH/$MAIN_ENV_NAME/bin/featureCounts\"" >> $WORKING_DIR/src/envs.py
    else 
        echo "FEATURECOUNTS_PATH = " >> $WORKING_DIR/src/envs.py
    fi
    # vContact3 
    if [ -f $CONDA_ENVS_PATH/$MAIN_ENV_NAME/bin/featureCounts ]; then
        echo -e "VCONTACT3_PATH = \"$CONDA_ENVS_PATH/vcontact3/bin/vcontact3\"" >> $WORKING_DIR/src/envs.py
    else 
        echo "VCONTACT3_PATH = " >> $WORKING_DIR/src/envs.py
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
    if [ -d $WORKING_DIR/dependencies/checkvdb ] && compgen -d $WORKING_DIR/dependencies/checkvdb/ > /dev/null; then
        echo -e "CHECKV_DB_PATH = \"$(compgen -d $WORKING_DIR/dependencies/checkvdb/)\"" >> $WORKING_DIR/src/envs.py
    else 
        echo "CHECKV_DB_PATH = " >> $WORKING_DIR/src/envs.py
    fi
    # vcontact3_db
    if [ -d $WORKING_DIR/dependencies/vcontact3_db ] && compgen -d $WORKING_DIR/dependencies/vcontact3_db/ > /dev/null; then
        echo -e "VCONTACT3_DB_PATH = \"$WORKING_DIR/dependencies/vcontact3_db/\"" >> $WORKING_DIR/src/envs.py
        echo -e "VCONTACT3_DB_VERSION = \"220\"" >> $WORKING_DIR/src/envs.py
    else 
        echo "VCONTACT3_DB_PATH = " >> $WORKING_DIR/src/envs.py
        echo "VCONTACT3_DB_VERSION = " >> $WORKING_DIR/src/envs.py
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
    "--check_envs")
        check_envs
        ;;
    *)
        echo -e "options: '--all' '--tools' '--databases' '--check_envs'"
        ;;
esac