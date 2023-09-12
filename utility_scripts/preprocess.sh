PATH_TO_DATA=$1
DATA_NAME=$2
R_LIB_PATH=$3

datapath=${PATH_TO_DATA}
dataname=${DATA_NAME}
rlibpath=${R_LIB_PATH}

if ($dataname=='zheng'); then

python preprocess_tutorial_data.py --path_to_data tutorial_data/zheng2017/ --path_to_save tutorial_data/zheng2017/5000hvgs_zheng2017_preprocessed.h5ad --data_name "zheng" --n_hvgs 5000

elif ($dataname=='vanGalen'); then

Rscript utility_scripts/read_rds.R ${rlibpath}

python preprocess_tutorial_data.py --path_to_data tutorial_data/vanGalen2019/ --path_to_save tutorial_data/vanGalen2019/5000hvgs_vanGalen2019_preprocessed.h5ad --data_name "vanGalen" --n_hvgs 5000

elif ($dataname=='dominguezconde'); then

python preprocess_tutorial_data.py --path_to_data tutorial_data/dominguezconde2022/ --path_to_save tutorial_data/dominguezconde2022/5000hvgs_dominguezconde2022_preprocessed.h5ad --data_name "dominguezconde" --n_hvgs 5000

elif ($dataname=='raghavan'); then

...
fi