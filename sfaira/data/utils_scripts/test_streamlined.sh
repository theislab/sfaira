#!/bin/bash

CODE_PATH="/home/icb/${USER}/git"
OUT_PATH="/storage/groups/ml01/workspace/david.fischer/sfaira/cellxgene/processed_data/"
OUT_FN="${OUT_PATH}validation_summary.txt"
SCHEMA="cellxgene"
DOIS="10.1016/j.cell.2019.08.008,10.1016/j.celrep.2018.11.086,10.1016/j.cmet.2016.08.020,10.1016/j.neuron.2019.06.011,10.1038/s41422-018-0099-2,10.1038/s41467-018-06318-7,10.1038/s41590-020-0602-z,10.1084/jem.20191130,10.1101/2020.03.13.991455,10.1101/2020.10.12.335331,10.1126/science.aay3224,10.1126/science.aba7721,10.15252/embj.2018100811,no_doi_10x_genomics"

source "/home/${USER}/.bashrc"
echo "Summary of exports of data sets ${DOIS}" > OUT_FN
for doi in "${DOIS[@]}"; do
  echo "Summary of exports of data set ${doi}" >> OUT_FN
  cellxgene schema validate ${OUT_PATH}${doi}/ >> OUT_FN
done

CODE_PATH="/home/icb/${USER}/git"
OUT_PATH="/storage/groups/ml01/workspace/david.fischer/sfaira/cellxgene/processed_data/"
OUT_FN="${OUT_PATH}validation_summary.txt"
SCHEMA="cellxgene"
DOIS="10.1016/j.cell.2017.09.004,10.1016/j.cell.2018.02.001,10.1016/j.cell.2018.08.067,10.1016/j.cell.2019.06.029,10.1016/j.cels.2016.08.011,10.1016/j.cmet.2019.01.021,10.1016/j.devcel.2020.01.033,10.1038/nmeth.4407,10.1038/s41467-019-10861-2,10.1038/s41467-019-12464-3,10.1038/s41467-019-12780-8,10.1038/s41586-018-0698-6,10.1038/s41586-019-1373-2,10.1038/s41586-019-1631-3,10.1038/s41586-019-1652-y,10.1038/s41586-019-1654-9,10.1038/s41586-020-2157-4,10.1038/s41586-020-2922-4,10.1038/s41591-019-0468-5,10.1038/s41593-019-0393-4,10.1038/s41597-019-0351-8,10.1073/pnas.1914143116,10.1101/661728,10.1101/753806,10.1126/science.aat5031,10.1186/s13059-019-1906-x,no_doi_regev"

source "/home/${USER}/.bashrc"
python ${CODE_PATH}/sfaira/sfaira/data/utils_scripts/streamline_selected.py ${DATA_PATH} ${META_PATH} ${CACHE_PATH} ${OUT_PATH} ${SCHEMA} ${DOIS}
