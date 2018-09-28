### Manual Prep ##

### mkdir 00_data
### mkdir 00_data/00_raw

### Copy sample data into Sample_Folder/00_data/00_raw
### make sure 'cutadapt' is available from bash

TOOLS=$PWD/tools

BBMAP=$TOOLS/bbmap/bbduk.sh
PHIX=$TOOLS/PhiX/adapters_bbduk.fa

### fill information above till here

# Run this script as below
# bash script/bash_01_startscript.sh 1-P1 @D00132 &

######## SAMPLE_NAME=$1
SAMPLE_NAME="1-Trosvik-P1"
######## MISEQ_NAME=$2
MISEQ_NAME="@M02980"
echo
echo "Sample name provided: " $SAMPLE_NAME
echo "MiSeq name provided: " $MISEQ_NAME

######## cd 00_data/00_raw/
######## ln -s  $SAMPLE_NAME*R1* $SAMPLE_NAME-R1.fq.gz
######## ln -s  $SAMPLE_NAME*R2* $SAMPLE_NAME-R2.fq.gz
######## 
######## cd ../../
######## 
######## echo "Starting bbmap"
######## 
######## mkdir -p 10_preprocess/10_bbmap
######## cd 10_preprocess/10_bbmap
######## mkdir log
######## 
######## $BBMAP -Xmx10g ref=$PHIX ktrim=r k=23 mink=11 hdist=1 tbo tpe qtrim=r trimq=15 maq=15 minlen=36 forcetrimright=300 in=../../00_data/00_raw/$SAMPLE_NAME-R1.fq.gz in2=../../00_data/00_raw/$SAMPLE_NAME-R2.fq.gz outu=$SAMPLE_NAME-clean_R1.fq.gz outu2=$SAMPLE_NAME-clean_R2.fq.gz 2> log/bbmap.out2
######## 
######## cd ../../
########
######## echo "Completed bbmap. Starting python 02"
########
######## python script/python_02.py $SAMPLE_NAME $MISEQ_NAME
######## bash bash_python_02.sh
########
######## echo "Completed python 02. Starting python 03"
######## python script/python_03.py $SAMPLE_NAME $MISEQ_NAME
######## bash bash_python_03.sh
######## 
######## echo "Completed python 03. Starting python 04"
######## 
######## python script/python_04.py $SAMPLE_NAME $MISEQ_NAME
######## bash bash_python_04.sh
######## 
######## echo "Completed python 04. Starting python 05"
######## 
python script/python_05.py $SAMPLE_NAME $MISEQ_NAME
######## bash bash_python_05.sh
######## 
######## echo "Completed python 05. All done"
