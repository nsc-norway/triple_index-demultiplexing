### Manual Prep ##

### mkdir 00_data
### mkdir 00_data/00_raw

### Copy sample data into Sample_Folder/00_data/00_raw
### make sure 'cutadapt' is available from bash

TRIMMOMATIC=/data/extraProjects/arvindsu/Pal_3/tools/Trimmomatic-0.36/trimmomatic-0.36.jar
TRIMMOMATIC_ADAPTOR=/data/extraProjects/arvindsu/Pal_3/tools/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa

BBMAP=/data/extraProjects/arvindsu/Pal_3/tools/bbmap/bbmap.sh
PHIX=/data/extraProjects/arvindsu/Pal_3/tools/PhiX/ToRemove.fa

### fill information above till here

# Run this script as below
# bash script/bash_01_startscript.sh 1-P1 @D00132 &

SAMPLE_NAME=$1
HISEQ_NAME=$2

echo
echo "Sample name provided: " $SAMPLE_NAME
echo "HiSeq name provided: " $HISEQ_NAME

mkdir 00_data/10_cat
cd 00_data/10_cat
cat ../00_raw/$SAMPLE_NAME*R1* > $SAMPLE_NAME-R1.fastq.gz
cat ../00_raw/$SAMPLE_NAME*R2* > $SAMPLE_NAME-R2.fastq.gz
cd ../../

echo "Completed cat. Starting Trimmomatic"

mkdir 10_preprocess
mkdir 10_preprocess/10_trimmomatic
cd 10_preprocess/10_trimmomatic
mkdir log

java -Xmx10G -jar $TRIMMOMATIC PE -phred33 -threads 2 -trimlog log/trim.log ../../00_data/10_cat/$SAMPLE_NAME-R1.fastq.gz ../../00_data/10_cat/$SAMPLE_NAME-R2.fastq.gz $SAMPLE_NAME-trim_R1.fq $SAMPLE_NAME-trim_R1_un.fq $SAMPLE_NAME-trim_R2.fq $SAMPLE_NAME-trim_R2_un.fq ILLUMINACLIP:$TRIMMOMATIC_ADAPTOR:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 CROP:250 2> log/trim.log2

cd log
rm trim.log
cd ..

echo "Completed Trimmomatic. Starting bbmap"

cd ../../
mkdir 10_preprocess/20_bbmap
cd 10_preprocess/20_bbmap
mkdir log

$BBMAP -Xmx10g ref=$PHIX in=../10_trimmomatic/$SAMPLE_NAME-trim_R1.fq in2=../10_trimmomatic/$SAMPLE_NAME-trim_R2.fq outu=$SAMPLE_NAME-NonPhiX_R1.fq outu2=$SAMPLE_NAME-NonPhiX_R2.fq 2> log/bbmap.out2

rm -r ref

cd ../../

echo "Completed bbmap. Starting python 02"

python ../script/python_02.py $SAMPLE_NAME $HISEQ_NAME
bash bash_python_02.sh

echo "Completed python 02. Starting python 03"

python ../script/python_03.py $SAMPLE_NAME $HISEQ_NAME
bash bash_python_03.sh

echo "Completed python 03. Starting python 04"

python ../script/python_04.py $SAMPLE_NAME $HISEQ_NAME
bash bash_python_04.sh

echo "Completed python 04. Starting python 05"

python ../script/python_05.py $SAMPLE_NAME $HISEQ_NAME
bash bash_python_05.sh

echo "Completed python 05. All done"
