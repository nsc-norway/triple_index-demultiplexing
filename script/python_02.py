import sys

SAMPLE_NAME = sys.argv[1]
HISEQ_NAME = sys.argv[2]

out_file = open('bash_python_02.sh', 'w')

print >> out_file, 'mkdir 20_demultiplex'
print >> out_file, 'mkdir 20_demultiplex/10_cutadapt'
print >> out_file, 'mkdir 20_demultiplex/10_cutadapt/log'
print >> out_file, ''
    
import os

pwd = os.getcwd()

PCR1_fwd = '../tools/index/PCR1_fwd'
PCR1_rev = '../tools/index/PCR1_rev'

read1 = '10_preprocess/20_bbmap/' + SAMPLE_NAME + '-NonPhiX_R1.fq '
read2 = '10_preprocess/20_bbmap/' + SAMPLE_NAME + '-NonPhiX_R2.fq '

for item in open(PCR1_fwd):
    id = item.rsplit()[0]
    id_seq = item.rsplit()[1]
    cutadapt_string1 = 'cutadapt --discard-untrimmed --no-trim -g \'^'
    output_file = '20_demultiplex/10_cutadapt/' + SAMPLE_NAME + '_R1_' + id + '.fastq.gz'
    output_log = '20_demultiplex/10_cutadapt/log/cutadapt_' + SAMPLE_NAME + '_R1_' + id + '.log'
   
    cutadapt_process = cutadapt_string1 + id_seq + '\' -o ' + output_file + ' ' + read1 + ' > ' + output_log
    print >> out_file, cutadapt_process
    
print >> out_file, ''

for item in open(PCR1_rev):
    id = item.rsplit()[0]
    id_seq = item.rsplit()[1]
    cutadapt_string1 = 'cutadapt --discard-untrimmed --no-trim -g \'^'
    output_file = '20_demultiplex/10_cutadapt/' + SAMPLE_NAME + '_R2_' + id + '.fastq.gz'
    output_log = '20_demultiplex/10_cutadapt/log/cutadapt_' + SAMPLE_NAME + '_R2_' + id + '.log'
   
    cutadapt_process = cutadapt_string1 + id_seq + '\' -o ' + output_file + ' ' + read2 + ' > ' + output_log
    print >> out_file, cutadapt_process
    
out_file.close()
