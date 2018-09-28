import sys

SAMPLE_NAME = sys.argv[1]
HISEQ_NAME = sys.argv[2]

out_file = open('bash_python_03.sh', 'w')

print >> out_file, 'mkdir 20_demultiplex/20_header'
print >> out_file, 'mkdir 20_demultiplex/20_header/common'
print >> out_file, ''
    
import os

pwd = os.getcwd()

PCR1_fwd = 'tools/index/PCR1_fwd'
PCR1_rev = 'tools/index/PCR1_rev'

read1 = '10_preprocess/10_bbmap/' + SAMPLE_NAME + '-NonPhiX_R1.fq '
read2 = '10_preprocess/10_bbmap/' + SAMPLE_NAME + '-NonPhiX_R2.fq '
  
for item in os.listdir('20_demultiplex/10_cutadapt'):
    if '.fastq.gz' in item:
        header_process = 'gunzip -c 20_demultiplex/10_cutadapt/' + item + ' | grep ' + HISEQ_NAME + ' | sed "s/[12]:N:0/1:N:0/g" |sort > 20_demultiplex/20_header/' + item.replace('.fastq.gz','.header')
        print >> out_file, header_process

print >> out_file, '' 
print >> out_file, 'cd 20_demultiplex/20_header'
print >> out_file, ''
    
for item in open(PCR1_fwd):
    for item_r in open(PCR1_rev):
        print >> out_file, 'comm -12 ' + SAMPLE_NAME + '_R1_' + item.split()[0] + '.header ' + SAMPLE_NAME + '_R2_' + item_r.split()[0] + '.header > common/' + SAMPLE_NAME + '_R1_' + item.split()[0] + '-' + SAMPLE_NAME + '_R2_' + item_r.split()[0] 
    

print >> out_file, ''
print >> out_file, 'cd ../../'
print >> out_file, ''

out_file.close()
