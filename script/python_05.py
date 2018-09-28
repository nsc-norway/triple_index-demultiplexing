import sys

SAMPLE_NAME = sys.argv[1]
HISEQ_NAME = sys.argv[2]

out_file = open('bash_python_05.sh', 'w')

print >> out_file, 'mkdir 20_demultiplex/40_cutadapt'
print >> out_file, 'mkdir 20_demultiplex/40_cutadapt/log'
print >> out_file, ''

import os

pwd = os.getcwd()

PCR1_fwd_spacer = {}
for item in open('tools/index/PCR1_fwd_spacer', 'r'):
    PCR1_fwd_spacer[item.rstrip().split()[0]] = item.rstrip().split()[1]
    
PCR1_rev_spacer = {}
for item in open('tools/index/PCR1_rev_spacer', 'r'):
    PCR1_rev_spacer[item.rstrip().split()[0]] = item.rstrip().split()[1]

cutadapt_string1 = 'tools/cutadapt/bin/cutadapt -g "^'

for item in os.listdir('20_demultiplex/30_split/'):
    if 'R1.fq' in item:
        spacer = PCR1_fwd_spacer[item.replace(SAMPLE_NAME, '').split('-')[0].split('_')[2]]
        input_file = '20_demultiplex/30_split/' + item
        output_file = '20_demultiplex/40_cutadapt/' + item.replace('.fq','_cutadapt.fastq.gz')
        output_log = '20_demultiplex/40_cutadapt/log/cutadapt_' + item.replace('.fq','.log')
        
        print >> out_file, cutadapt_string1 + spacer + '" -o ' + output_file + ' ' + input_file + ' > ' + output_log 
        
    if 'R2.fq' in item:
        spacer = PCR1_rev_spacer[item.replace(SAMPLE_NAME, '').split('-')[1].split('_')[2]]
        input_file = '20_demultiplex/30_split/' + item
        output_file = '20_demultiplex/40_cutadapt/' + item.replace('.fq','_cutadapt.fastq.gz')
        output_log = '20_demultiplex/40_cutadapt/log/cutadapt_' + item.replace('.fq','.log')
        
        print >> out_file, cutadapt_string1 + spacer + '" -o ' + output_file + ' ' + input_file + ' > ' + output_log
            
out_file.close()
