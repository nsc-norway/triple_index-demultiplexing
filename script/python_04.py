import sys

SAMPLE_NAME = sys.argv[1]
HISEQ_NAME = sys.argv[2]

out_file = open('bash_python_04.sh', 'w')

print >> out_file, 'mkdir 20_demultiplex/30_split'
print >> out_file, ''
    
import os

pwd = os.getcwd()

PCR1_fwd = 'tools/index/PCR1_fwd'
PCR1_rev = 'tools/index/PCR1_rev'

read1 = '10_preprocess/10_bbmap/' + SAMPLE_NAME + '-clean_R1.fq.gz '
read2 = '10_preprocess/10_bbmap/' + SAMPLE_NAME + '-clean_R2.fq.gz '
  
for item in os.listdir('20_demultiplex/20_header/common/'):
    if SAMPLE_NAME + '_' in item:
        item_path = '20_demultiplex/20_header/common/' + item
        out_path = '20_demultiplex/30_split/' + item

        print >> out_file, 'grep -A3 -Fxf ' + item_path + ' <(gunzip -c ' + read1 + ') | grep -v "^\-\-$" > ' + out_path + '_R1.fq'
        print >> out_file, 'grep -A3 -Fxf <(sed "s/1:N:0/2:N:0/g" ' + item_path + ') <(gunzip -c ' + read2 + ') | grep -v "^\-\-$" > ' + out_path + '_R2.fq'
        print >> out_file, ''
        
out_file.close()
