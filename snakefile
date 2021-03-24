# directory
workdir: '/projectnb/bf528/users/lava_lamp/project_3/'

# samples
samples = []
with open(('alec_pr3/samples.txt'),'r') as f:
  for line in f:
    samples.append(line.strip('\n'))

# output files 
bams = []
fastqcs = []
with open(('alec_pr3/samples.txt'),'r') as f:
  for line in f:
    bams.append('alec_pr3/bams/'+line.strip('\n')+'Aligned.sortedByCoord.out.bam')
    fastqcs.append('alec_pr3/fastqc/'+line.strip('\n')+'_1_fastqc.html')
    fastqcs.append('alec_pr3/fastqc/'+line.strip('\n')+'_2_fastqc.html')


rule All: 
  input: 
    bams,
    fastqcs

rule fastqc:
  input: 
    'alec_pr3/raw_reads/{sample}_{number}.fastq.gz'
  output: 
    'alec_pr3/fastqc/{sample}_{number}.fastqc'
  shell: 
   'fastqc {input} --outdir=./alec_pr3/fastqc '

rule STAR:
  input: 
    reads_1 = ('alec_pr3/raw_reads/{sample}_1.fastq.gz'),
    reads_2 = ('alec_pr3/raw_reads/{sample}_2.fastq.gz'),
    genome = '/project/bf528/project_3/reference/rn4_STAR'
  shell: 
    'STAR --genomeDir {input.genome} --runThreadN 16 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix alec_pr3/bams/{wildcards.sample}  --readFilesIn {input.reads_1} {input.reads_2}'


