import yaml
import pandas as pd

configfile: "config.yml"

path_to_data = config['data']['base_path']

def get_cohort_data(cohort):
    return(pd.read_csv(config['data']['cohorts'][cohort], comment='#'))

def get_fastq_path(sample, read_in_pair=1):
    sample_line = all_samples[
        (all_samples.sample_name == sample) &
        (all_samples.read_in_pair == read_in_pair)
    ]
    return(sample_line.path.to_list()[0])

all_samples = pd.concat([get_cohort_data(c) for c in config['data']['cohorts']])

def clean(command):
    command = command.replace('\n', ' ').replace('\t', ' ')
    while '  ' in command:
        command = command.replace('  ', ' ')
    return(command.strip())

rule all:
    input:
        expand(
            path_to_data + '/samples/{sample}/bin_metrics/{filename}',
            sample = all_samples['sample_name'].unique().tolist(),
            filename = ['COUNTwiggle.wig.txt', 'medestrand.binmethyl.txt']
        )

rule tag_to_header:
    input:
        R1=lambda wildcards: get_fastq_path(wildcards.sample, 1),
        R2=lambda wildcards: get_fastq_path(wildcards.sample, 2)
    output:
        path_to_data + '/samples/{sample}/fastq_tag_to_header/data.seq1.smi.fq.gz',
        path_to_data + '/samples/{sample}/fastq_tag_to_header/data.seq2.smi.fq.gz'
    params:
        outdir=lambda wildcards, output: '/'.join(output[0].split('/')[0:-1]) + '/data'
        # tag_to_header.py requires you to specify an output prefix, and appends seq[#].smi.fq.gz to it
    resources: cpus=1, mem_mb=8000, time_min='72:00:00'
    shell:
        clean(r'''
        module load picard &&
        python src/external/tag_to_header.py --infile1 {input.R1} --infile2 {input.R2} --outprefix {params.outdir} --taglen 2 --spacerlen 1
        ''')

rule bwa_mem:
    input:
        path_to_data + '/samples/{sample}/fastq_tag_to_header/data.seq1.smi.fq.gz',
        path_to_data + '/samples/{sample}/fastq_tag_to_header/data.seq2.smi.fq.gz'
    output:
        path_to_data + '/samples/{sample}/bwa_mem/aligned.hg19.sam'
    resources: cpus=4, mem_mb=16000, time_min='72:00:00'
    shell:
        clean(r"""
        module load picard &&
        module load bwa/0.7.15 &&
        module load igenome-human/hg19 &&
        module load samtools/1.3.1 &&
        bwa mem -M -t4 -R'@RG\tID:1\tSM:{wildcards.sample}\tPL:illumina\tPU:.\tLB:NA' $BWAINDEX {input} > {output}
        """)

rule sam_to_sorted_bam:
    input:
        path_to_data + '/samples/{sample}/bwa_mem/aligned.hg19.sam'
    output:
        bam=path_to_data + '/samples/{sample}/bwa_mem/aligned.hg19.sorted.bam',
        index=path_to_data + '/samples/{sample}/bwa_mem/aligned.hg19.sorted.bam.bai'
    resources: cpus=4, mem_mb=8000, time_min='72:00:00'
    shell:
        # Try it without fixmate for replication purposes
        #"samtools view -buS -f 2 -F 4 -@4 {input} | samtools fixmate -m - - | samtools sort -@4 -o {output.bam} && samtools index {output.bam}"
        clean(r'''
        module load picard &&
        module load bwa/0.7.15 &&
        module load igenome-human/hg19 &&
        module load samtools/1.3.1 &&
        samtools view -buS -f 2 -F 4 -@4 {input} | samtools sort -@4 -o {output.bam} && samtools index {output.bam}
        ''')

rule bam_markdup:
    input:
        path_to_data + '/samples/{sample}/bwa_mem/aligned.hg19.sorted.bam',
    output:
        bam=path_to_data + '/samples/{sample}/bwa_mem/aligned.hg19.sorted.markdup.bam',
        index=path_to_data + '/samples/{sample}/bwa_mem/aligned.hg19.sorted.markdup.bam.bai',
    resources: mem_mb=8000, time_min='72:00:00'
    shell:
        # Try it with rmdup for replication purposes
        #"samtools markdup -r {input} {output.bam} && samtools index {output.bam}"
        clean(r'''
        module load picard &&
        module load bwa/0.7.15 &&
        module load igenome-human/hg19 &&
        module load samtools/1.3.1 &&
        samtools rmdup {input} {output.bam} && samtools index {output.bam}
        ''')

rule bam_to_wig:
    input:
        path_to_data + '/samples/{sample}/bwa_mem/aligned.hg19.sorted.markdup.bam',
    output:
        path_to_data + '/samples/{sample}/bin_metrics/COUNTwiggle.wig.txt',
    resources: mem_mb=30000, time_min='72:00:00'
    shell:
        clean(r'''
        Rscript src/bam_to_wig.R -b {input} -o {output}
        ''')

rule bam_to_binmethyl:
    input:
        path_to_data + '/samples/{sample}/bwa_mem/aligned.hg19.sorted.markdup.bam',
    output:
        path_to_data + '/samples/{sample}/bin_metrics/medestrand.binmethyl.txt',
    resources: mem_mb=30000, time_min='72:00:00'
    shell:
        clean(r'''
        Rscript src/correct_cpg.R -b {input} -o {output}
        ''')

