import yaml

configfile: "config.yml"

path_to_data = config['data']['base_path']

rule all:
    input:
        expand(
            path_to_data + '/samples/{sample}/wig_file/cfmedip.hg19.wig',
            sample = glob_wildcards('/cluster/projects/scottgroup/people/eric/projects/cfMeDIPseq_pipeline_test/data/input/fastq_files/{sample}_R1_001.fastq.gz').sample
        )

rule tag_to_header:
    input:
        R1=path_to_data + '/input/fastq_files/{sample}_R1_001.fastq.gz',
        R2=path_to_data + '/input/fastq_files/{sample}_R2_001.fastq.gz',
    output:
        path_to_data + '/samples/{sample}/fastq_tag_to_header/data.seq1.smi.fq.gz',
        path_to_data + '/samples/{sample}/fastq_tag_to_header/data.seq2.smi.fq.gz'
    params:
        outdir=lambda wildcards, output: '/'.join(output[0].split('/')[0:-1]) + '/data'
        # tag_to_header.py requires you to specify an output prefix, and appends seq[#].smi.fq.gz to it
    resources: cpus=1, mem_mb=8000, time_min='72:00:00'
    shell:
        'python src/external/tag_to_header.py --infile1 {input.R1} --infile2 {input.R2} --outprefix {params.outdir} --taglen 2 --spacerlen 1'

rule bwa_mem:
    input:
        path_to_data + '/samples/{sample}/fastq_tag_to_header/data.seq1.smi.fq.gz',
        path_to_data + '/samples/{sample}/fastq_tag_to_header/data.seq2.smi.fq.gz'
    output:
        path_to_data + '/samples/{sample}/bwa_mem/aligned.hg19.sam'
    resources: cpus=4, mem_mb=16000, time_min='72:00:00'
    shell:
        r"bwa mem -M -t4 -R'@RG\tID:1\tSM:{wildcards.sample}\tPL:illumina\tPU:.\tLB:NA' $BWAINDEX {input} > {output}"

rule sam_to_sorted_bam:
    input:
        path_to_data + '/samples/{sample}/bwa_mem/aligned.hg19.sam'
    output:
        bam=path_to_data + '/samples/{sample}/bwa_mem/aligned.hg19.sorted.bam',
        index=path_to_data + '/samples/{sample}/bwa_mem/aligned.hg19.sorted.bam.bai'
    resources: cpus=4, mem_mb=8000, time_min='72:00:00'
    shell:
        "samtools view -buS -f 2 -F 4 -@4 {input} | samtools fixmate -m - - | samtools sort -@4 -o {output.bam} && samtools index {output.bam}"

rule bam_markdup:
    input:
        path_to_data + '/samples/{sample}/bwa_mem/aligned.hg19.sorted.bam',
    output:
        bam=path_to_data + '/samples/{sample}/bwa_mem/aligned.hg19.sorted.markdup.bam',
        index=path_to_data + '/samples/{sample}/bwa_mem/aligned.hg19.sorted.markdup.bam.bai',
    resources: mem_mb=8000, time_min='72:00:00'
    shell:
        "samtools markdup -r {input} {output.bam} && samtools index {output.bam}"

rule bam_to_wig:
    input:
        path_to_data + '/samples/{sample}/bwa_mem/aligned.hg19.sorted.markdup.bam',
    output:
        path_to_data + '/samples/{sample}/wig_file/cfmedip.hg19.wig',
    resources: mem_mb=8000, time_min='72:00:00'
    shell:
        "Rscript src/bam_to_wig.R -b {input} -o {output}"

