import yaml
import pandas as pd
import gzip

configfile: "config.yml"

path_to_data = config['data']['base_path']

def get_cohort_data(cohort):
    samplesheet = pd.read_csv(config['data']['cohorts'][cohort]['samplesheet'], comment='#').drop_duplicates()
    samplesheet = samplesheet.sort_values(by=['sample_name', 'read_in_pair'])
    return(samplesheet)

def get_fastq_path(sample, library, read_in_pair=1):
    library = int(library)
    sample_line = all_samples[
        (all_samples.sample_name == sample) &
        (all_samples.library_index == library) &
        (all_samples.read_in_pair == read_in_pair)
    ]
    return(sample_line.path.to_list()[0])

all_samples = pd.concat([
    get_cohort_data(cohort_name)
    for cohort_name
    in config['data']['cohorts']
    if config['data']['cohorts'][cohort_name]['active']
    ])

all_samples_list = all_samples.sample_name.tolist()

def clean(command):
    command = command.replace('\n', ' ').replace('\t', ' ')
    while '  ' in command:
        command = command.replace('  ', ' ')
    return(command.strip())

rule all:
    input:
        #expand(path_to_data + '/samples/{sample}/merged/bin_stats/bin_stats.tsv', sample=all_samples.sample_name.unique()),
        #expand('/cluster/projects/pughlab/projects/ezhao/pipelines/cfmedipseq_pipeline/samples/CMP-01-02-cfDNA-02/merged/bin_stats/bin_stats_fit_{method}.Rds', method=['LBFGS', 'VB', 'MCMC'])
        expand(
            path_to_data + '/samples/{sample}/merged/bin_stats/bin_stats_fit_nbglm.tsv',
            sample = all_samples_list
        ),
        expand(
            path_to_data + '/samples/{sample}/merged/bin_stats/bin_stats_EPIC_positions.tsv',
            sample = all_samples_list
        ),

rule gunzip_fastq:
    input:
        lambda wildcards: get_fastq_path(wildcards.sample, int(wildcards.lib), int(wildcards.read)),
    output:
        path_to_data + '/samples/{sample}/libraries/lib_{lib}/fastq/R{read}.fastq',
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    shell:
        'gunzip -dc {input} > {output}'

rule extract_barcodes:
    input:
        R1 =    path_to_data + '/samples/{sample}/libraries/lib_{lib}/fastq/R1.fastq',
        R2 =    path_to_data + '/samples/{sample}/libraries/lib_{lib}/fastq/R2.fastq',
    output:
        R1 =    path_to_data + '/samples/{sample}/libraries/lib_{lib}/fastq/extract_barcode_R1.fastq',
        R2 =    path_to_data + '/samples/{sample}/libraries/lib_{lib}/fastq/extract_barcode_R2.fastq'
    params:
        outprefix = lambda wildcards, output: output.R1.split('_barcode_')[0]
    resources: cpus=1, mem_mb=16000, time_min='5-00:00:00'
    shell:
        clean(r'''
        python src/ConsensusCruncher/ConsensusCruncher/extract_barcodes.py
            --read1 {input.R1}
            --read2 {input.R2}
            --outfile {params.outprefix}
            --blist ''' + config['paths']['barcodes']
        )


def get_read_group_from_fastq(fastq_file, sample_name):
    with gzip.open(fastq_file, 'rt') as fastq:
        header = next(fastq)
        (instrument, run_number, flowcell, lane, tile, xpos, ypos) = header.split(' ')[0].split(':')
        lib_value = header.strip().split(' ')[1].split(':')[3]
        rg_line = r"@RG\tID:{flowcell}_{lane}\tSM:{sample}\tPL:Illumina\tPU:.\tLB:{lib_value}".format(
            flowcell = flowcell,
            lane = lane,
            sample = sample_name,
            lib_value = lib_value
        )
        return(rg_line)

rule bwa_mem:
    input:
        path_to_data + '/samples/{sample}/libraries/lib_{lib}/fastq/extract_barcode_R1.fastq',
        path_to_data + '/samples/{sample}/libraries/lib_{lib}/fastq/extract_barcode_R2.fastq',
    output:
        path_to_data + '/samples/{sample}/libraries/lib_{lib}/bwa_mem/aligned.sam'
    resources: cpus=4, mem_mb=16000, time_min='72:00:00'
    params:
        read_group = lambda wildcards, input: get_read_group_from_fastq(
            fastq_file = get_fastq_path(wildcards.sample, wildcards.lib),
            sample_name = wildcards.sample
        )
    shell:
        clean(r"""
        bwa mem -M -t4
        -R'{{params.read_group}}' 
        {bwa_index}
        {{input}} > {{output}}""".format(
            bwa_index = config['paths']['bwa_index']
        ))

rule sam_to_sorted_bam:
    input:
        path_to_data + '/samples/{sample}/libraries/lib_{lib}/bwa_mem/aligned.sam'
    output:
        bam=path_to_data + '/samples/{sample}/libraries/lib_{lib}/bwa_mem/aligned.sorted.bam',
        index=path_to_data + '/samples/{sample}/libraries/lib_{lib}/bwa_mem/aligned.sorted.bam.bai'
    resources: cpus=32, mem_mb=30000, time_min='72:00:00'
    shell:
        # Try it without fixmate for replication purposes
        #"samtools view -buS -f 2 -F 4 -@4 {input} | samtools fixmate -m - - | samtools sort -@4 -o {output.bam} && samtools index {output.bam}"
        clean(r'''
        samtools view -buS -f 2 -F 4 -@32 {input} |
        samtools fixmate -m - - |
        samtools sort -@32 -o {output.bam} && samtools index {output.bam}
        ''')

def get_libraries_of_sample(sample):
    filtered_table = all_samples[all_samples.sample_name == sample]
    return(list(set(filtered_table.library_index.to_list())))

rule merge_bam:
    input:
        lambda wildcards: expand(
                path_to_data + '/samples/' + wildcards.sample + '/libraries/lib_{lib}/bwa_mem/aligned.sorted.bam',
                lib=get_libraries_of_sample(wildcards.sample)
        )
    output:
        path_to_data + '/samples/{sample}/merged/bwa_mem/aligned.sorted.bam'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    shell:
        'samtools merge {output} {input} && samtools index {output}'

rule bam_markdup:
    input:
        path_to_data + '/samples/{sample}/merged/bwa_mem/aligned.sorted.bam'
    output:
        bam=path_to_data + '/samples/{sample}/merged/bwa_mem/aligned.sorted.markdup.bam',
        index=path_to_data + '/samples/{sample}/merged/bwa_mem/aligned.sorted.markdup.bam.bai'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    resources: mem_mb=8000, time_min='72:00:00'
    shell:
        "samtools markdup -r {input} {output.bam} && samtools index {output.bam}"

rule bam_to_wig:
    input:
        path_to_data + '/samples/{sample}/bwa_mem/aligned.sorted.markdup.bam',
    output:
        path_to_data + '/samples/{sample}/bin_metrics/COUNTwiggle.wig.txt',
    resources: mem_mb=30000, time_min='72:00:00'
    shell:
        clean(r'''
        Rscript src/bam_to_wig.R -b {input} -o {output}
        ''')

rule bam_to_binmethyl:
    input:
        path_to_data + '/samples/{sample}/bwa_mem/aligned.sorted.markdup.bam',
    output:
        path_to_data + '/samples/{sample}/bin_metrics/medestrand.binmethyl.txt',
    resources: mem_mb=30000, time_min='72:00:00'
    shell:
        clean(r'''
        Rscript src/correct_cpg.R -b {input} -o {output}
        ''')

rule binmethyl_to_wig:
    input:
        binmethyl=path_to_data + '/samples/{sample}/bin_metrics/medestrand.binmethyl.txt',
        wig=path_to_data + '/samples/{sample}/bin_metrics/COUNTwiggle.wig.txt'
    output:
        path_to_data + '/samples/{sample}/bin_metrics/medestrand.binmethyl.wig',
    shell:
        'Rscript src/binmethyl_to_wig.R -b {input.binmethyl} -m {input.wig} -o {output}'

def get_bsgenome_chrom(species, chrom):
    chrom_map = {'Arabidopsis1': 'Chr1', 'Arabidopsis3': 'Chr3'}
    if species == 'human':
        return(chrom)
    elif species == 'arabidopsis':
        return(chrom_map[chrom])

rule bam_bin_stats:
    input:
        path_to_data + '/samples/{sample}/merged/bwa_mem/aligned.sorted.markdup.bam',
    output:
        binstat=path_to_data + '/samples/{sample}/merged/bin_stats/by_chromosome/bin_stats_{species}_{chrom}.tsv',
        filtered=path_to_data + '/samples/{sample}/merged/bin_stats/filtered_out/bin_stats_{species}_{chrom}.tsv',
    params:
        bsgenome = lambda wildcards:config['paths']['bsgenome'][wildcards.species],
        bsgenome_chr = lambda wildcards: get_bsgenome_chrom(wildcards.species, wildcards.chrom)
    resources: cpus=1, mem_mb=16000, time_min='24:00:00'
    shell:
        clean('''
        Rscript src/R/bin_stats.R
            -b {input}
            -g {params.bsgenome}
            -c {wildcards.chrom}
            -o {output.binstat}
            --filtered {output.filtered}
            --bsgchr {params.bsgenome_chr}
        ''')

MAJOR_HUMAN_CHROMOSOMES = ['chr' + str(i) for i in range(1, 22)] + ['chrX', 'chrY']
all_chromosome_tuples = [('human', c) for c in MAJOR_HUMAN_CHROMOSOMES] + [('arabidopsis', 'Arabidopsis1'), ('arabidopsis', 'Arabidopsis3')]

rule merge_bin_stats:
    input:
        [path_to_data + '/samples/{{sample}}/merged/bin_stats/by_chromosome/bin_stats_{species}_{chrom}.tsv'.format(species=a[0], chrom=a[1]) for a in all_chromosome_tuples],
    output:
        path_to_data + '/samples/{sample}/merged/bin_stats/bin_stats.tsv'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    run:
        for i, input_file in enumerate(input):
            input_data = pd.read_csv(input_file, delimiter='\t', comment='#')
            if i == 0:
                input_data.to_csv(output[0], header=True, sep='\t', index=False)
            else:
                input_data.to_csv(output[0], header=False, sep='\t', index=False, mode='a')

# Below is the full bayesian approach for fitting, which remains under development
rule fit_bin_stats:
    input:
        path_to_data + '/samples/{sample}/merged/bin_stats/bin_stats.tsv'
    output:
        path_to_data + '/samples/{sample}/merged/bin_stats/bin_stats_fit_{method}.Rds'
    resources: cpus=1, mem_mb=30000, time_min='5-00:00:00'
    shell:
        'Rscript src/R/fit_cpg_bins.R -i {input} -o {output} --method {wildcards.method}'

# This is the currently used negative binomial GLM approach to fitting
rule cfmedip_nbglm:
    input:
        path_to_data + '/samples/{sample}/merged/bin_stats/bin_stats.tsv'
    output:
        fit=path_to_data + '/samples/{sample}/merged/bin_stats/bin_stats_fit_nbglm.tsv',
        model=path_to_data + '/samples/{sample}/merged/bin_stats/bin_stats_model_nbglm.Rds',
    resources: cpus=1, mem_mb=16000, time_min='5-00:00:00'
    shell:
        'Rscript src/R/cfmedip_nbglm.R -i {input} -o {output.fit} --modelout {output.model}'

rule cfmedip_array_positions:
    input:
        fit=path_to_data + '/samples/{sample}/merged/bin_stats/bin_stats_fit_nbglm.tsv',
        manifest='data/methylation_array_manifest/{array}.hg38.manifest.tsv'
    output:
        path_to_data + '/samples/{sample}/merged/bin_stats/bin_stats_{array}_positions.tsv',
    shell:
        'Rscript ~/git/cfMeDIP-seq-analysis-pipeline/src/R/cfmedip_to_array_sites.R -c {input.fit} -m {input.manifest} -o {output}'


