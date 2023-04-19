base_dir = '/specify/working/directory'

import pandas as pd
p1 = pd.read_csv('data/20211109_sample_sheet_pool_1.csv')
p2 = pd.read_csv('data/20211109_sample_sheet_pool_2.csv')
p3 = pd.read_csv('data/20211109_sample_sheet_pool_3.csv')
#Same sheets as the ones used for demultiplexing. Should have 3 columns: Sample_ID, Sample_Name and Sample_Barcode
ms = pd.concat([p1,p2,p3])

#IMPORTANT, select here samples to process, or all samples
#S = ['M13', 'I13','K13']
#S = ['M13']
S = ms.Sample_Name.tolist()

def getR1(S):
    return [f'{base_dir}demux_pool_{ms.loc[ms["Sample_Name"]==s,"Library_ID"].item()}/{ms.loc[ms["Sample_Name"]==s,"Sample_ID"].item()}-{s}-{ms.loc[ms["Sample_Name"]==s,"Sample_Barcode"].item()}_R1.fastq.gz' for s in S]
#For each sample identified by Sample_Name, function returns full path to file
#Files correspond to demultiplexing output with no trimming

def getR2(S):
    return [f'{base_dir}demux_pool_{ms.loc[ms["Sample_Name"]==s,"Library_ID"].item()}/{ms.loc[ms["Sample_Name"]==s,"Sample_ID"].item()}-{s}-{ms.loc[ms["Sample_Name"]==s,"Sample_Barcode"].item()}_R2.fastq.gz' for s in S]
#Same for corresponding Read2 file

rule all:
    input:
        expand(base_dir+'aligned/{sample}.sorted.{ext}', sample=S, ext=['depth.out', 'stat', 'calmd.bam']),
        expand(base_dir+'filtered/{sample}.sorted.calmd.bam', sample=S)
#Should be changed whenever a new rule is added, corresponds to final target
#depth.out files can be used to infer copy number variations
#depth.out and all bam files are very heavy, prefer temp files as output, especially when running the final script on all samples

rule bwa_mem2:
    input:
        read1 = getR1,
        read2 = getR2
    params:
        genome = base_dir+'genome/S288C_chr'
#Ref here is path + prefix used for indexing, full name of reference is indicated below just in case
#_reference_sequence_R64-3-1_20210421_chr.fsa
#Command to index: bwa-mem2 index -p S288C /specify/working/directory/genome/S288C_reference_sequence_R64-3-1_20210421_chr.fsa
    output:
        base_dir+'aligned/{sample}.sorted.bam'
        #temp(base_dir+'aligned/{sample}.sorted.bam')
        #Directory and content are deleted when all rules which use it as input have used it
        #Allows to gain some space by deleting heavy bam files
    shell:
        'bwa-mem2 mem {params.genome} {input.read1} {input.read2} | samtools sort -T tpm.{wildcards.sample} -o {output}; '
        'samtools index {output}'

rule bam_stats:
    input:
        rules.bwa_mem2.output
    output:
        base_dir+'aligned/{sample}.sorted.stat'
    shell:
        'samtools stats {input} > {output}'

rule depth:
    input:
        rules.bam_stats.input
    output:
        base_dir+'aligned/{sample}.sorted.depth.out'
    shell:
        'samtools depth -a {input} > {output}'

#Remove reads that do not map uniquely
rule unique_align:
    input:
        rules.bwa_mem2.output
    output:
        base_dir+'filtered/{sample}.filtered.bam'
    shell:
        'samtools view -F 256 {input} -o {output}; '
        'samtools index {output}'

#Remove PCR duplicates
rule rmdup:
    input:
        base_dir+'filtered/{sample}.filtered.bam'
    output:
        rmdup = base_dir+'filtered/{sample}.rmdup.bam',
        #rmdup = temp(base_dir+'filtered/{sample}.rmdup.bam'),
        metrics = base_dir+'filtered/{sample}.sorted.rmdup.metrics'
    shell:
        'java -jar /path/to/picard.jar MarkDuplicates -I {input} -O {output.rmdup} -M {output.metrics} -REMOVE_DUPLICATES true'

#Realign reads around indels
rule calmd:
    input:
        base_dir+'filtered/{sample}.rmdup.bam'
    params:
        ref = base_dir+'genome/S288C_reference_sequence_R64-3-1_20210421_chr.fsa'
    output:
        base_dir+'filtered/{sample}.sorted.calmd.bam'
    shell:
        'samtools calmd -bAr -@2 {input} {params.ref} > {output}'