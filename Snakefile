## April Lo
## Berger Lab, FHCRC
## Seattle, WA, USA
## 2019-08-14
## Updated 2021-05-21
## snakemake 5.2.4

import os
import sys
import re
import pandas as pd

configfile: "config.yaml"

dataName = config["dataName"]            
bamDir = os.path.abspath(config["bamDir"]) + "/"
covDir = os.path.abspath(config["covDir"]) + "/"
metricsDir = os.path.abspath(config["metricsDir"]) + "/"
projectDir = os.getcwd()
rMATSdir = os.path.abspath(config["rMATSdir"]) + "/"

sampleToFastq = {}

## Read in fastq filenames and parse for attributes
fastqfofn = open(config["fastqfofn"], "r")
for fastq in fastqfofn:
    fastq = fastq.strip()
    assert os.path.exists(fastq)
    fastqfile = os.path.basename(fastq)
    construct, rep, sample, read, n = fastqfile.strip().split("_")
    key = construct + "_" + rep
    if (key not in sampleToFastq):
        sampleToFastq[key] = []
    sampleToFastq[key].append(fastq)
fastqfofn.close()

## Generate sampleInfo for RNASeQC and downstream analysis
sampleInfo = open("sampleInfo.txt", "w")
sampleInfo.write("SampleId\tBam File\tNotes\n")
for sample in sampleToFastq.keys():
    group = sample.split("_")[0]
    sampleInfo.write(sample + "\t" + sample + ".bam\t" + group + "\n")

geneToVariant = {}
geneVariantPairs = []

## An input function helper to separate control and test samples
# SetsCoord: these are bams sorted by coord
controlSetsCoord = {}
perturbSetsCoord = {}
for sample in sampleToFastq.keys():
    construct, rep = sample.strip().split("_")
    key = construct
    bam = bamDir + sample + ".sortedByCoord.bam"
    if (construct in ["Renilla"]):
        if (key not in controlSetsCoord):
            controlSetsCoord[key] = []
        controlSetsCoord[key].append(bam)
    else:
        gene, variant = construct.split("-")
        if (gene not in geneToVariant):
            geneToVariant[gene] = []
        if (variant not in geneToVariant[gene] and variant != "WT"):
            geneToVariant[gene].append(variant)
        if (key not in perturbSetsCoord):
            perturbSetsCoord[key] = []
        perturbSetsCoord[key].append(bam)

totalSetsCoord = {**controlSetsCoord, **perturbSetsCoord}

perturbvscontroldirs = []
for perturb in perturbSetsCoord.keys():
    for control in controlSetsCoord.keys():
        resDir = str(perturb) + "_vs_" + str(control)
        perturbvscontroldirs.append(resDir)

variantvsWTdirs = []
for gene in geneToVariant.keys():
    for variant in geneToVariant[gene]:
        resDir = str(gene) + "/" + str(variant) + "_vs_WT"
        variantvsWTdirs.append(resDir)
        
rule all:
    input:
        expand(bamDir + "{sample}.bam",
               sample = sampleToFastq.keys()),
        expand(projectDir + "/QC/readcounts/{sample}.txt",
               sample = sampleToFastq.keys()),
        expand(projectDir + "/QC/rRNAcounts/{sample}.in.bam",
               sample = sampleToFastq.keys()),
        expand(projectDir + "/QC/gene-body-coverage/{sample}.geneBodyCoverage.txt",
               sample = sampleToFastq.keys()),
        expand(metricsDir + "{sample}.insert_size_metrics.txt",
               sample = sampleToFastq.keys()),        
        projectDir + "/QC/gene-body-coverage_all.pdf",        
        "featureCounts/hg19/counts_all.txt",
        spliceresults = expand(rMATSdir + "{perturbvscontroldir}/{AS_type}.MATS.{JC}.txt",
                               perturbvscontroldir = perturbvscontroldirs,
                               AS_type = ["A3SS", "A5SS", "MXE", "RI", "SE"],
                               JC = ["JC", "JCEC"]),
        spliceresultsWT = expand(rMATSdir + "vsWT/{variantvsWTdir}/{AS_type}.MATS.{JC}.txt",
                                 variantvsWTdir = variantvsWTdirs,
                                 AS_type = ["A3SS", "A5SS", "MXE", "RI", "SE"],
                                 JC = ["JC", "JCEC"])
        
## Input function to get fastqs of read group 1 for STAR alignment
def sample_wildcard_to_fastqR1s(wildcards):
    key = str(wildcards.sample)
    fastqs = sampleToFastq[key]
    fastqR1s = []
    for fastq in fastqs:
        readgroup = fastq.strip().split("_")[-2]
        if readgroup == "R1":
            fastqR1s.append(fastq)
    return(ancient(fastqR1s))

## Input function to get fastqs of read group 2 for STAR alignment
def sample_wildcard_to_fastqR2s(wildcards):
    key = str(wildcards.sample)
    fastqs = sampleToFastq[key]
    fastqR2s = []
    for fastq in fastqs:
        readgroup = fastq.strip().split("_")[-2]
        if readgroup == "R2":
            fastqR2s.append(fastq)
    return(ancient(fastqR2s))

## Align RNA-seq paired reads using STAR
## Note this is the most memory intensive step of the pipeline
rule STARalign:
    input:
        fastqR1s = sample_wildcard_to_fastqR1s,
        fastqR2s = sample_wildcard_to_fastqR2s
    output:
        bam = bamDir + "{sample}.sortedByCoord.bam", 
    resources:
        mem = 64,
    params:
        outDir = bamDir + "{sample}/",
        STAR = config["STAR"],
        genomeDir = config["STARgenome"]
    log: "logs/STAR/{sample}.log"
    threads: 8
    shell:
        """
        mkdir -p {params.outDir}
        R1list=$(echo {input.fastqR1s} | tr ' ' ,)
        R2list=$(echo {input.fastqR2s} | tr ' ' ,)
        echo $R1list
        echo $R2list
        {params.STAR} --runThreadN {threads} --genomeDir {params.genomeDir} --readFilesIn "$R1list" "$R2list" --readFilesCommand zcat --outFileNamePrefix {params.outDir} --outSAMtype BAM SortedByCoordinate --twopassMode Basic 2>{log}
        mv {params.outDir}Aligned.sortedByCoord.out.bam {output.bam}
        """

## Picard mark duplicate reads
rule markdups:
    input:
        bam = bamDir + "{sample}.sortedByCoord.bam"
    output:
        bam = temp(bamDir + "{sample, [^.]+}.markdup.bam")
    resources:
        mem = 64,
    params:
        samplename = "{sample}"
    log:
        "logs/picard/{sample}.log"
    threads: 8
    shell:
        """
        java -jar /home/solexa/apps/picard/picard-tools-1.114/MarkDuplicates.jar ASSUME_SORTED=TRUE SORTING_COLLECTION_SIZE_RATIO=0.1 I={input.bam} O={output.bam} M=logs/picard/{params.samplename}.marked_dup_metrics.txt
        """


## Add read group info to bam file headers, using picard
rule addreadgroups:
    input:
        bam = bamDir + "{sample}.markdup.bam"
    output:
        bam = temp(bamDir + "{sample, [^.]+}.info.bam")
    resources:
        mem = 28,
    params:
        samplename = "{sample}",
        flowCell = config["flowCell"]
    log:
        "logs/picard/{sample}.log"
    threads: 4
    shell:
        """
        java -jar /home/solexa/apps/picard/picard-tools-1.114/AddOrReplaceReadGroups.jar I={input.bam} O={output.bam} LB={params.samplename} PU={params.flowCell} PL=Illumina SM={params.samplename} 2>log
        """

## Reorder bam file to match reference using (picard)
rule reorderbam:
    input:
        bam = bamDir + "{sample}.info.bam",
        fasta = "/home/alo2/bergerlab_shared/Projects/RNA_eVIP/reference/hg19/Sequence/genome.fa"
    output:
        bam = bamDir + "{sample, [^.]+}.bam"
    resources:
        mem = 14,
    params:
        tmpDir = bamDir + "{sample}.tmp"
    log:
        "logs/picard/{sample}.log"
    threads: 2
    shell:"""
java -jar /home/solexa/apps/picard/picard-tools-1.114/ReorderSam.jar I={input.bam} O={output.bam} R={input.fasta} TMP_DIR={params.tmpDir} 2>log
rm {input.bam}
"""

## Use picard tool to determine distribution (mean and standard deviation) of insert sizes of paired end reads
rule getInsertSize:
    input:
        bam = bamDir + "{sample}.bam"
    output:
        metrics = metricsDir + "{sample}.insert_size_metrics.txt",
        histogram = metricsDir + "{sample}.insert_size_histogram.pdf"
    log:
        "logs/picard/{sample}.log"
    threads: 1
    shell:"""
java -jar /home/solexa/apps/picard/picard-tools-1.114/CollectInsertSizeMetrics.jar \I={input.bam} \
O={output.metrics} \
H={output.histogram}
"""
    
## Generate index for bam file, using samtools
rule indexbam:
    input:
        bam = bamDir + "{sample}.bam"
    output:
        bai = bamDir + "{sample, [^.]+}.bam.bai"
    resources:
        mem = 14,
    log:
        "logs/samtools/{sample}.log"
    threads: 2
    shell:"""
samtools index {input.bam} 2>log
"""

## RSeQC function for getting read stats from bam file        
rule getReadCounts:
    input:
        bam = bamDir + "{sample}.bam"
    output:
        countfile = projectDir + "/QC/readcounts/{sample}.txt"
    threads: 1
    log: "logs/RSeQC/{sample}_bam_stat.log"
    shell:"""
bam_stat .py -i {input.bam} > {output.countfile} 2>{log}
"""

## RSeQC function for getting the gene body coverages across the experiment
rule getGeneBodyCoverage:
    input:
        bams = expand(bamDir + "{sample}.bam", sample=sampleToFastq.keys()),
        bed = "/home/alo2/bergerlab_shared/Projects/RNA_eVIP/reference/hg19/Annotation/hg19_UCSC_knownGene.bed"
    output:
        projectDir + "/QC/gene-body-coverage_all.pdf"
    params:
        outprefix = projectDir + "/QC/gene-body-coverage_all"
    threads:1
    log: "logs/RSeQC/geneBody_coverage.log"
    shell:"""
inputbams=$(echo {input.bams} | tr ' ' ,)
echo $inputbams
geneBody_coverage.py -i $inputbams -r {input.bed} -f pdf -o {params.outprefix}
"""

## RSeQC function for getting the gene body coverages for each sample
rule getGeneBodyCoverageIndividual:
    input:
        bam = bamDir + "{sample}.bam",
        bed = "/home/alo2/bergerlab_shared/Projects/RNA_eVIP/reference/hg19/Annotation/hg19_UCSC_knownGene.bed"
    output:
        projectDir + "/QC/gene-body-coverage/{sample}.geneBodyCoverage.txt"
    params:
        outprefix = projectDir + "/QC/gene-body-coverage/{sample}"
    threads:2
    log: "logs/RSeQC/{sample}_gene-body-coverage.log"
    conda: "RNAseqQC_env.yml"
    shell:"""
geneBody_coverage.py -i {input.bam} -r {input.bed} -f pdf -o {params.outprefix}
"""

## RSeQC function to count ribosomal RNAs (a useful QC measure)
rule getRibosomalRNAcounts:
    input:
        bam = bamDir + "{sample}.bam",
        bed = "/home/alo2/bergerlab_shared/Projects/RNA_eVIP/reference/hg19/Annotation/hg19_rRNA.bed"
    output:
        rRNAfile = projectDir + "/QC/rRNAcounts/{sample}.in.bam",
        otherfile = projectDir + "/QC/rRNAcounts/{sample}.ex.bam",
        junkfile = projectDir + "/QC/rRNAcounts/{sample}.junk.bam",
        counts = projectDir + "/QC/rRNAcounts/{sample}.counts.txt"
    threads:2
    params:
        outprefix = projectDir + "/QC/rRNAcounts/{sample}"
    log: "logs/RSeQC/{sample}_rRNAcounts.log"
    conda: "RNAseqQC_env.yml"         
    shell:"""
split_bam.py -i {input.bam} -r {input.bed} -o {params.outprefix} > {output.counts}
"""
        
## Run featurecounts from Subread's package to count reads mapped to genes/features
rule featureCounts:
    input:
        bams = expand(bamDir + "{sample}.bam", sample=sampleToFastq.keys()),
        gtf = "/home/alo2/bergerlab_shared/Projects/RNA_eVIP/reference/hg19/Annotation/genes.gtf"
    output:
        counts = "featureCounts/hg19/counts_all.txt"
    resources:
        mem = 4
    params:
        subreadDir = "/home/alo2/apps/subread/subread-1.5.3-Linux-x86_64",
        outDir = "featureCounts"
    threads: 1
    log:
        "logs/featureCounts/" + dataName + ".log"
    shell:
        "{params.subreadDir}/bin/featureCounts -p -t exon -g gene_id -a {input.gtf} -o {output.counts} {input.bams} 2>log"

## Generate coverage pileup files for downstream use
rule mpileup:
    input:
        bams = expand(bamDir + "{sample}.bam", sample=sampleToFastq.keys())
    output:
        pileups = expand(covDir + "{sample}.bam", sample=sampleToFastq.keys())
        
def totalBamsCoord(wildcards):
	key = str(wildcards.perturb)
	return(totalSetsCoord[key])
    
def perturbBamsCoord(wildcards):
	key = str(wildcards.perturb)
	return(perturbSetsCoord[key])
         
def controlBamsCoord(wildcards):
	key = str(wildcards.control)
	return(controlSetsCoord[key])


def repgroup_wildcard_to_bam(wildcards):
    return

## Write pertubation names to file for later use
f = open("perturbations.txt", "w")
for perturb in perturbSetsCoord.keys():
    f.write(perturb)
    f.write("\n")
f.close()

## Run rMATS splicing analysis comparing all genetic perturbations to control vector samples
rule rMATSvsControls:
    input:
        totalBams = totalBamsCoord,
        controlBams = controlBamsCoord
    output:
        expand(rMATSdir + "{{perturb}}_vs_{{control}}/{AS_type}.MATS.{JC}.txt",
               AS_type = ["A3SS", "A5SS", "MXE", "RI", "SE"],
               JC = ["JC", "JCEC"])
    params:
        rMATS = config["rMATS"],
        outDir = rMATSdir + "/{perturb}",
        gtf = config["gtf"],
        readLength = config["readLength"],
        cmpDir = rMATSdir + "{perturb}_vs_{control}"
    resources:
        mem = 16
    log: "logs/rMATS/{perturb}_vs_{control}.log"
    threads: 4
    shell: """
mkdir -p {params.cmpDir}
echo {input.totalBams} | tr ' ' , > {params.cmpDir}/b1.txt
echo {input.controlBams} | tr ' ' , > {params.cmpDir}/b2.txt
python {params.rMATS} --b1 {params.cmpDir}/b1.txt --b2 {params.cmpDir}/b2.txt --gtf {params.gtf} --od {params.cmpDir} -t paired --readLength {params.readLength} --libType fr-unstranded --cstat 0.0001 --nthread {threads} --tstat {threads} &>{log}
"""

def wtBamsCoord(wildcards):
    gene = str(wildcards.gene)
    key = gene + "-WT"
    return(perturbSetsCoord[key])

def variantBamsCoord(wildcards):
    gene = str(wildcards.gene)
    variant = str(wildcards.variant)
    key = gene + "-" + variant
    return(perturbSetsCoord[key])

## Prep files needed for running rMATS
rule preprMATSvsWT:
    input:
        wtBams = wtBamsCoord,
        variantBams = variantBamsCoord
    output:
        b1 = rMATSdir + "vsWT/{gene}/{variant}_vs_WT/b1.txt",
        b2 = rMATSdir + "vsWT/{gene}/{variant}_vs_WT/b2.txt"
    threads: 1
    params:
        outDir = rMATSdir + "vsWT/{gene}/{variant}_vs_WT/"
    resources:
        mem = 4
    log: "logs/rMATS/vsWT/{gene}/{variant}.log"
    run:
        subprocess.run("mkdir -p " + params.outDir, shell=True)        
        # Generate b1 and b2 txt files for rMATS input
        # 1 = variant
        # 2 = WT
        b1 = open(output.b1, "w")
        b1.write(",".join(input.variantBams))
        b1.close()
        b2 = open(output.b2, "w")
        b2.write(",".join(input.wtBams))
        b2.close()
        
## Run rMATS splicing analysis comparing variant perturbations to wild-type perturbations                
rule rMATSvsWT:
    input:
        wtBams = wtBamsCoord,
        variantBams = variantBamsCoord,
        b1 = rMATSdir + "vsWT/{gene}/{variant}_vs_WT/b1.txt",
        b2 = rMATSdir + "vsWT/{gene}/{variant}_vs_WT/b2.txt"        
    output:
        expand(rMATSdir + "vsWT/{{gene}}/{{variant}}_vs_WT/{AS_type}.MATS.{JC}.txt",
               AS_type = ["A3SS", "A5SS", "MXE", "RI", "SE"],
               JC = ["JC", "JCEC"])
    params:
        rMATS = config["rMATS"],
        gtf = config["gtf"],
        readLength = config["readLength"],
        outDir = rMATSdir + "vsWT/{gene}/{variant}_vs_WT",
        cstat = 0.0001
    resources:
        mem = 16
    log: "logs/rMATS/vsWT/{gene}/{variant}.log"
    threads: 4
    shell: """
python {params.rMATS} --b1 {input.b1} --b2 {input.b2} --gtf {params.gtf} --od {params.outDir} -t paired --readLength {params.readLength} --libType fr-unstranded --cstat {params.cstat} --nthread {threads} --tstat {threads} &>{log}
"""        


