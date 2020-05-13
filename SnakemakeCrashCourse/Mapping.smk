""" Snakefile with rules to run the Mapping workflow."""

# module import
""" Here we import a python module that we need for the following two
   python functions """
import os


# Two python functions that just allow simpler notation in rules
""" It is good keep the code and the results separate, so
    we run the workflow in a separate working directory. However,
    we sometimes need to access files in the code folder. This
    function prepends a file path with the current path to the workflow
    base directory (kept in snakemake variable workflow.basedir) """
def prependWfd(path):
    return(os.path.normpath(os.path.join(workflow.basedir, path)))

""" Mostly, snakemake uses relative paths (to the working directory)
    when referring to, e.g., input and output files. However, when
    soft-linking (creating aliases) from external file, full paths is
    required by the system. This function prepends the path to the current
    working directory (accessed by the python function os.getcwd()). """
def prependPwd(path):
    return(os.path.normpath(os.path.join(os.getcwd(), path)))


# Configure files
""" Read the yaml config file in the present analysis directory. If
    not existing, warn user and copy the one from the workflow source
    directory. """
if not os.path.exists("MappingConfig.yaml"):
    print("No user-defined config file, MappingConfig.yaml, found in the"
          "present directory. Copying the config file from the workflow"
          "source directory. You probably want to edit this and add your"
          "own changes.")
    os.system("cp {} MappingConfig.yaml".format(prependWfd('MappingConfig.yaml')))
configfile: "MappingConfig.yaml"

""" Read sample info from the tabular config file 'samples.tsv' (using
    the python module pandas. This can be a tab-separated or
    comma-separated spreadsheet text file. If not existing, warn user
    and copy the one from the workflow source directory. columns can
    be accessed as 'samples['<colname>']' and cells can be accessed as
    'samples['<col>']['<sample name>']'. """
import pandas as pd
if not os.path.exists("MappingSamples.tsv"):
    print("No user-defined tabular config file, MappingSamples.yaml, \n"
          "found in the present directory. Copying the config file from \n"
          "the workflow source directory. You probably want to edit this \n
          "and add your own changes.")
    os.system("cp {} MappingSamples.tsv".format(prependWfd('MappingConfigSamples.tsv')))
samples = pd.read_csv("MappingSamples.tsv", sep="[\t,]",
                      comment="#", engine="python").set_index("sample")

# Cluster, groups, and localrules
""" The rules that softlinks external files) are ridiculously fast;
    run these on the login node by listing them as 'localrules'. """
localrules: linkFastq, linkReference, linkRegions


# Target rule
""" By listing the final count files for all samples, this will run the
    whole workflow for all samples """
rule all:
    input:
        # expand creates a list by replacing {sample} by all
        # sample names in samples
        count = expand(
            "count/{sample}.bam.bai",
            sample = samples.index
        )

# Rule order
""" The remaining rules are written in the order of expected execution"""

# It is good practice to softlink all external files into the working directory
rule linkFastq:
    output:
        fastq = "fastq/{sample}.FASTQ.gz"
    input:
        # see further MappingConfig.yaml
        fastq = lambda wc: samples["fastq"][wc.sample]
    params:
        # full path needed for ln -s, neater to redefine in params.
        fastq = lambda wc: prependPwd("fastq/{s}.FASTQ.gz".format(s=wc.sample))
    log:
        # must have (all) same wildcards as output
        log = "log/linkFastq_{sample}.log"
    shell:
        """
        # everything inside curly braces is expanded to the corresponding variable defined above,
        # e.g., {log.log} becomes "log/linkFastq_<theactualSampleName>.log"
        exec &> {log.log}

        echo "Linking {input.fastq} to {output.fastq}"
        ln -s {input.fastq} {params.fastq}
        echo "Done!"
        """

rule linkReference:
    output:
        reference = "metadata/reference.fasta"
    input:
        reference = config["reference"]
    params:
        reference = prependPwd("metadata/reference.fasta")
    log:
        log = "log/linkReference.log"
    shell:
        """
        exec &> {log.log}

        echo "Linking {input.reference} to {output.reference}"
        ln -s {input.reference} {params.reference}
        echo "Done!"
        """

rule map:
    output:
        bam = "mapped/{sample}.FASTQ_Mapping_mapped.bam",
    input:
        reference = "metadata/reference.fasta",
        fastq = "fastq/{sample}.FASTQ.gz"
    params:
        outdir = "mapped"
    # The conda directive points to a conda environment file (see
    # documentation in that file). If snakemake is run with option
    # `--use-conda`, it will create the conda environment and run
    # the present rule within the environment.
    conda: "envs/Mapping.yaml"
    # The group directive enables sending several rules as on cluster
    # jobs (deafult is one rule: one job). Rules with the same group
    # name will be run in the same cluster job.
    group: "MapAndIndex"
    log:
        log = "log/map_{sample}.log"
    shell:
        """
        exec &> {log.log}

        echo "Running bowtie version:"
        bwowtie --version
        echo "to map{input.fastq} to {input.reference}"
        zcat {input.reads} | bowtie -m 1000 -a -v 0 -p {threads} -S {params.refprefix} - {output.sam} > {output.logfile}
        ls {params.outdir}
        echo "Done"
        """

rule index:
    output:
        bai = "mapped/{prefix}.bam.bai"
    input:
        bam = "mapped/{prefix}.bam",
    params:
        tmp = "mapped/{prefix}.bam.bkp"
    conda: "envs/samtools.yaml"
    group: "MapAndIndex"
    log:
        log = "log/index_{prefix}.log"
    shell:
        """
        exec &> {log.log}
        echo "Create backup"
        mv {input.bam} {params.tmp}

        echo "Sorting and indexing bam"
        samtools sort {params.tmp} -o {input.bam}
        samtools index {input.bam}

        echo "Removing backup"
        rm -f {params.tmp}
        echo "Done"
        """

rule count:
    output:
        tcount = "count/{sample}.FASTQ_Mapping_mapped_filtered_tcount.tsv"
    input:
        snp = "snp/{sample}.FASTQ_Mapping_mapped_filtered_snp.vcf",
        reference = "metadata/reference.fasta",
        regions = "metadata/regions.bed",
        bam = "filtered/{sample}.FASTQ_Mapping_mapped_filtered.bam",
        bai = "filtered/{sample}.FASTQ_Mapping_mapped_filtered.bam.bai"
    params:
        outdir = "count/",
        snpdir = "snp/",
        readlength = lambda wc: samples["readlength"][wc.sample]
    conda: "envs/Mapping.yaml"
    group: "snpAndCount"
    log:
        log = "log/count_{sample}.log"
    shell:
        """
        exec &> {log.log}

        echo "Running Mapping version:"
        Mapping --version
        echo "to perform dunk count on {input.bam}"
        time Mapping count -o {params.outdir} -s {params.snpdir}  -r {input.reference} -b {input.regions} -t 2 {input.bam} -l {params.readlength}
        # additional possible options:
        # -m
        # -q <minimum base quality>

        echo "Done"
        """
