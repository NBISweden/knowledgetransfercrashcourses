""" Snakefile with rules to run the Mapping workflow."""

# Comments:
""" This file is extensively commented for course purposes.
    In a real case, comments should be aimed at enhancing reproducibility,
    but should not impair readability. They could include, e.g., brief
    descriptions of a rule's purpose, explaining specific, unintuitive
    usages in code, or to structure large snakefiles into logical blocks. """


# module import
""" Here we import a python module that we need for the following two
   python functions """
import os, re


# Three python functions that just allow simpler notation in rules
""" It is good keep the code and the results separate, so
    we run the workflow in a separate working directory. However,
    we sometimes need to access files in the code folder. This
    function prepends a file path with the current path to the workflow
    base directory (kept in snakemake variable workflow.basedir).
    Functions from the python module os are used to ascertain as
    a correct path."""
def prependWfd(path):
    return(os.path.normpath(os.path.join(workflow.basedir, path)))

""" Mostly, snakemake uses relative paths (to the working directory)
    when referring to, e.g., input and output files. However, when
    soft-linking (creating aliases) from external file, full paths is
    required by the system. This function prepends the path to the current
    working directory (accessed by the python function os.getcwd()). """
def prependPwd(path):
    return(os.path.normpath(os.path.join(os.getcwd(), path)))

""" Special function for the purpose of this crash course.
    *** Do not use in rules! ***
    Updates the dummy paths (beginning with ''/path/to/'in the given configuration file to a correct absolute path corresponding to
    the Example data. NB! In a real case, the data would not be stored
    in a git repo. Instead the known absolute paths to the data would
    be given directly in the config file. """
def correctPathsInConfigFiles(configfile):
    input = open(configfile, 'r')
    ret = []
    for row in input:
        ret.append(re.sub("/path/to/(ExampleData/\S+)", prependWfd("../\g<1>"), row))
    input.close()
    output = open(configfile,'w')
    for row in ret:
        output.write(row)
    output.close()

    # For python geeks: Same function can be written more compactly as:
    # def correctPathsInConfigFiles(configfile):
    #     with open(configfile, 'r') as input: a=[ re.sub("/path/to/(ExampleData/\S+)", prependWfd("../\g<1>"), row) for row in input ]
    #     with open(configfile,'w') as output: [ output.write(i) for i in a ]


# Configure files
""" Read the yaml config file in the present analysis directory. If
    not existing, warn user and copy the one from the workflow source
    directory. """
if not os.path.exists("mappingConfig.yaml"):
    print("No user-defined config file, mappingConfig.yaml, found in the\n"
          "present directory. Copying the config file from the workflow\n"
          "source directory. You probably want to edit this and add your\n"
          "own changes.")
    # copyAndUpdateConfigFiles('mappingConfig.yaml')
    os.system("cp {} mappingConfig.yaml".format(prependWfd('mappingConfig.yaml')))
    # Special path correction for the purpose of this course
    correctPathsInConfigFiles('mappingConfig.yaml')
configfile: "mappingConfig.yaml"

""" Read sample info from the tabular config file 'samples.tsv' (using
    the python module pandas). This can be a tab-separated or
    comma-separated spreadsheet text file. If not existing, warn user
    and copy the one from the workflow source directory. columns can
    be accessed as 'samples['<colname>']' and cells can be accessed as
    'samples['<col>']['<sample name>']'. """
import pandas as pd
if not os.path.exists("mappingSamples.tsv"):
    print("No user-defined tabular config file, mappingSamples.yaml, \n"
          "found in the present directory. Copying the config file from \n"
          "the workflow source directory. You probably want to edit this \n"
          "and add your own changes.")
    os.system("cp {} mappingSamples.tsv".format(prependWfd('mappingSamples.tsv')))
    # Special path correction for the purpose of this course
    correctPathsInConfigFiles('mappingsamples.tsv')
samples = pd.read_csv("mappingSamples.tsv", sep="[\t,]",
                      comment="#", engine="python").set_index("sample")


# Cluster, groups, and localrules
""" The rules that softlinks external files) are ridiculously fast;
    run these on the login node by listing them as 'localrules'. """
localrules: linkFastq, linkReference


# Target rule
""" By listing the final count files for all samples, this will run the
    whole workflow for all samples """
rule all:
    input:
        # expand creates a list by replacing {sample} by all
        # sample names in samples
        mapped = expand(
            "mapped/{sample}.bam.bai",
            sample = samples.index
        )


# Workflow
""" The remaining rules are written in the order of expected execution"""


# Soft-linking
""" It is good practice to softlink all external files into the
    working directory """
rule linkFastq:
    output:
        fastq_1 = "fastq/{sample}_R1.fastq.gz",
        fastq_2 = "fastq/{sample}_R2.fastq.gz"
    input:
        # see further mappingConfig.yaml
        fastq_1 = lambda wc: samples["fastq_1"][wc.sample],
        fastq_2 = lambda wc: samples["fastq_2"][wc.sample]
    params:
        # full path needed for ln -s, neater to redefine in params.
        fastq_1 = lambda wc: prependPwd("fastq/{s}_R1.fastq.gz".format(s=wc.sample)),
        fastq_2 = lambda wc: prependPwd("fastq/{s}_R2.fastq.gz".format(s=wc.sample))
    log:
        # must have (all) same wildcards as output
        log = "log/linkFastq_{sample}.log"
    shell:
        """
        # everything inside curly braces is expanded to the corresponding
        # variable defined above, e.g., {log.log} becomes
        # "log/linkFastq_<theactualSampleName>.log"
        exec &> {log.log}

        echo "Linking {input.fastq_1} to {params.fastq_1}"
        ln -s {input.fastq_1} {params.fastq_1}

        echo "Linking {input.fastq_2} to {output.fastq_2}"
        ln -s {input.fastq_2} {params.fastq_2}

        echo "Done!"
        """

rule linkReference:
    output:
        reference = "metadata/reference.fasta.gz"
    input:
        reference = config["reference"]
    params:
        reference = prependPwd("metadata/reference.fasta.gz")
    log:
        log = "log/linkReference.log"
    shell:
        """
        exec &> {log.log}

        echo "Linking {input.reference} to {output.reference}"
        ln -s {input.reference} {params.reference}
        echo "Done!"
        """

""" temporarily unzip the gzip fasta or fastq when needed by another rule """
rule unzip:
    output:
        # The temp(...) indicate a temporary file required by another rule.
        # It will be removed automatically when that rule is run
        unzipped = temp("{prefix}.{suffix, fast[aq]}")
    input:
        zipped = "{prefix}.{suffix}.gz"
    shell:
        """
        gunzip -c {input.zipped} > {output.unzipped}
        """

""" convert genome fasta file to bowtie format """
rule bowtieIndex:
    output:
        index ="metadata/contigs.1.ebwt",
    input:
        reference = "metadata/reference.fasta"
    params:
        refprefix = "metadata/contigs"
    conda: "envs/bowtie.yaml"
    log: "log/bowtieReference.log"
    shell:
        """
        exec &> {log}
        echo "Constructing bowtie reference"

        bowtie-build -f {input.reference} {params.refprefix}

        echo "Done"
        """

""" perform the mapping using bowtie """
rule bowtie:
    output:
        sam = temp("mapped/{sample}.sam")
    input:
        # input could be a python list of several files (we could have
        # used `expand`, here as in the `all` rule -- how?)
        fastq_1 = "fastq/{sample}_R1.fastq",
        fastq_2 = "fastq/{sample}_R2.fastq",
        reference ="metadata/contigs.1.ebwt"
    params:
        refprefix = "metadata/contigs",
    # The conda directive points to a conda environment file (see
    # documentation in that file). If snakemake is run with option
    # `--use-conda`, it will create the conda environment and run
    # the present rule within the environment.
    conda: "envs/bowtie.yaml"
    # The group directive enables sending several rules as on cluster
    # jobs (deafult is one rule: one job). Rules with the same group
    # name will be run in the same cluster job.
    group: "MapAndIndex"
    log:
        log = "log/map_{sample}.log"
    shell:
        """
        exec &> {log.log}

        echo "Running bowtie to map {input.fastq_1} and {input.fastq_2} to {input.reference}"

        bowtie -a -S {params.refprefix} -1 {input.fastq_1} -2 {input.fastq_2} {output.sam}

        echo "Done"
        """

rule sortAndIndex:
    output:
        bam = "mapped/{prefix}.bam",
        bai = "mapped/{prefix}.bam.bai"
    input:
        sam = "mapped/{prefix}.sam",
    conda: "envs/samtools.yaml"
    group: "MapAndIndex"
    log:
        log = "log/index_{prefix}.log"
    shell:
        """
        exec &> {log.log}
        echo "Create backup"

        echo "Creating sorted and indexed bam"
        samtools sort {input.sam} -o {output.bam}
        samtools index {output.bam}

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
    conda: "envs/featureCounts.yaml"
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
