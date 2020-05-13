# Snakemake Crash Course

<!-- @import "[TOC]" {cmd="toc" depthFrom=1 depthTo=6 orderedList=false} -->
<!-- code_chunk_output -->

- [Snakemake Crash Course](#snakemake-crash-course)
  - [Prerequisites](#prerequisites)
  - [Introduction to Snakemake](#introduction-to-snakemake)
  - [The Snakefile](#the-snakefile)
    - [Directives:](#directives)
    - [Wildcards and connecting rules](#wildcards-and-connecting-rules)
    - [Rule order and the Target rule](#rule-order-and-the-target-rule)
    - [Cluster jobs, groups, and localrules](#cluster-jobs-groups-and-localrules)
    - [Using `conda` environments](#using-conda-environments)
    - [Comments in the code](#comments-in-the-code)
  - [The configure files](#the-configure-files)
    - [_Standard_ configure files in yaml or json format](#_standard_-configure-files-in-yaml-or-json-format)
    - [_Tabular_ configure files in spreadsheet text format](#_tabular_-configure-files-in-spreadsheet-text-format)
  - [Running Snakemake](#running-snakemake)
    - [Running on a Cluster](#running-on-a-cluster)
  - [Very useful options](#very-useful-options)
    - [Wrapper script for Snakemake](#wrapper-script-for-snakemake)

<!-- /code_chunk_output -->

## Prerequisites
Please complete the [Git Crash Course](/GitCrashCourse/README.md), the [Conda crash course](TBA), and the
[Tmux Crash Course](/TmuxCrashCourse/README.md) before continuing.

## Introduction to Snakemake
`snakemake` is a tool for designing and performing a workflows i.e.,
a sequence of individual steps, each running programs or scripts, that is
required to convert some input, say, sequencing reads, to a desired output,
say, read counts for individual genomic regions of interest.

In one sense, a snakemake workflow can be viewed as just a glorified bash
script performing the desired steps in sequence, but it provides a more
structured layout and is more flexible and easier to rerun/reuse.

In another sense, snakemake could be described as "_python_ with inspiration
from the _make_ program". However, the basic snakemake rule syntax is not so
very much python... and you definitely don't need to be a python expert to
start using snakemake. However, you can use python code snippets _almost_
anywhere in the Snakefile (see examples in [Mapping example Snakefile](/MappingMappingsml)).

## The Snakefile
The cornerstone of snakemake is the Snakefile, which defines the workflow
steps. Each step is coded as a (usually) named `rule` and comprise a
number of `directives`, described below. Additionally, commands defining
configure-files, import of python modules or definition of python
functions can be included in the Snakefile.

### Directives:
Directives The most important directives are the `output`, `input` and `shell`
directives and are defined as follows (although here heavily documented):

```
# Notice the colon (:) after the rule and the directives
# (python style definitions)
rule linkFastq:
    # the output directive is a (possibly nested) list of desired
    # output files, the list items can be named. The output directive
    # defines wildcards inside curly brackets {}; the value of the
    # wildcards is determined when the rule is called either from
    # command-line or by another rule (more about this below).
    output:
        fastq = "fastq/{sample}.fastq"
    # the input directive is a (possibly nested) list of required
    # input files, the list items can be named. It often contains
    # the wildcards defined in the output directive.       
    input:
        fastq = "path/to/originals/{sample}.fastq"
    # The shell minimally contains the bash commands needed to
    # convert input to output.  The code should be inside quotes.
    shell:
        """
        ln -s {input.fastq} {output.fastq}
        """
```

Other important directives include

- `params:` which can define other important parameters, typically
  _not_ files.
- `log:` which defines a log file that can be used to capture
  `std output` and `std error`.
- `conda:` tells snakemake how create a conda environment providing
  programs required by the rule (see examples in MappingSnakefile).
- `group:` assigns a "cluster job group" to a rule, see below under [].

Examples on how these are used can be found in the MappingSnakefile.

### Wildcards and connecting rules
If the a Snakefile with the rule above is run as

`snakemake fastq/mysample.fastq`

snakemake will first check that the file does not already exists,
in which case it reports this and stops. Otherwise, it will parse
the Snakefile to find a rule with output that
matches `fastq/mysample.fastq`. It will see that the rule
`linkFasta` matches if wildcard `"{sample}"` is defined to
be  `"mysample"`, so it defines the wildcard and run the rule.

It will then look for its input file, which after expanding the
wildcard becomes `"path/to/originals/mysample.fastq"`; if this
file exists, snakemake will use this file as the input, Otherwise, it
will look for another rule with an output matching
`"path/to/originals/mysample.fastq"`and run this other rule
to create the input. If you study the workflow in MappingSnakefile,
you will see that different rules have matching output and input;
this creates a sequence of rules to be run -- a *workflow*.

### Rule order and the Target rule
If snakemake is run without an explicit file argument (see below
under [Running Snakemake](#Running\ Snakemake)), the first rule in the Snakefile is by default run. Therefore, the first
rule is often designed as a, so-called, _target_ rule that only
has an input directive comprising the desired final output files
(see example in MappingSnakefile).

With the exception of the Target rule, rules can be written in any
order in the snakefile. It makes sense to have them in 'chronological'
order, i.e., the order they are expected to be executed.

The workflow in [Mapping.smk](Mapping.smk) is set up to ananlyze a single sample
by first soft-link the external files (fastq-file and reference genome fasta
file) needed into the analysis folder, map all reads in the fastq file the
reference genome, and finally sorting and indexing the resulting bam-file. We
can run the full workflow for a single sample, _s1_, by setting the final
indexed bam-file  (`mapped/s1.bam.bai`) as the  argument of `snakemake`.
However, if we list the output files from the final rule (i.e., the indexed
bam-files) for all samples in the Target rule, just running `snakemake` without
argument will go through the full workflow for all samples.


### Cluster jobs, groups, and localrules
When run on a cluster (see below and in `doMapping.sh`), by default,
each rule will be run as a separate job. However, if a few sequential
rules each takes less than, say, 15 min and all uses the same
amount of resources (memory or cpus/threads), then these could
be sent as one single job by assigning them the same *group name*
using the directive `group`, see example below. Moreover, if some jobs
are ridiculously short (e.g., softlinking external files), these
can be run on the login node by listing them as `localrules` (see
example in [Mapping.smk](Mapping.smk)).

### Using `conda` environments
A very nice feature of `snakemake` is that it automatizes the use
of software environments, either using `conda`  or as container
(not treated here), on a rule-based basis. To obtain this

1. first create a conda envirnment yaml-file, including the software,
required by a rule, as dependencies and the corresponding channels
(see [Conda crash course](TBA)).
2. Add the `conda` directive to the rule definition, e.g.,  
```
rule index:
    output:
        bai = "mapped/{prefix}.bam.bai"
    input:
        bam = "mapped/{prefix}.bam",
    params:
        tmp = "mapped/{prefix}.bam.bkp"
    conda: "envs/samtools.yaml"
    log:
        log = "log/index_{prefix}.log"
    shell:
        """
        samtools sort {params.tmp} -o {input.bam}
        samtools index {input.bam}
        """
```
When running this rule with the option `--use-conda` snakemake
will first create the conda environment (if not already existing)
and then activate it before executing `shell` commands of the rule.

### Comments in the code
Comments can be added in either of the following ways:

- `#...` Everything after a hash on a line is ignored by snakemake
- `"""..."""` Everything between triplicated double quotes will be ignored
by snakemake. However, snakemake will still try to expand strings
between curly brackets as variables, which can lead to errors.
- `'''...'''` Everything between triplicated single quotes will be
ignored by snakemake. This will prevent variable expansions,
so use this to comment out code, etc.



## The configure files
It is convenient to define variables, including external input data,
in a configure file. There are two variants of config files.

### _Standard_ configure files in yaml or json format
The _standard_ configure files can be in yaml-format or json format.
An advantage of the yaml format is that it is quite similar to the
directive syntax in snakemake; however, a really annoying drawback
is that indentations _must be as spaces and not as ASCII tab-characters_.
Many, but not all, text editors handles this correctly, however,
i.e., when you hit the tab-key it will insert a sequence of spaces
instead of the ASCII tab character (I use emacs with a yaml mode;
  atom is quite nice and at least have a function for converting
  tabs to spaces).
The json format is more strongly typed; it uses brackets and commas
instead of indentations and requires all strings to be quoted;
also, it does not support comments.

- yaml example

```yaml
# to assign a single value to variable
variable1: 17
# to assign a list (indentations are important)
variable2:
  - value1
  - value2
# to assigna dictionary (indentations are important)
variable3:
  key1 : value1
  key2 : value2
```

- json example

```json
{
  "variable1": 17,
  "variable2": [
    "value1",
    "value2"
  ],
  "variable3": {
  "key1" : "value1",
  "key2" : "value2"
  }
}
```

If the configure file location is defined in the Snakefile, as, e.g.,

```
configfile: "MappingConfig.yaml"
```

either of the configure file examples above will be expand to the same python
dictionary, which can be accessed by `config["variablename"]`. Thus,
`config["variable1"]` will expand to 17, `config["variable2"]` expands to the
python list `["value1", "value2"]`, and `config["variable3"]["key1"]` expand
to `"value1"` (see further _MappingConfig.yaml_").

### _Tabular_ configure files in spreadsheet text format

If we have a lot of samples, for which we each want to include a number
of configuration variables, the standard configure file can become rather
messy and cumbersome. In such a case, it may be better to use a
_tabular_ configure file. This is essentially a spreadsheet, but it must be
in a text file format (not MS Excel, MacOSX Numbers, etc.). This text file
format could be tab-delimited (`.tsv`) or comma-demimited (`.csv`); `snakemake`
uses pythons pandas module which determines the delimiter automatically.
(*__NB!__ Unfortunately, this does not, currently, seem to work seamlessly
  in pandas. Below, I have, instead, manually set the separator to be
  either tab or comma. Unfortunately, this precludes using, e.g.,
  commas inside fields.*)

A tabular configure file could look as follows.

```
# These are my samples
sample  fastq                   readlength  Treatment
s1      /crex/path/to/s1.fastq  71          case
s2      /crex/path/to/s2.fastq  70          control
...
```
This is typicallly read into the Snakefile as follows:

```
import pandas as pd
samples = pd.read_csv("tabularSamples.tsv", sep = "[\t,]", comment="#").set_index("sample")
```
each column can then be accessed in the Snakefile as a dictionary indexed
by the sample name, e.g.,:

`samples["fastq"]["sq"]` expands to `"/crex/path/to/s1.fastq"` and
`samples["reads"]` expands to `"[70,71]"`.

## Running Snakemake
Simply typing `snakemake` on the command line (given that snakemake
is installed on the system), will start snakemake and look for a
Snakefile in the present directory. However, it is good practice to
keep the code and the results in separate folders, so we usually want to
tell snakemake where to find the Snakefile. This is done by the
option `-s`:
```
snakemake -s path/to/Snakefile
```

As mentioned above, when run like this, i.e., without a file argument,
snakemake will execute the first rule in the Snakefile. If this is
designed as a _target file_ (see above), typically the whole workflow
will be executed.

However, often one only wants to perform a single rule (or all rules
up to that rule) for a single sample. This can be achieved by running
snakemake with the output of the rule and sample as a file argument:
```
snakemake -s path/to/Mappingfile path/to/requested/outputfile
```

### Running on a Cluster
__Important__: Snakemake should always be started in a _virtual terminal_
(see further the TmuxCrashCourse.md) on the _login node_. Snakemake will
then, based on the Snakefile, determine jobs and send them to the cluster.
You do not need to manually submit jobs!

To run snakemake in parallel on UPPMAX, you use the following options:
- `--cluster` this defines the sbatch command needed to commit jobs to the
cluster, and uses rule-specific variables defined in a cluster-config file.
- `--cluster-config` this tells snakemake where to look for the
cluster-config file (see further MappingCluster.yaml)

Details on how to run snakemake can be found in the wrapper script `doMapping.sh`

## Very useful options
There are numerous options to snakemake, but I will only list some of
the most useful ones here:

- `-n` dry-run; _always_ run a dry-run before the real run to avoid
unpleasant surprises!
- `-p` print shell commands that will be performed -- good for pre-checking
and debugging rules
- `-r` state reason for rerunning a rule -- again good for pre-checking
and debugging
- `-f` force a rerun of the (last) rule that that has the requested
file as output -- sometimes needs doing!

### Wrapper script for Snakemake
If using a lot of options, the command line call to snakemake can be
quite long and cumbersome, especially when running on a cluster.
Therefore, I usually create a wrapper script that hides all
the absolutely necessary options. I design the wrapper script so that
I can pass other options (e.g., `-n`) to the wrapper script and it
will be passed on to snakemake. See further the script doMapping.sh.
