# Snakemake Crash Course

**TODO**

- Maybe add rule to snakefile that downloads the reference (using wget) to
reduce repo size

## Prerequisites

Please complete the [Git Crash Course](/GitCrashCourse/README.md), the [Conda
Crash Course](TBA), and the [Tmux Crash Course](/TmuxCrashCourse/README.md)
before continuing.

## Introduction

In many projects, the process of going raw data to final analysis results
comprise several sequential analysis steps. For example, in a RNAseq analysis,
this can include:

- adaptor trimming.
- mapping reads to genome,
- estimating read coverage (counts) of featureCounts,
- various QC steps, and
- differential expression analysis

This sequence of steps is usually called a *workflow* or a *pipeline*.
Typically, a majority of the steps are performed in the same way for the
different input samples, which often can constitute a large number. It is
therefore desirable to automatize the workflow so that that type of common
steps can be easily repeated for all samples.

**Snakemake** is a tool for designing and performing automated workflows,
i.e., a sequence of individual steps (each running dedicated programs or
scripts) that is required to convert some input, say, sequencing reads,
to a desired output, say, read counts for individual genomic regions of
interest. The workflow can be branching, coalescing, and parts of it can
be run in parallel, typically on a high-performance computer (HPC) cluster.

If you already shave used dedicated `bash` scripts that run various scripts
on your data in the right order, then, in one sense, a snakemake workflow can
be viewed as just a glorified version of such of`bash` script. However,
snakemake provides a more structured layout and is more flexible and easier
to rerun/reuse.

The structure of the workflow, that is, the different analysis steps to be run,
is defined in [The Snakefile](#The\ Snakefile). In the Snakefile, each step
in the pipeline is performed by, so-called, *rules*. Each rule can be viewed
as corresponding to a good cooking recipe; it states what the *output* should
be (e.g., a meringue), what *input* is required (egg-white and sugar) and what
actions are needed to obtain the output from the input (use a mixer
to whip the egg-whites, add the sugar and bake in oven). In a bioinformatics
pipeline, the output and the input would both be files and the actions would,
e.g., be shell commands (the mixer and the oven could be thought of as
representing other programs or scripts that are called).

By using a consistent system of naming the input and output of the different
steps, these are connected so that the output of one analysis step is the input
of the next analysis step. Moreover, by using *wildcards*  allows the same step
be run on different inputs. The workflow can then be run all steps in on go or
it is possible to run just part of the workflow by specifying the desired output
as arguments when [Running snakemake](#Running\ Snakemake). It can also be
modified to, e.g., different inuts, by the variable values set in
[The configure files](#The\ configure\ files).

Sometimes, snakemake are described as "_python_ with inspiration from the
_make_ program". However, the basic snakemake rule syntax is not so very much
python... and you definitely don't need to be a python expert to start using
snakemake. However, you _can_ use python code snippets _almost_ anywhere in
the Snakefile (see examples in the Mapping example Snakefile
[mapping.smk](./mapping.smk)).

The different parts that constitute a snakemake workflow are described below.

## The Snakefile

The cornerstone of snakemake is the Snakefile, which defines the workflow steps.
Each step is coded as a (usually) named `rule` and comprise a number of
`directives` (e.g., `output`, `input` and `shell`), described below.
Additionally, commands defining configure-files, import of python modules or
definition of python functions can be included in the Snakefile.

### Directives

The most important directives are the `output`, `input` and `shell` directives.
These descibe which output files will be produced, which input files are needed
and what `bash` commands (e.g., calls to programs) are needed to produce the
output from the input. Their use in the Snakefile is described in the code
snippet and the associated texts below. In the code snippet, extensive 
*comments* (starting with a `#`) are used to explain the *following* code
line(s) (e.g., a directive). Notice that comments are ignored (not executed)
by snakemake when running the workflow (see further
[Comments in the code](#Comments\ in\ the\ code)).

```
# This is how a rule is written in the Snakefile. Notice that a
# colon (:) is required after the rule name and after the directives
# (compare python style function definitions)
rule sortAndIndex:

    # The output directive is a (possibly nested) list of desired
    # output files, the list items can be named (notice that commas
    # are required between list items). The output directive defines
    # wildcards inside curly brackets {}; the value of the wildcards
    # is determined when the rule is called either from command-line
    # or by another rule (more about this below). If an output folder
    # (here "mapped" is missing, snakemake creates it automatically.)
    output:
        bam = "mapped/{prefix}.bam",
        bai = "mapped/{prefix}.bam.bai"

    # The input directive is a (possibly nested) list of required
    # input files, the list items can be named (notice that commas
    # are required between list items). It often contains the wildcards
    # defined in the output directive.       
    input:
        sam = "mapped/{prefix}.sam",

    # The shell minimally contains the bash commands needed to
    # convert input to output.  The code should be quoted (use triple
    # quotes for code blocks and single quotes for single lines).
    shell:
        """
        echo "Create backup"

        echo "Creating sorted and indexed bam"
        samtools sort {input.sam} -o {output.bam}
        samtools index {output.bam}

        echo "Done"
        """
```

- The `output` directive describes what files the rule can produce, which, in
the example above, expands to the python dictionary
`{"bam": "mapped/{prefix}.bam", "bai": "mapped/{prefix}.bam.bai"}`.   
- The `input` directive describes the input files the rule require to do it job,
which, in the example above, expands to the (single-itemed) python dictionary
`{"sam": "mapped/{prefix}.sam"}`.  
- The `shell` directive, describes the bash commands used to produce the output
from the input. (Note, the `shell` directive can be replaced by the `script` or
`run` directives, which are not covered here.)

Other important directives include

- `params:` which can define other important parameters, typically _not_ files.

- `log:` which defines a log file that can be used to capture `std output`
and`std error`.

- `conda:` tells snakemake how create a conda environment providing programs
required by the rule.

- `group:` assigns a "cluster job group" to a rule.

Examples on how these are used can be found below and in the Snakefile
[mapping.smk](/mapping.smk).


### Wildcards and connecting rules

If `snakemake` is run, with a Snakefile including the rule above, and requesting
the output `mapped/mysample.bam` , then `snakemake` will first check that the
requested output file does not already exists, in which case it reports this and
stops. Otherwise, it will parse the Snakefile looking for a rule with output
that matches the required output `mapped/mysample.bam`. It will find that the
rule `sortAndIndex` matches if wildcard `"{sample}"` is defined to be
`"mysample"`, so it defines the wildcard `"{sample}" = "mysample"`and prepares
to run the rule.

It will then look for its input file required by the rule, which after expanding
the wildcard becomes `mapped/mysample.sam`; if this file exists, snakemake will
use this file as the input, Otherwise, it will look for another rule with an
output matching `mapped/mysample.sam`and run this other rule to create the
input. If you study the workflow in [mapping.smk](./mapping.smk), you will see
that different rules have matching output and input; this behaviour creates a
sequence of rules to be run to create `mapped/mysample.bam` -- a *workflow*!.
(Can you see which rule in [mapping.smk](./mapping.smk) matches the output
`"mapped/mysample.sam"`? What is its input after expanding the wildcard?).

### Rule order and the Target rule

If snakemake is run without an explicit file argument (see below under [Running
Snakemake](#Running\ Snakemake)), the first rule in the Snakefile is by default
run. Therefore, the first rule is often designed as a, so-called, _target_ rule
that only has an input directive comprising the desired final output files (see
example in [mapping.smk](./mapping.smk)).

With the exception of the Target rule, rules can be written in any order in the
snakefile. It makes sense to have them in 'chronological' order, i.e., the order
they are expected to be executed.

The workflow in [mapping.smk](./mapping.smk) is set up to analyze a
single sample by first soft-link the external files (fastq-file and reference
genome fasta file) needed into the analysis directory, map all reads in the
fastq file the reference genome, and finally sorting and indexing the resulting
bam-file. We can run the full workflow for a single sample, _s1_, by setting the
final indexed bam-file  (`mapped/s1.bam.bai`) as the  argument of `snakemake`.
However, if we list the output files from the final rule (i.e., the indexed
bam-files) for all samples in the Target rule, just running `snakemake` without
argument will go through the full workflow for all samples.


### Cluster jobs, groups, and localrules

When run on a cluster (see below and in `doMapping.sh`), by default, each rule
will be run as a separate job. However, if a few sequential rules each takes
less than, say, 15 min and all uses the same amount of resources (memory or
cpus/threads), then these could be sent as one single job by assigning them the
same *group name* using the directive `group`, see example below. Moreover, if
some jobs are ridiculously short (e.g., soft-linking external files), these can
be run on the login node by listing them as `localrules` (see example in
[mapping.smk](./mapping.smk)).

### Using `conda` environments

A very nice feature of `snakemake` is that it automatizes the use of software
environments, either using `conda`  or as container (not treated here), on a
rule-based basis. To obtain this:

1. first create a conda environment yaml-file, including the software, required
by a rule, as dependencies and the corresponding channels (see [Conda crash
course](TBA)).

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

When running this rule with the option `--use-conda` snakemake will first create
the conda environment (if not already existing) and then activate it before
executing `shell` commands of the rule.

### Comments in the code

Comments can be added in either of the following ways:

- `#...` Everything after a hash on a line is ignored by snakemake

- `"""..."""` Everything between triplicated double quotes will be ignored by
  snakemake. However, snakemake will still try to expand strings between curly
  brackets as variables, which can lead to errors.

- `'''...'''` Everything between triplicated single quotes will be ignored by
  snakemake. This will prevent variable expansions, so use this to comment out
  code, etc.


## The configure files

It is convenient to define variables, including external input data, in a
configure file. There are two variants of config files.

### _Standard_ configure files in yaml or json format

The _standard_ configure files can be in yaml-format or json format. An
advantage of the yaml format is that it is quite similar to the directive syntax
in snakemake; however, a really annoying drawback is that indentations _must be
as spaces and not as ASCII tab-characters_. Many, but not all, text editors
handles this correctly, however, i.e., when you hit the tab-key it will insert a
sequence of spaces instead of the ASCII tab character (I use emacs with a yaml
mode; atom is quite nice and at least have a function for converting tabs to
spaces). The json format is more strongly typed; it uses brackets and commas
instead of indentations and requires all strings to be quoted; also, it does not
support comments.

- yaml example

```yaml
# to assign a single value to variable
variable1: 17
# to assign a list (indentations are important)
variable2:
  - value1
  - value2
# to assign a dictionary (indentations are important)
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
configfile: "mappingConfig.yaml"
```

either of the configure file examples above will be expand to the same python
dictionary, which can be accessed by `config["variablename"]`. Thus,
`config["variable1"]` will expand to 17, `config["variable2"]` expands to the
python list `["value1", "value2"]`, and `config["variable3"]["key1"]` expand to
`"value1"` (see further [mappingConfig.yaml](./mappingConfig.yaml)).

### _Tabular_ configure files in spreadsheet text format

If we have a lot of samples, for which we each want to include a number of
configuration variables, the standard configure file can become rather messy and
cumbersome. In such a case, it may be better to use a _tabular_ configure file.
This is essentially a spreadsheet, but it must be in a text file format (not MS
Excel, MacOSX Numbers, etc.). This text file format could be tab-delimited
(`.tsv`) or comma-demimited (`.csv`); `snakemake` uses pythons pandas module
which determines the delimiter automatically. (*__NB!__ Unfortunately, this does
not, currently, seem to work seamlessly in pandas. Below, I have, instead,
manually set the separator to be either tab or comma. Unfortunately, this
precludes using, e.g., commas inside fields.*)

A tabular configure file could look as follows.

```
# These are my samples
sample  fastq                   readlength  Treatment
s1      /crex/path/to/s1.fastq  71          case
s2      /crex/path/to/s2.fastq  70          control
...
```
This is typicallly read into the Snakefile as follows:

```python
import pandas as pd
samples = pd.read_csv("tabularSamples.tsv", sep = "[\t,]", comment="#").set_index("sample")
```

each column can then be accessed in the Snakefile as a dictionary indexed by the
sample name, e.g.,:

`samples["fastq"]["sq"]` expands to `"/crex/path/to/s1.fastq"` and
`samples["reads"]` expands to `"[70,71]"`.

## Running Snakemake

Simply typing `snakemake` on the command line (given that snakemake is installed
on the system), will start snakemake and look for a Snakefile in the present
directory. However, it is good practice to keep the code and the results in
separate directories, so we usually want to tell snakemake where to find the
Snakefile. This is done by the option `-s`:

```
snakemake -s path/to/Snakefile
```

As mentioned above, when run like this, i.e., without a file argument, snakemake
will execute the first rule in the Snakefile. If this is designed as a _target
file_ (see above), typically the whole workflow will be executed.

However, often one only wants to perform a single rule (or all rules up to that
rule) for a single sample. This can be achieved by running snakemake with the
output of the rule and sample as a file argument:

```
snakemake -s path/to/Snakefile path/to/requested/outputfile
```

### Using `conda` environments

Activating the use of `conda` environments, for rules with `conda` directives,
is done with the option `--use-conda`, e.g.,

```
snakemake -s path/to/Snakefile --use-conda path/to/requested/outputfile
```

### Running on a Cluster

**Important:** Snakemake should always be started in a *virtual terminal*
(see further the TmuxCrashCourse.md) on the *login node*. Snakemake will
then, based on the Snakefile, determine jobs and send them to the cluster.
*You do not need to manually submit jobs!*

To run snakemake in parallel on UPPMAX, the following options are used:

- `--cluster`: this defines the sbatch command needed to commit jobs to the cluster, and uses rule-specific variables defined in a *cluster-config file*.

- `--cluster-config`: this tells snakemake where to look for the *cluster-config file*.


#### The cluster config file

The yaml-formatted cluster-config file defines values for the various `sbatch`
options (see [UPPMAX `slurm` user guide](https://www.uppmax.uu.se/support/user-guides/slurm-user-guide/)), e.g.,

```
__default__ :
  name: default            # name of job shown in, e.g., jobinfo
  account: snicXXX-YYY     # Uppmax account to use
  time: 0-02:00:00         # Allocated time
  partition: core          # Type of computational unit allocated, core or node
  n: 1                     # Allocated number of computational units
  other: "--mail-user your@email --mail-type=FAIL,STAGEOUT" # email are sent to this address on failed jobs

mapAndFilter:
  name: mapAndFilter
  time: 0-01:30:00
```

The `__default__` values are implicitly used for all non-local rules, while
explicit rule- or group-specific settings change one or more individual values
for the corresponding rule or group (see further examples in
[mappingCluster.yaml](mappingCluster.yaml)), These values can be accessed in the
`--cluster` `sbatch` command as wildcards of the form `{cluster.variable}`.
These wildcards are expanded on a rule basis, so that for rules in the group
`mapAndFilter` `{cluster.time}` will expand to 1.5h, while for all other rules
it will expand to 2h.

#### The `sbatch` command

The `sbatch` command is given as a quoted string with `snakemake` style
wildcards for the actual values of the `sbatch` option values:

```
" sbatch -J {cluster.name} -A {cluster.account} -p {cluster.partition} -n {cluster.n} -t {cluster.time} {cluster.other}"
```

### Other very useful options

There are numerous options to snakemake, but I will only list some of the most
useful ones here:

- `-n` dry-run; _always_ run a dry-run before the real run to avoid unpleasant
surprises!

- `-p` print shell commands that will be performed -- good for pre-checking and
debugging rules

- `-r` state reason for rerunning a rule -- again good for pre-checking and
debugging

- `-f` force a rerun of the (last) rule that that has the requested file as
output -- sometimes needs doing!

### Wrapper script for Snakemake

When a lot of the `snakemake` options above are used, the command-line call to
 `snakemake` becomes quite cumbersome to type each time, e.g.,

```
snakemake -s path/to/Snakefile --use-conda --cluster-config pathe/to/cluster-config-file -- cluster " sbatch -J {cluster.name} -A {cluster.account} -p {cluster.partition} -n {cluster.n} -t {cluster.time} {cluster.other}" -p -r -n path/to/requested/outputfile
```

Therefore it is convenient and recommended to use a convenience *wrapper script*
that encapsulates, i.e., hides, all the absolutely necessary options. The
wrapper script is designed so that other, not encapsulated, options (e.g., `-n`)
can be passed to the wrapper script and it will be passed on to snakemake. This
wrapper script is then called instead of `snakemake`. An example wrapper script
can be found in [doMapping.sh](./doMapping.sh), which also contain some
additional convenience items.

If we use this wrapper script, the call above simplifies to

```
bash doMapping.sh -n path/to/requested/outputfile
```


## Exercises

In the exercises, the workflow defined in the Snakefile
[mapping.smk](/mapping.smk) is used together with this [Example
data](../ExampleData/README.md).

Throughout all exercises, take a look at the Snakefile
[mapping.smk](mapping.smk), the config files
[mappingConfig.yaml](mappingConfig.yaml) and
[mappingSamples.tsv](mappingSamples.tsv), and the env files, and try to
understand what they do and what happens when you run the commands in the
exercises.

0. If you have not already done so, fork the Crash Course Bitbucket repository
and then clone a git working directory (_gwd_), as described in the
[Git Crash course](../GitCrashCourse/README.md). To also appreciate
`snakemake`'s cluster capabilities you should do this on an UPPMAX cluster
login-node (e.g., rackham), but you can run on your laptop as well.
1. Create a *Analysis working directory* (`awd`) --- this should be different
and outside the git working directory (`gwd`) -- and `cd` into `awd`.

2. Run `snakemake` from `awd` to create the output file `fastq/s1_R1.fastq.gz`.
Do not run it as a cluster job or use the wrapper script `doMapping.sh` yet, but
think about what minimum options are needed. (Tip: It's good practice to do a
dry-run first to check check what happens and if all options are set correctly.)  
What files were created? Why?

3. Use the `doMapping.sh` wrapper script and run `snakemake` without a output
file argument.  
What samples did you get final output files for? Why?

4. Update the workflow so that final output files also for a sample named `s2`,
is produced; use input fastq `SRR3222412-19_1.fq.gz` and `SRR3222412-19_1.fq.gz`
in the [Example Data directory](../ExampleData/README.md) for `s_2`. Then rerun
the command used in 3. (Tip: Edit the *local* config files in `awd` -- not those
in `gwd`; use a dry-run to check if you got things right.)  
Did it work? Were all samples run?


### Extra-curricular exercise

5. Add a rule that uses the program `featureCount` to summarize the read counts
on the different features of the genome annotation file
[../ExampleData/Mus_musculus.GRCm38.99.chromosome.19.gtf.gz](../ExampleData/annotation/README.md).
The relevant shell command for featureCount would be:

        featureCounts -g gene_id -t exon -s 1 -R BAM -a {input.gtf} -o {output.counts} {input.bam}  
    Remember to also change the final output files in the target rule.  
    Did it work?
