# Introduction to Snakemake
`snakemake` is a tool for designing and performing a workflows i.e.,
a sequence of individual steps, each running programs or scripts, that is
required to convert some input, say, sequencing reads, to a desired output,
say, read counts for individual genomic regions of interest.

In one sense, a snakemake workflow can be viewed as just a glorified bash
script performing the desired steps in sequence, but it provides a more
structured layout and is more flexible and easier to rerun/reuse.

In another sense, snakemake could be described as "_python_ with inspiration
from the _make_ program". However, the basic snakemake syntax is not so very
much python... and you definitely don't need to be a python expert to start
using snakemake. However, you can use python code snippets _almost_
anywhere in the Snakefile (see examples in SlamdunkSnakefile).

# The Snakefile
The cornerstone of snakemake is the Snakefile, which defines the workflow
steps. Each step is coded as a (usually) named `rule` and comprise a
number of `directives`, described below. Additionally, commands defining
configure-files, import of python modules or definition of python
functions can be included in the Snakefile.

## Directives:
Directives The most important directives are the `output`, `input` and `shell`
directives and are defined as follows (although here heavily documented):

```
# Notice the colon (:) after the rule and the directives
# (python style definitions)
rule linkFasta:
    # the output directive is a (possibly nested) list of desired
    # output files, the list items can be named. The output directive
    # defines wildcards inside curly brackets {}; the value of the
    # wildcards is determined when the rule is called either from
    # command-line or by another rule (more about this below).
    output:
        fastq = "Slamdunk/fastq/{sample}.fastq"
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
  programs required by the rule (see examples in SlamdunkSnakefile).
- `group:` assigns a "cluster job group" to a rule, see below.

Examples on how these are used can be found in the SlamdunkSnakefile.

## Wildcards and connecting rules
If the a Snakefile with the rule above is run as

`snakemake Slamdunk/fastq/mysample.fastq`

snakemake will first check that the file does not already exists,
in which case it reports this and stops. Otherwise, it will parse
the Snakefile to find a rule with output that
matches `Slamdunk/fastq/mysample.fastq`. It will see that the rule
`linkFasta` matches if wildcard `"{sample}"` is defined to
be  `"mysample"`, so it defines the wildcard and run the rule.

It will then look for its input file, which after expanding the
wildcard becomes `"path/to/originals/mysample.fastq"`; if this file
exists, snakemake will use this file as the input, Otherwise, it will
look for another rule with a matching output and run this other rule
to create the input. If you study the workflow in SlamdunkSnakefile,
you will see that different rules have matching output and input;
this creates a sequence of rules to be run -- a workflow.

## Rule order and the Target rule
If snakemake is run without file argument, i.e., `snakemake` the
first rule in the Snakefile is by default run. Therefore, the first
rule is often designed as a, so-called, _target_ rule that only
has an input directive comprising the desired final output files
(see example in SlamdunkSnakefile).

With the exception of the Target rule, rules can be written in
any order in the snakefile. I tend to have them in 'chronological'
order, i.e., the order they are expected to be executed.

The workflow in SlamdunkSnakefile is set up to perform all steps
(dunks) of the Slamdunk program for one single sample. We can
run the full workflow for a single sample, _s1_, by typing
`snakemake Slamdunk/filtered/s1_slamdunk_mapped_filtered_snp_tcount.csv`,
i.e, with the output file of rule `count` for sample s1.
However, by listing the output files from the count rule for
all samples in the Target rule, just running `snakemake` will
run the full workflow for all samples.

## Cluster runs, groups, and localrules
When run on a cluster (see below and in `doSlamdunk.sh`), by default,
each rule will be run as a separate job. However, if a few sequential
rules each takes less than, say, 15 min and all uses the same
amount of resources (memory or cpus/threads), then these could
be sent as one single job by assigning the same `group` using
the directive `group`, see example below. Moreover, if some jobs
are ridiculously short (e.g., softlinking external files), these
can be run on the login node by listing them as `localrules`.

## Comments
Comments can be added in either of the following ways:

- `\#...` Everything after a hash on a line is ignored by snakemake
- `"""..."""` Everything between triplicated double quotes will be ignored
by snakemake. However, snakemake will still try to expand strings
between curly brackets as variables, which can lead to errors.
- `'''...'''` Everything between triplicated double quotes will be
ignored by snakemake. This will prevent variable expansions,
so use this to comment out code etc.



# The configure files
It is convenient to define variables, including external input data,
in a configure file.

## _Standard_ configure files in yaml or json format
The _standard_ configure files can be in yaml-format or json format.
An advantage of the yaml format is that it is quite similar to the
directive syntax in snakemake; however, a really annoying drawback
is that indentations _must be as spaces and not as ASCII tab-characters_.
Many, but not all, text editors handles this correctly, however,
i.e., when you hit the tab-key it will insert a sequence of spaces
instead of the ASCII tab character (I use emacs with a yaml mode;
  atom is quite nice and at least have a function for converting tabs to spaces).
The json format is more strongly typed; it uses brackets and commas
instead of indentations and requires all strings to be quoted;
also, it does not support comments.

- yaml example

```yaml
# to assign a single value to variable
variable1: 17
# to assign a list
variable2:
  - value1
  - value2
# to assigna dictionary
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

If the configure file location is defined in the Snakefile, either
of the examples will be expand to the same python dictionary, which
can be accessed by `config["variablename"]` (see further
  _SlamdunkConfig.yaml_").

## _Tabular_ configure files in spreadsheet text format

If we have a lot of samples, for which we each want to include a number
of configuration variables, the standard configur file can become rather
messy and cumbersome. In such a case, it may be better to use a
_tabular_ configure file. This is essentially a spreadsheet, but it most be
in a text file format (not MS Excel, MacOSX numbers, etc.). This text file
format could be tab-delimited (tsv) or comma-demimited (csv); `snakemake`
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

# Running Snakemake
Simply typing `snakemake` on the command line (given that snakemake
is installed on the system), will start snakemake and look for a
Snakefile in the present directory. However, it is good keep the
code and the results in separate folders, so we usually want to
tell snakemake where to find the Snakefile. This is done by the
option `-s`:
```
snakemake -s path/to/Slamdunkfile
```

As mentioned above, when run like this, i.e., without a file argument,
snakemake will execute the first rule in the Snakefile. If this is
designed as a _target file_ (see above), typically the whole workflow
will be executed.

However, often one only wants to perform a single rule (or all rules
up to that rule) for a single sample. This can be achieved by running
snakemake with the output of the rule and sample as a file argument:
```
snakemake -s path/to/Slamdunkfile path/to/requested/outputfile
```

## Running on a Cluster
__Important__: Snakemake should always be started in a _virtual terminal_
(see further the TmuxCrashCourse.md) on the _login node_. Snakemake will
then, based on the Snakefile, determine jobs and send them to the cluster.
You do not need to manually submit jobs!

To run snakemake in parallel on UPPMAX, you use the following options:
- `--cluster` this defines the sbatch command needed to commit jobs to the
cluster, and uses rule-specific variables defined in a cluster-config file.
- `--cluster-config` this tells snakemake where to look for the
cluster-config file (see further SlamdunkCluster.yaml)

Details on how to run snakemake can be found in the wrapper script `doSlamdunk.sh`

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

## Wrapper script for Snakemake
If using a lot of options, the command line call to snakemake can be
quite long and cumbersome, especially when running on a cluster.
Therefore, I usually create a wrapper script that hides all
the absolutely necessary options. I design the wrapper script so that
I can pass other options (e.g., `-n`) to the wrapper script and it
will be passed on to snakemake. See further the script doSlamdunk.sh.
