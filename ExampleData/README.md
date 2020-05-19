#Example data to be used in the Snakemake Crash Course Exercises
The Example data is taken from the [RNAseq Analysis](https://nbisweden.github.io/workshop-ngsintro/2001/lab_rnaseq.html) lab of the NBIS NGS intro course (as given Jan 2020); The data was originally compiled by Roy Francis, NBIS. The following data description is taken from the RNAseq Analysis lab.

## Data description

The data used in this exercise is from the paper: **Poitelon, Yannick, *et al*. YAP and TAZ control peripheral myelination and the expression of laminin receptors in Schwann cells. [Nature neuroscience 19.7 (2016): 879](https://www.nature.com/articles/nn.4316)**. In this study, YAP and TAZ genes were knocked-down in Schwann cells to study myelination, using the sciatic nerve in mice as a model.

Myelination is essential for nervous system function. Schwann cells interact with neurons and the basal lamina to myelinate axons using receptors, signals and transcription factors. Hippo pathway is a conserved pathway involved in cell contact inhibition, and it acts to promote cell proliferation and inhibits apoptosis. The pathway integrates mechanical signals (cell polarity, mechanotransduction, membrane tension) and gene expression response. In addition to its role in organ size control, the Hippo pathway has been implicated in tumorigenesis, for example its deregulation occurs in a broad range of human carcinomas. Transcription co-activators YAP and TAZ are two major downstream effectors of the Hippo pathway, and have redundant roles in transcriptional activation.

The material for RNA-seq was collected from 2 conditions (**Wt** and **KO**), each with 3 biological replicates.

| Accession  | Condition | Replicate |
-------------|-----------|-----------|
| SRR3222409 | KO        | 1         |
| SRR3222410 | KO        | 2         |
| SRR3222411 | KO        | 3         |
| SRR3222409 | Wt        | 1         |
| SRR3222413 | Wt        | 2         |
| SRR3222414 | Wt        | 3         |


For the purpose of this tutorial, to shorten the time needed to run various bioinformatics steps, we have picked reads for a single chromosome (Chr 19) and downsampled the reads. We randomly sampled, without replacement, 25% reads from each sample, using `fastq-sample` from the toolset [fastq-tools](https://homes.cs.washington.edu/~dcjones/fastq-tools/).

In this crash course, the data has been further reduced by only using a single replicate of the first samples for each condition (i.e., SRR3222409 and SRR3222412).

The fastq files can be found in the [fastq](./fastq/README.md) directory.

## Reference genome files

It is best if the reference genome (.fasta) and annotation (.gtf) files come from the same source to avoid potential naming conventions problems.  For this exercise, we have choose ensembl version 99, and limited both reference and annotation to chromosome 19.

The reference fasta and annotation gtf files can be found in the [reference](reference/README.md) and [annotation](annotation/README.md) directories, respectively.
