# Reproducibility Crash Course

**WORK IN PROGRESS**  
**Todo: see also subfolder issues**  

- Maybe a separate course/page on yaml, since it pops up both for conda and
snakemake?

**Feedback: please add suggestions to this and subfolder issues**

## Prerequisites
Experience of working in a bash terminal. Some of the material refers to usage
on UPPMAX, Hence, while not absolutely necssary, you will benefit best if you
have some experience of working on UPPMAX.

## Introduction

This repo contains a couple of crash courses to common tools used in NBIS
projects to ensure reproducibility and documentation. The target audience are
the recipient users in NBIS projects.

This set of crash courses are focused around the workflow tool `snakemake`,
 which together with `git`, `conda` and `tmux` help researchers to:

- Document the different data handling and analysis steps in a project, and the
order they are performed, in a structured way as a workflow.
- Avoid manual changes to data files, and  use the workflow (or an associated
script) to perform the required changes instead.
- Document the software and system environment needed for the workflow.
- Store workflow code, environment information and other documentation
documentation and  and track changes made to them in a code repository.
- Download existing workflow code to another computer system, e.g. a cluster,
and rerunning it on data.

The crash course covers the bare essential commands and knowledge to start using
 `snakemake` workflows and is hoped to constitute a basis for further own
 learning.


### 1. Environment management

Reproducing a workflow may require installing several software tools and
libraries. `conda` is an open source environment management system that greatly
enhances (re)creating a software environment, e.g., required by a workflow.
Take the [Conda Crash Course](/CondaCrashCourse/README.md) (to be added) to
learn how to use it.

An alternative way to manage environments, that provide an even higher level of
reproduciblility by providing also the operating system, are *containers* such
as `docker` and `singularity`. These are *not* covered in this crash course.

### 2. Version control

A version control system stores and keep track of your workflow code and forms
the basis for reproducible research in NBIS. The most commonly used version
control system is `git`. Take the [Git Crash Course](/GitCrashCourse/README.md)
to learn how to use it.

### 3. Virtual terminals

When running workflows on a cluster, a *virtual terminal*, such as `tmux`, that
allows you to start the workflow and then log out is extremely useful. Take the
[Tmux Crash Course](/TmuxCrashCourse/README.md) to learn how to use it.

### 4. Workflow tools

A workflow tool, such as `snakemake` allows a structured way to construct,
document and run (and rerun) an analysis project.`snakemake` automatizes running
multiple samples, running parallel jobs on a cluster, using conda, and more.
Take the [Snakemake Crash course](/SnakemakeCrashCourse/README.md) to learn how
to use it.

An alternative workflow tool is `nextflow`, which is *not* covered here. For
both `snakemake` and `nextflow` collections of ready-made workflows are
available at [snakemake-workflows](https://github.com/snakemake-workflows) and
[nf-core](https://nf-co.re/pipelines).
