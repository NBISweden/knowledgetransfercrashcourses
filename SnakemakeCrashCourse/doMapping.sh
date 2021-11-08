#!/usr/bin/env bash
set -e

# get the workflow base directory (=where this script reside)
DIR=${0%doMapping.sh}

# Allow for either running on cluster or on local computer/laptop
if [ "$CLUSTER" = "rackham" ]; then # Change/add cluster name if needed
    echo "Running on rackham"
    # load relevant conda modules
    module load bioinfo-tools
    # This makes snakemake use the uppmax local conda repos.
    # These seems, currently, to be slower than conda over the net,
    # (possibly due to a smakemake issue), but are essential when
    # running local only (i.e., on bianca). Must be loaded before
    # the snakemake module.
#    module load conda # disabled for now, as it appears slow
    module load snakemake/5.10.0

    # To run on cluster, a cluster config file is needed. If not present,
    # copy from workflow base directory. Update the contact email at
    # the same time.
    if [[ ! -f MappingCluster.yaml ]]; then
	echo "No MappingCluster.yaml found in present directory"
	echo "Do you want me to set this directory up as a analysis directory (y/n)?"
	read answer
	if [[ "$answer" != "y" ]]; then
	    exit
	fi
	echo "What's your email (for UPPMAX job error notifications)?"
	read email
	awk -v email=$email '{sub("your@email", email, $0); print $0}' $DIR/MappingCluster.yaml > MappingCluster.yaml
	echo "MappingCluster.yaml copied with email updated. Your can make further changes to this file manually."
    fi

    time snakemake \
      --snakefile $DIR/MappingSnakefile.smk \
      --use-conda \
      --cores 1 \
      --cluster-config MappingCluster.yaml \
      --cluster " sbatch -J {cluster.name} -A {cluster.account} -p {cluster.partition} -n {cluster.n} -t {cluster.time} {cluster.other} " \
      -j \
      $@
   # Explanation of command and options:
   # `time`           just gives execution time report
   # -s               tells which Snakefile to use
   # --use-conda      tells snakemake to create and use rule-specific conda envs
   # --cluster-config tells what cluster-config file to use
   # --cluster " sbatch -J {cluster.name} ..."
   #                  This is the sbatch command with variables
   #                  taken from the cluster-config
   # -j               Tells snakemake all cores available for the job
   # $@               pass additional arguments (requested outout file, options)
   #                  to snakemake, e.g. `doMapping.sh -n myref.fasta` will add
   #                  `-n myref.fasta` directly to the end of the snakemake call
   # \                `line break`, marks that the command continues on the
   #                  next line (NB! no space or text after the `\`)
else
    # When run locally, we don't need --cluster-config or --cluster
    time snakemake -s $DIR/MappingSnakefile.smk --use-conda --cores 1 $@
fi


# How snakemake handles the cluster calls:
# 1. Snakemake should be run on the login node in a virtual terminal (use tmux -- recommended -- or screen)
# 2. Snakemake first parses the Snakefile and create the workflow as a directed acyclic graph (DAG)
# 3. If needed it creates the requested conda environments
# 4. It parses the --cluster-config file and reads a python dictionary indexed by rule names and __default__
# 5. For each rule to be run, it the --cluster sbatch command and substitutes the {cluster.xxx}
#    variables for corresponding values from the --cluster-config for the rule
# 6. It sends the job to the slurm queue -- when possible it will run rules/jobs in parallel
# 7. When job is started, it first creates any conda-environment that is needed, and then runs the rule.
# 8. When the job returns, snakemake checks it it has failed or succeeded and then decide this
#    means that it should start another job/rule that depends on the output of the one just run.
