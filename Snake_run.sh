#!/usr/bin/env bash

#SBATCH --time=0-12:00:00
#SBATCH --nodes=2
#SBATCH --mem=190000M
#SBATCH -o logs/PK_EDTA_Step1-%j.out
#SBATCH -e logs/PK_EDTA_Step1-%j.err
#SBATCH --account=rpp-rieseber

eval "$(conda shell.bash hook)"

conda activate snakemake

snakemake --cluster "sbatch -A {cluster.account} -t {cluster.time} --mem {cluster.mem} -N {cluster.nodes} --ntasks-per-node {cluster.ntask} --cpus-per-task {cluster.cpus} -o {cluster.output} -e {cluster.error}" --cluster-config cluster_config.yaml --jobs 24

conda deactivate
