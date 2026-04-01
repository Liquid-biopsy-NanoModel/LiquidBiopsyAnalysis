#!/bin/sh
#SBATCH -J "Nanopore_exp1.baseCall.align"
#SBATCH -o Nanopore_exp1.baseCall.align.out
#SBATCH -e Nanopore_exp1.baseCall.align.err
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100:1
#SBATCH --time=3-00:00:00
#SBATCH --mem=200G
#SBATHC --cpus-per-task=55
#module load cuda/11.7
module load samtools
export PATH=/home/ssingh6/dorado-0.7.3-linux-x64/bin:$PATH
export PATH=/home/ssingh6/hmmcopy_utils/bin:$PATH
bash Liquid_biopsy_pipeline_renewed.sh 28_Zurich
