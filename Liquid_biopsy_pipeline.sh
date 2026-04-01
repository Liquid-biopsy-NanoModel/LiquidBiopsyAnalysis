#!/bin/bash
# pipeline.sh
#
# Usage: ./pipeline.sh <sample_id>
#
# This script executes a pipeline that includes basecalling, alignment,
# filtering, read counting, CNV analysis, variant calling, Nano dx processing,
# and neural network classification.
#
# All file paths are dynamically created using the supplied sample ID.
# If a required conda environment is missing, it is built from a YAML file found
# in the sample ID base folder (named as <env_name>.yml).
# The process is logged with timestamps.

# Exit immediately on any error.
set -e

# Ensure a sample ID is provided.
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <sample_id>"
  exit 1
fi

SAMPLE_ID=$1
BASE_DIR="/home/ssingh6/ssingh-AWS/Seqdata/${SAMPLE_ID}"
WORK_DIR="${BASE_DIR}/Pod5/pod5"
LOGFILE="${BASE_DIR}/pipeline_${SAMPLE_ID}.log"

# Redirect stdout and stderr to the logfile and console.
exec > >(tee -a "$LOGFILE") 2>&1

echo "-------------------------------------------"
echo "Pipeline started at: $(date)"
echo "Sample ID: ${SAMPLE_ID}"
echo "Base directory: ${BASE_DIR}"
echo "Working directory: ${WORK_DIR}"
echo "-------------------------------------------"

# Change to the working directory.
cd "$WORK_DIR" || { echo "Failed to change directory to ${WORK_DIR}"; exit 1; }

############################################
## Basecalling, Alignment, and Sorting    ##
############################################

echo "[ $(date) ] Starting Basecalling step..."
dorado basecaller hac,5mCG_5hmCG ${BASE_DIR}/Pod5/pod5/ \
  > ${BASE_DIR}/${SAMPLE_ID}.hac_hmcg.bam
echo "[ $(date) ] Basecalling completed."

echo "[ $(date) ] Starting Alignment step..."
dorado aligner /home/ssingh6/reference/hg38.fa \
  ${BASE_DIR}/${SAMPLE_ID}.hac_hmcg.bam \
  > ${BASE_DIR}/${SAMPLE_ID}.hac_hmcg.aligned.bam
echo "[ $(date) ] Alignment completed."

echo "[ $(date) ] Sorting BAM file..."
samtools sort -m 3G -O BAM \
  -o ${BASE_DIR}/${SAMPLE_ID}.hac_hmcg.aligned.sort.bam \
  ${BASE_DIR}/${SAMPLE_ID}.hac_hmcg.aligned.bam
samtools index ${BASE_DIR}/${SAMPLE_ID}.hac_hmcg.aligned.sort.bam
echo "[ $(date) ] Sorting completed."

echo "[ $(date) ] Basecalling process ends!"
echo "Process Ends!"

#################################################
## Alignment Filtering and Read Count Analysis ##
#################################################

echo "[ $(date) ] Filtering alignments..."
# Extract headers separately, then filter alignments
samtools idxstats  ${BASE_DIR}/${SAMPLE_ID}.hac_hmcg.aligned.sort.bam | \
  cut -f1 | grep -v -e 'chrM' -e 'chrUn_KI' -e 'chrEBV' -e 'chrUn_GL' \
             -e 'chr1_KI' -e 'chr2_KI' -e 'chr3_KI' -e 'chr4_KI' \
             -e 'chr5_KI' -e 'chr9_KI' -e 'chr11_KI' -e 'chr14_KI' \
             -e 'chr15_KI' -e 'chr16_KI' -e 'chr17_KI' -e 'chr22_KI' \
             -e 'chrY_KI' -e 'chr17_GL' > include_contigs.txt

samtools view -@ 50 -b  ${BASE_DIR}/${SAMPLE_ID}.hac_hmcg.aligned.sort.bam $(cat include_contigs.txt) \
    >  ${BASE_DIR}/${SAMPLE_ID}.hac_hmcg.aligned.filtered.sort.bam
samtools index ${BASE_DIR}/${SAMPLE_ID}.hac_hmcg.aligned.filtered.sort.bam

echo "[ $(date) ] Alignment filtering completed."

echo "[ $(date) ] Alignment assessment begins. "
samtools view -h ${BASE_DIR}/${SAMPLE_ID}.hac_hmcg.aligned.sort.bam | sed -E 's/^@SQ\tSN:([0-9]+|X|Y)/@SQ\tSN:chr\1/; s/^([^\t]+\t[0-9]+\t)([0-9]+|X|Y)/\1chr\2/' | samtools view -b - \ > ${BASE_DIR}/${SAMPLE_ID}.hac_hmcg.aligned.sort.chr.bam

echo "[ $(date) ] Running readCounter analysis..."
readCounter --build ${BASE_DIR}/${SAMPLE_ID}.hac_hmcg.aligned.filtered.sort.chr.bam
readCounter --window 1000000 --quality 20 --c "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY" \
  ${BASE_DIR}/${SAMPLE_ID}.hac_hmcg.aligned.filtered.sort.chr.bam \
  > ${BASE_DIR}/${SAMPLE_ID}.hac_hmcg.aligned.filtered.sort.wig
echo "[ $(date) ] readCounter analysis completed."

######################################
## ichorCNV Analysis (with liquid_biopsy_CNV) ##
######################################

# Source conda so that conda environments can be activated.
#source /home/ec2-user/anaconda3/etc/profile.d/conda.sh
source ~/miniconda3/etc/profile.d/conda.sh
# Function: Activate a conda environment if it exists;
# otherwise, create it from a YAML file in the base directory.
activate_env() {
  ENV_NAME=$1
  if conda info --envs | awk '{print $1}' | grep -q "^${ENV_NAME}$"; then
    echo "[ $(date) ] Conda environment '${ENV_NAME}' found. Activating..."
    conda activate "${ENV_NAME}"
  else
    echo "[ $(date) ] Conda environment '${ENV_NAME}' not found. Creating from YAML and activating..."
    # Expect the YAML file to be located in BASE_DIR and named as <ENV_NAME>.yml
    conda env create -n "${ENV_NAME}" -f "${BASE_DIR}/${ENV_NAME}.yml"
    conda activate "${ENV_NAME}"
  fi
}

# Activate the 'liquid_biopsy_CNV' environment before running ichorCNV.
activate_env "Liquid_biopsy_CNV"

echo "[ $(date) ] Running ichorCNV analysis..."
Rscript /home/ssingh6/ichorCNA/scripts/runIchorCNA.R --id tumor_sample \
  --WIG ${BASE_DIR}/${SAMPLE_ID}.hac_hmcg.aligned.filtered.sort.wig \
  --ploidy "c(2,3)" --normal "c(0.5,0.6,0.7,0.8,0.9)" --maxCN 5 \
  --gcWig /home/ssingh6/ichorCNA/inst/extdata/gc_hg38_1000kb.wig \
  --mapWig /home/ssingh6/ichorCNA/inst/extdata/map_hg38_1000kb.wig \
  --centromere /home/ssingh6/ichorCNA/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt \
  --normalPanel /home/ssingh6/ichorCNA/inst/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds \
  --includeHOMD False --chrs "c(1:22, \"X\")" --chrTrain "c(1:22)" \
  --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
  --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 \
  --outDir ${BASE_DIR}

echo "[ $(date) ] ichorCNV analysis completed."


##########################
## Variant Calling      ##
##########################

echo "[ $(date) ] Performing Variant Calling..."
sniffles -i ${BASE_DIR}/${SAMPLE_ID}.hac_hmcg.aligned.filtered.sort.bam \
  -v ${BASE_DIR}/${SAMPLE_ID}.hac_hmcg.aligned.filtered.sort.vcf -t 50
#python3 -m sniffles2_plot \
#  -i ${BASE_DIR}/${SAMPLE_ID}.hac_hmcg.aligned.filtered.sort.vcf \
#  -o ./tumor_sample
echo "[ $(date) ] Variant Calling completed."

# Deactivate the 'liquid_biopsy_CNV' environment.

conda deactivate
##########################
## Nano dx Analysis     ##
##########################

echo "[ $(date) ] Starting Nano dx analysis..."

# Activate the 'nanodx' environment.
activate_env "nanodx"

echo "[ $(date) ] Running modkit pileup for Nano dx..."
modkit pileup ${BASE_DIR}/${SAMPLE_ID}.hac_hmcg.aligned.sort.chr.bam \
  - --ref /home/ssingh6/reference/hg38.fa --preset traditional --only-tabs --edge-filter 0,27 --threads 50 \
  | gzip > /home/ssingh6/ssingh-AWS/Seqdata/nanodx/nanoDx/results/${SAMPLE_ID}.CpG.bed.gz

echo "[ $(date) ] Processing pileup output and filtering with bedtools..."
zcat /home/ssingh6/ssingh-AWS/Seqdata/nanodx/nanoDx/results/${SAMPLE_ID}.CpG.bed.gz | cut -f1-11 | \
  bedtools intersect -a - -b /home/ssingh6/ssingh-AWS/Seqdata/nanodx/nanoDx/static/450K_hg38.bed -wa -wb | \
  awk -v OFS='\t' '{$4=$15; print}' | cut -f1-11 \
  > /home/ssingh6/ssingh-AWS/Seqdata/nanodx/nanoDx/results/${SAMPLE_ID}.CpG.bed

echo "[$(BASE_DIR)] Generating t-SNE plots..."
Rscript /home/ssingh6/ssingh-AWS/Seqdata/nanodx/nanoDx/plot_tSNE_standalone.R \
  /home/ssingh6/ssingh-AWS/Seqdata/nanodx/nanoDx/results/${SAMPLE_ID}.CpG.bed \
  /home/ssingh6/ssingh-AWS/Seqdata/nanodx/nanoDx/training_set/Capper_et_al.h5 \
 /home/ssingh6/ssingh-AWS/Seqdata/nanodx/nanoDx/static/colorMap_Capper_et_al.txt \
  ${BASE_DIR}/${SAMPLE_ID}-tSNE-Capper.pdf \
  ${BASE_DIR}/${SAMPLE_ID}-tSNE-Capper.html \
  tsne 94 30 2500 10 0.5 100
echo "[ $(date) ] t-SNE plot generation completed."

# Deactivate the 'nanodx' environment.
conda deactivate

# Activate the 'NN_model' environment for neural network classification.
activate_env "NN_model"

echo "[ $(date) ] Running neural network classifier..."
python /home/ssingh6/ssingh-AWS/Seqdata/nanodx/nanoDx/nn_classifier.py \
  /home/ssingh6/ssingh-AWS/Seqdata/nanodx/nanoDx/training_set/Capper_et_al_NN.pkl \
  /home/ssingh6/ssingh-AWS/Seqdata/nanodx/nanoDx/results/${SAMPLE_ID}.CpG.bed \
  /home/ssingh6/ssingh-AWS/Seqdata/nanodx/nanoDx/results/${SAMPLE_ID}_votes.tsv \
  /home/ssingh6/ssingh-AWS/Seqdata/nanodx/nanoDx/results/${SAMPLE_ID}_summary.txt

#Rscript ${BASE_DIR}/${SAMPLE_ID}/focal_amplifier.R /home/ec2-user/Seqdata/nanodx/nanoDx/results/${SAMPLE_ID}.CpG.bed ${BASE_DIR}/${SAMPLE_ID}/tumor_sample.cna.seg 
conda deactivate
echo "[ $(date) ] Deactivated 'NN_model' environment."

##########################
## Pipeline Completed   ##
##########################

echo "-------------------------------------------"
echo "Pipeline finished at: $(date)"
echo "Log file saved as: ${LOGFILE}"
echo "-------------------------------------------"

