#!/bin/bash
# pipeline.sh
# Usage: ./pipeline.sh <sample_id>

set -eo pipefail   # <--- FIXED: removed '-u' to prevent ADDR2LINE crashes

export CONDA_HOME="/home/ssingh6/miniconda3"
export PATH="$CONDA_HOME/bin:$PATH"

if [ -f "$CONDA_HOME/etc/profile.d/conda.sh"]; then
	source "$CONDA_HOME/etc/profile.da/conda.sh"
fi

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <sample_id>"
  exit 1
fi

SAMPLE_ID="$1"
BASE_DIR="/home/ssingh6/ssingh-AWS/Seqdata/${SAMPLE_ID}"
WORK_DIR="${BASE_DIR}/Pod5/pod5"
LOGFILE="${BASE_DIR}/pipeline_${SAMPLE_ID}.log"

REF_FA="/home/ssingh6/reference/hg38.fa"
REF_FA_FALLBACK="/home/ssingh6/reference/hg38.fa"

THREADS=50
SORT_THREADS=16

mkdir -p "${BASE_DIR}"

# Log stdout+stderr to console + file
exec > >(tee -a "$LOGFILE") 2>&1

echo "-------------------------------------------"
echo "Pipeline started at: $(date)"
echo "Sample ID: ${SAMPLE_ID}"
echo "Base directory: ${BASE_DIR}"
echo "Working directory: ${WORK_DIR}"
echo "Log file: ${LOGFILE}"
echo "-------------------------------------------"

cd "$WORK_DIR" || { echo "Failed to cd to ${WORK_DIR}"; exit 1; }

############################################
# Conda helper
############################################
source ~/miniconda3/etc/profile.d/conda.sh

activate_env() {
  local ENV_NAME="$1"
  echo "[ $(date) ] Activating env: ${ENV_NAME}"
  conda activate "${ENV_NAME}"
}

############################################
# Basecalling (skip if already exists)
############################################
BASECALL_BAM="${BASE_DIR}/${SAMPLE_ID}.hac_hmcg.bam"

echo "[ $(date) ] Basecalling step..."
if [ -s "${BASECALL_BAM}" ]; then
  echo "[ $(date) ] Found existing basecalled BAM: ${BASECALL_BAM}"
  echo "[ $(date) ] Skipping basecalling."
else
  echo "[ $(date) ] No basecalled BAM found; running dorado basecaller..."
  dorado basecaller hac,5mCG_5hmCG "${BASE_DIR}/Pod5/pod5/" > "${BASECALL_BAM}"
  echo "[ $(date) ] Basecalling completed."
fi

echo "[ $(date) ] Verifying basecalling output: ${BASECALL_BAM}"

if [ ! -s "${BASECALL_BAM}" ]; then
  echo "[ $(date) ] ERROR: Basecalling BAM missing or empty: ${BASECALL_BAM}"
  exit 10
fi

if ! samtools view -H "${BASECALL_BAM}" >/dev/null 2>&1; then
  echo "[ $(date) ] ERROR: Cannot read BAM header (corrupt?): ${BASECALL_BAM}"
  exit 11
fi

READS="$(samtools view -c "${BASECALL_BAM}")"
if [ "${READS}" -le 0 ]; then
  echo "[ $(date) ] ERROR: Basecalling BAM has zero reads: ${BASECALL_BAM}"
  exit 12
fi

echo "[ $(date) ] Basecalling OK (reads=${READS})."

############################################
# Reference select
############################################
if [ -f "${REF_FA}" ]; then
  USE_REF="${REF_FA}"
else
  echo "[ $(date) ] WARNING: ${REF_FA} not found; falling back to ${REF_FA_FALLBACK}"
  USE_REF="${REF_FA_FALLBACK}"
fi

############################################
# Alignment + Sorting + Indexing (skip if already exists)
############################################
ALIGNED_BAM="${BASE_DIR}/${SAMPLE_ID}.hac_hmcg.aligned.bam"
ALIGNED_SORT_BAM="${BASE_DIR}/${SAMPLE_ID}.aligned.sort.bam"

echo "[ $(date) ] Alignment+sorting step..."
if [ -s "${ALIGNED_SORT_BAM}" ] && [ -s "${ALIGNED_SORT_BAM}.bai" ]; then
  echo "[ $(date) ] Found existing aligned+sorted+indexed BAM: ${ALIGNED_SORT_BAM}"
  echo "[ $(date) ] Skipping alignment+sorting."
else
  echo "[ $(date) ] Starting Alignment (ref: ${USE_REF})..."
  dorado aligner "${USE_REF}" "${BASECALL_BAM}" > "${ALIGNED_BAM}"
  echo "[ $(date) ] Alignment completed."

  echo "[ $(date) ] Sorting aligned BAM..."
  samtools sort -@ "${SORT_THREADS}" -m 3G -O BAM \
    -o "${ALIGNED_SORT_BAM}" \
    "${ALIGNED_BAM}"

  samtools index -@ "${SORT_THREADS}" "${ALIGNED_SORT_BAM}"
  echo "[ $(date) ] Sorting+index completed."
fi

############################################
# Canonical contigs list (chr1-22,X,Y)
############################################
INCLUDE_CONTIGS="${BASE_DIR}/include_contigs.${SAMPLE_ID}.txt"
printf "chr%s\n" {1..22} X Y > "${INCLUDE_CONTIGS}"
echo "[ $(date) ] Using contigs list: ${INCLUDE_CONTIGS}"
head "${INCLUDE_CONTIGS}"

############################################
# Filter to canonical contigs
############################################
FILTERED_BAM="${BASE_DIR}/${SAMPLE_ID}.filtered.sort.chr.bam"

echo "[ $(date) ] Filtering alignments to canonical contigs..."
samtools view -@ "${THREADS}" -b "${ALIGNED_SORT_BAM}" \
  $(cat "${INCLUDE_CONTIGS}") \
  | samtools sort -@ "${SORT_THREADS}" -o "${FILTERED_BAM}" -

samtools index -@ "${SORT_THREADS}" "${FILTERED_BAM}"
echo "[ $(date) ] Alignment filtering completed."

############################################
# readCounter
############################################
echo "[ $(date) ] Running readCounter..."
readCounter --build "${FILTERED_BAM}"

readCounter --window 1000000 --quality 20 \
  --c "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY" \
  "${FILTERED_BAM}" \
  > "${BASE_DIR}/${SAMPLE_ID}.filtered.sort.wig"

echo "[ $(date) ] readCounter completed."

############################################
# ichorCNA
############################################
activate_env "Liquid_biopsy_CNV"

echo "[ $(date) ] Running ichorCNA..."
Rscript /home/ssingh6/ichorCNA/scripts/runIchorCNA.R --id tumor_sample \
  --WIG "${BASE_DIR}/${SAMPLE_ID}.filtered.sort.wig" \
  --ploidy "c(2,3)" --normal "c(0.5,0.6,0.7,0.8,0.9)" --maxCN 5 \
  --gcWig /home/ssingh6/ichorCNA/inst/extdata/gc_hg38_1000kb.wig \
  --mapWig /home/ssingh6/ichorCNA/inst/extdata/map_hg38_1000kb.wig \
  --centromere /home/ssingh6/ichorCNA/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt \
  --normalPanel /home/ssingh6/ichorCNA/inst/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds \
  --includeHOMD False --chrs "c(1:22, \"X\")" --chrTrain "c(1:22)" \
  --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
  --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 \
  --outDir "${BASE_DIR}"

echo "[ $(date) ] ichorCNA completed."
conda deactivate

############################################
# Variant Calling (Sniffles)
############################################
echo "[ $(date) ] Variant calling with Sniffles..."
sniffles -i "${FILTERED_BAM}" \
  -v "${BASE_DIR}/${SAMPLE_ID}.filtered.sort.vcf" \
  -t "${THREADS}"
echo "[ $(date) ] Variant calling completed."

############################################
# NanoDx
############################################
activate_env "nanodx"

echo "[ $(date) ] Running modkit pileup..."
modkit pileup "${ALIGNED_SORT_BAM}" \
  - --ref "${USE_REF}" --preset traditional --only-tabs --edge-filter 0,27 --threads "${THREADS}" \
  | gzip > "/home/ssingh6/ssingh-AWS/Seqdata/nanodx/nanoDx/results/${SAMPLE_ID}.CpG.bed.gz"

echo "[ $(date) ] Intersecting with 450K and formatting..."
zcat "/home/ssingh6/ssingh-AWS/Seqdata/nanodx/nanoDx/results/${SAMPLE_ID}.CpG.bed.gz" \
  | cut -f1-11 \
  | bedtools intersect -a - -b /home/ssingh6/ssingh-AWS/Seqdata/nanodx/nanoDx/static/450K_hg38.bed -wa -wb \
  | awk -v OFS='\t' '{$4=$15; print}' \
  | cut -f1-11 \
  > "/home/ssingh6/ssingh-AWS/Seqdata/nanodx/nanoDx/results/${SAMPLE_ID}.CpG.bed"

echo "[ $(date) ] Generating t-SNE plots..."
Rscript /home/ssingh6/ssingh-AWS/Seqdata/nanodx/nanoDx/plot_tSNE_standalone.R \
  "/home/ssingh6/ssingh-AWS/Seqdata/nanodx/nanoDx/results/${SAMPLE_ID}.CpG.bed" \
  /home/ssingh6/ssingh-AWS/Seqdata/nanodx/nanoDx/training_set/Capper_et_al.h5 \
  /home/ssingh6/ssingh-AWS/Seqdata/nanodx/nanoDx/static/colorMap_Capper_et_al.txt \
  "${BASE_DIR}/${SAMPLE_ID}-tSNE-Capper.pdf" \
  "${BASE_DIR}/${SAMPLE_ID}-tSNE-Capper.html" \
  tsne 94 30 2500 10 0.5 100

echo "[ $(date) ] NanoDx plotting completed."
conda deactivate

activate_env "NN_model"
echo "[ $(date) ] Running neural network classifier..."
python /home/ssingh6/ssingh-AWS/Seqdata/nanodx/nanoDx/nn_classifier.py \
  /home/ssingh6/ssingh-AWS/Seqdata/nanodx/nanoDx/training_set/Capper_et_al_NN.pkl \
  "/home/ssingh6/ssingh-AWS/Seqdata/nanodx/nanoDx/results/${SAMPLE_ID}.CpG.bed" \
  "/home/ssingh6/ssingh-AWS/Seqdata/nanodx/nanoDx/results/${SAMPLE_ID}_votes.tsv" \
  "/home/ssingh6/ssingh-AWS/Seqdata/nanodx/nanoDx/results/${SAMPLE_ID}_summary.txt"
conda deactivate

echo "-------------------------------------------"
echo "Pipeline finished at: $(date)"
echo "Log file saved as: ${LOGFILE}"
echo "-------------------------------------------"

