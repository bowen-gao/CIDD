#!/bin/bash

# =============================
# Script to run UniMAP Feature Extraction
# =============================

# Variables
CONDA_ENV_PATH="/mnt/nfs-ssd/data/conda_envs/UniMOL"
CUDA_LIB_PATH="/usr/local/cuda/lib64"
CONDA_LIB_PATH="/mnt/nfs-ssd/data/conda_envs/UniMOL/lib"
PROJECT_DIR="/mnt/nfs-ssd/data/fengshikun/UniMAP/iupac-pretrain"
MODEL_PATH="/mnt/nfs-ssd/data/fengshikun/UniMAP/model/train_1kw_gcn_3_8gpu_check_frag/final/pytorch_model.bin"
DATASET_PATH="all_feats/temp_frags.smi"
ATOM_VOCAB_FILE="./Merge_vocab.pkl"
TOKENIZER_PATH="./iupac_regex"
SMILES_TOKENIZER_PATH="./smiles_tokenizer"
OUTPUT_DIR="temp_smiles_only"
LOG_FILE="extract_feature_uni.log"
OUTPUT_PATH="/mnt/nfs-ssd/data/gaobowen/CIDD/temp_frags/"

# Hyperparameters
PER_DEVICE_EVAL_BATCH_SIZE=256
POOLER_TYPE="avg"
GNN_NUMBER_LAYER=3

# CUDA settings
CUDA_DEVICE=2

# Activate Conda environment
echo "Activating Conda environment: $CONDA_ENV_PATH"
source "$(conda info --base)/etc/profile.d/conda.sh"  # Ensure `conda` command is available
conda activate $CONDA_ENV_PATH

if [[ $? -ne 0 ]]; then
  echo "Failed to activate Conda environment: $CONDA_ENV_PATH"
  exit 1
fi

# Update LD_LIBRARY_PATH
export LD_LIBRARY_PATH="$CUDA_LIB_PATH:$CONDA_LIB_PATH:$LD_LIBRARY_PATH"

# Navigate to project directory
cd "$PROJECT_DIR" || {
  echo "Failed to change directory to $PROJECT_DIR"
  exit 1
}

# Execute Python script
echo "Running feature extraction..."
CUDA_VISIBLE_DEVICES=$CUDA_DEVICE python -u extract_feature_uni.py \
    --run_name MultiUNI_modal_Mol \
    --model_path "$MODEL_PATH" \
    --dataset_path "$DATASET_PATH" \
    --smiles_only \
    --get_frag \
    --atom_vocab_file "$ATOM_VOCAB_FILE" \
    --atom_vocab_size 10535 \
    --tokenizer_path "$TOKENIZER_PATH" \
    --smiles_tokenizer_path "$SMILES_TOKENIZER_PATH" \
    --output_dir "$OUTPUT_DIR" \
    --per_device_eval_batch_size $PER_DEVICE_EVAL_BATCH_SIZE \
    --pooler_type $POOLER_TYPE \
    --gnn_number_layer $GNN_NUMBER_LAYER | tee "$LOG_FILE" \
    #--out_lmdb_path $OUTPUT_PATH

# Check execution status
if [[ $? -eq 0 ]]; then
  echo "Feature extraction completed successfully. Logs saved to $LOG_FILE"
else
  echo "Feature extraction failed. Check logs in $LOG_FILE for details."
  exit 1
fi