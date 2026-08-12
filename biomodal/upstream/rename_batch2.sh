# biomodal/upstream/rename_batch2.sh
# Rename batch 2 FASTQ files and move to evoC-run input directory

SEQ_DATA="/expanse/lustre/projects/csd940/zalibhai/biomodal/seq-data"
INPUT_DIR="/expanse/lustre/projects/csd940/zalibhai/biomodal/evoC-run/input"

echo "===== Rename & Move Batch 2 FASTQs ====="
echo "From: $SEQ_DATA"
echo "To:   $INPUT_DIR"
echo ""

# ctrl-M (single underscore in source filename)
echo "Moving ctrl-M..."
mv "$SEQ_DATA/260311_evoC_Bap1_ctrl_M_S1_L004_R1_001.fastq.gz" \
   "$INPUT_DIR/evoC-Bap1-ctrl-M-B2_S1_L004_R1_001.fastq.gz"
mv "$SEQ_DATA/260311_evoC_Bap1_ctrl_M_S1_L004_R2_001.fastq.gz" \
   "$INPUT_DIR/evoC-Bap1-ctrl-M-B2_S1_L004_R2_001.fastq.gz"

# mut-M (single underscore in source filename)
echo "Moving mut-M..."
mv "$SEQ_DATA/260311_evoC_Bap1_mut_M_S2_L004_R1_001.fastq.gz" \
   "$INPUT_DIR/evoC-Bap1-mut-M-B2_S2_L004_R1_001.fastq.gz"
mv "$SEQ_DATA/260311_evoC_Bap1_mut_M_S2_L004_R2_001.fastq.gz" \
   "$INPUT_DIR/evoC-Bap1-mut-M-B2_S2_L004_R2_001.fastq.gz"

# ctrl-F (double underscore in source filename)
echo "Moving ctrl-F..."
mv "$SEQ_DATA/260311__evoC_Bap1_ctrl_F_S3_L004_R1_001.fastq.gz" \
   "$INPUT_DIR/evoC-Bap1-ctrl-F-B2_S3_L004_R1_001.fastq.gz"
mv "$SEQ_DATA/260311__evoC_Bap1_ctrl_F_S3_L004_R2_001.fastq.gz" \
   "$INPUT_DIR/evoC-Bap1-ctrl-F-B2_S3_L004_R2_001.fastq.gz"

# mut-F (double underscore in source filename)
echo "Moving mut-F..."
mv "$SEQ_DATA/260311__evoC_Bap1_mut_F_S4_L004_R1_001.fastq.gz" \
   "$INPUT_DIR/evoC-Bap1-mut-F-B2_S4_L004_R1_001.fastq.gz"
mv "$SEQ_DATA/260311__evoC_Bap1_mut_F_S4_L004_R2_001.fastq.gz" \
   "$INPUT_DIR/evoC-Bap1-mut-F-B2_S4_L004_R2_001.fastq.gz"

echo ""
echo "===== Verifying input directory ====="
echo ""
echo "Batch 1 (L000):"
ls -lh "$INPUT_DIR"/*L000*.fastq.gz 2>/dev/null
echo ""
echo "Batch 2 (L004):"
ls -lh "$INPUT_DIR"/*L004*.fastq.gz 2>/dev/null
echo ""
echo "Total FASTQ files:"
ls "$INPUT_DIR"/*.fastq.gz | wc -l
echo ""
echo "Done."
