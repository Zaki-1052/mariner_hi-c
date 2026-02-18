LAB_BL="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/tads/250123blacklist.bed"
ENCODE_BL="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/tads/mm10-blacklist.v2.bed"
MERGED_BL="/expanse/lustre/projects/csd940/zalibhai/abc/reference/mm10_merged_blacklist.bed"

# Verify both exist
ls -la "${LAB_BL}" "${ENCODE_BL}"

# Check column format consistency
echo "Lab blacklist format (first 3 lines):"
head -3 "${LAB_BL}"
echo "ENCODE blacklist format (first 3 lines):"
head -3 "${ENCODE_BL}"

# Both must have at minimum chr/start/end in columns 1-3
# Extract BED3 from both, cat, sort, merge
cat <(cut -f1-3 "${LAB_BL}") <(cut -f1-3 "${ENCODE_BL}") \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - \
  > "${MERGED_BL}"

echo "Merged blacklist: $(wc -l < "${MERGED_BL}") regions"
echo "  Lab-only:   $(bedtools subtract -a "${LAB_BL}" -b "${ENCODE_BL}" | wc -l) unique to lab"
echo "  ENCODE-only: $(bedtools subtract -a "${ENCODE_BL}" -b "${LAB_BL}" | wc -l) unique to ENCODE"
echo "  Overlapping: $(bedtools intersect -a "${LAB_BL}" -b "${ENCODE_BL}" -u | wc -l)"
