#!/usr/bin/env python3
"""Generate ubiquitous genes list from DESeq2 baseMean.
Depends on Step 0A (needs mm10_genes.bed).

Bundled format: one gene name per line, no header.
"""

import pandas as pd
import sys

RNASEQ = "/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx"
GENES_BED = "/expanse/lustre/projects/csd940/zalibhai/abc/reference/mm10_genes.bed"
OUTPUT = "/expanse/lustre/projects/csd940/zalibhai/abc/reference/mm10_ubiquitous_genes.tsv"

# Read RNA-seq results (16,572 genes)
rna = pd.read_excel(RNASEQ)
print(f"RNA-seq entries: {len(rna)}")

# Column 'ensembl_gene_id' actually contains gene symbols (verified)
rna = rna.rename(columns={"ensembl_gene_id": "gene_symbol"})

# Read gene annotation (has header starting with #)
genes = pd.read_csv(GENES_BED, sep="\t", comment="#", header=None,
    names=["chr", "start", "end", "gene_symbol", "score", "strand", "ensembl_id", "gene_type"])
valid_genes = set(genes["gene_symbol"])
print(f"Genes in annotation: {len(valid_genes)}")

# Filter to genes present in our annotation
rna_matched = rna[rna["gene_symbol"].isin(valid_genes)].copy()
print(f"RNA-seq genes matching annotation: {len(rna_matched)}")

# Top 5% by baseMean = ubiquitously expressed
threshold = rna_matched["baseMean"].quantile(0.95)
ubiq = rna_matched[rna_matched["baseMean"] >= threshold]["gene_symbol"]
print(f"baseMean threshold (95th percentile): {threshold:.1f}")
print(f"Ubiquitous genes: {len(ubiq)}")

# Save: one gene per line, no header (matches bundled format)
ubiq.to_csv(OUTPUT, index=False, header=False)
print(f"Written to: {OUTPUT}")

# Spot-check known housekeeping genes
known_hk = ["Actb", "Gapdh", "Rpl13a", "Rps18", "Hsp90ab1"]
for g in known_hk:
    found = g in ubiq.values
    print(f"  {g}: {'FOUND' if found else 'MISSING'} in ubiquitous list")
