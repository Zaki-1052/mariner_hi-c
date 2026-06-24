## Section 09: Summary Statistics

**Key numbers:**
- Genes tested = 20,969; 5mC sig = 10,775 (51.4%, 7,513 hyper / 3,262 hypo, mean +1.72%) (source: `summary_statistics.txt`)
- 5hmC sig = 11,484 (54.8%, 1,963 up / 9,521 down, mean -1.66%) (source: `summary_statistics.txt`)
- Co-significant = 8,371; coordinated mC↑/hmC↓ = 6,589 (78.7%) (source: `summary_statistics.txt`)
- Top coordinated genes: Tmem238, Syt1, Gclm, Sap30, Prxl2b (source: `summary_statistics.txt`)

**What this shows:** Section 09 is a machine-generated aggregation of all core DMR statistics (Sections 03/05/07). Its numbers are the canonical run-5 (8-sample) values, populated live from the run-5 DMR objects. Text-only output; this section produces no image (the report is `tables/summary_statistics.txt`).

**Figures:**
- (none — text report only; folder is a placeholder)

> Data-quality note: the `Samples: 4 (2 Control, 2 Mutant)` line in `summary_statistics.txt` is a hardcoded template string in `section_09_summary.R` (line 50) — a cosmetic bug. The numeric content is genuine run-5 8-sample data (independently verified against the run-5 BED files). Recommend correcting line 50 to "Samples: 8 (4 Control, 4 Mutant)".
