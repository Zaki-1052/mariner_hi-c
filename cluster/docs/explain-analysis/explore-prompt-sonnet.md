  Explore(Explore data file formats) Sonnet 4.6
  ⎿  Prompt:
       I'm planning Phase 4 of a Hi-C clustering analysis at /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/cluster/. I need
        to verify the actual on-disk format of all input data files for Phase 4 so my plan doesn't assume schemas that don't
       match reality.

       Please inspect each file with head (or wc -l + awk for big ones) and report exact column structure + sample values. Key
       facts I need:

       1. Loop cluster file —
       /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/cluster/bap1_late/cluster3/k-6/data/combined-clusters.txt
         - First 3 lines, column count, separator, header format
         - Confirm columns: GROUP, chr1, x1, x2, chr2, y1, y2, ctrl_merge, mut_merge
         - Total row count
         - Sorted? Group values present (clust1..clust6)?
       2. Loop metadata sidecar — /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/cluster/data/late_merged_loop_metadata.tsv
         - First 3 lines, column count, header
         - Confirm direction column values are exactly up_in_mutant / down_in_mutant / unchanged
         - Confirm coordinate cols match combined-clusters.txt for join
       3. DiffBind result files — three files in /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/peaks/diffbind/:
         - K27ac_diffbind_results_summit_appended_ap.txt
         - K27me3_diffbind_results_summit_appended_ap.txt
         - K119ub_diffbind_results_summit_appended_ap.txt
         - First 2 lines (header + 1 data row) of each
         - Confirm column names: Peak_Chr, Peak_Start, Peak_End, Fold, FDR (which positions/indices?)
         - What separator? Any header rows or skiprows needed?
       4. ChromHMM artifacts —
         - /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/cluster/bap1_late/chromHMM/learned_model/cerebellum_late_12_segmen
       ts.bed — first 3 lines, confirm 4-col format chr start end E{N}
         - /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/cluster/bap1_late/chromHMM/12state_rename_cerebellum.txt — full
       file (12 lines), confirm format E{N}<TAB>name
       5. mm10 gene annotation — /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/cluster/data/mm10_knownGene_pp.bed
         - First 3 lines, column count
         - Confirm 7-col gencode-style format
       6. Peak BEDs for loop classification (4.2) —
         - /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/peaks/CTCF.bed — first 2 lines, col count
         - /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/peaks/beds/H3K27acCerebellumLate2.bed — first 2 lines, col count
       7. BigWigs for ChIP signal at anchors (4.3) — verify these EXIST and report file size:
         - /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/peaks/bigwigs/macs2.narrow.aug18.dedup/index_13_ctrl_1_H3K27ac_S25
       _L001_aligned_reads.sorted.bw
         - /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/peaks/bigwigs/macs2.narrow.aug18.dedup/index_25_ctrl_1_H3K27me3_S3
       7_L001_aligned_reads.sorted.bw
         - /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/peaks/bigwigs/macs2.narrow.aug18.dedup/index_19_mut_1_H3K27ac_S46_
       L002_R2_001.fastq.gz_aligned_reads.sorted.bw
         - /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/peaks/bigwigs/macs2.narrow.aug18.dedup/index_23_mut_1_H3K27me3_S50
       _L002_R1_001_aligned_reads.sorted.bw
         - /Users/zakiralibhai/sdsc/bigwigs/H2AK119ubCtrl.bw and H2AK119ubMut.bw
         - /Users/zakiralibhai/sdsc/bigwigs/H3K27me1Ctrl.bw and H3K27me1Mut.bw
         - List the full directory /Users/zakiralibhai/sdsc/bigwigs/ to see if there are other useful files
       8. Existing output directory: ls /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/cluster/bap1_late/figures/ to confirm
        the subdir tree (annotation, ChIP_intersect, deeptools, deeptools_input) created in Phase 1 still exists.

       Report all this concisely under 600 words. Use exact file paths so I can copy them. Flag any anomalies (e.g., diffbind
       FDR column might be named FDR.X or padj instead of plain FDR — that would matter).
