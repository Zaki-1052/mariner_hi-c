```python
$ cd /home/user/mariner_hi-c/biomodal/downstream/modality/outputs/run-5/outputs_CG/Results/gencode.vM25.mouse.cpg_islands.annotation/DMR_20260402_191006 && python3 <<'EOF'
import csv

mc_file = "DMR_mc_control__mutant_20260402_191006.bed"
hmc_file = "DMR_hmc_control__mutant_20260402_191006.bed"
out_file = "/home/user/mariner_hi-c/biomodal/downstream/modality/exports/cpg_islands/cpg_islands_CG_run-5_mc_hmc.tsv"

def load(path):
    rows = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            key = (r["Chromosome"], r["Start"], r["End"])
            rows[key] = r
    return rows

mc = load(mc_file)
hmc = load(hmc_file)

keys = sorted(set(mc) | set(hmc), key=lambda k: (k[0], int(k[1]), int(k[2])))

cols = ["Chromosome", "Start", "End", "Annotation",
        "mc_num_contexts", "mc_mean_coverage",
        "mc_mean_mod_control", "mc_mean_mod_mutant",
        "mc_mod_fold_change", "mc_mod_difference",
        "mc_test_statistic", "mc_pvalue", "mc_qvalue",
        "hmc_num_contexts", "hmc_mean_coverage",
        "hmc_mean_mod_control", "hmc_mean_mod_mutant",
        "hmc_mod_fold_change", "hmc_mod_difference",
        "hmc_test_statistic", "hmc_pvalue", "hmc_qvalue"]

with open(out_file, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(cols)
    for k in keys:
        m = mc.get(k, {})
        h = hmc.get(k, {})
        ref = m or h
        row = [
            ref["Chromosome"], ref["Start"], ref["End"], ref["Annotation"],
            m.get("num_contexts",""), m.get("mean_coverage",""),
            m.get("mean_mod_group_1",""), m.get("mean_mod_group_2",""),
            m.get("mod_fold_change",""), m.get("mod_difference",""),
            m.get("test_statistic",""), m.get("dmr_pvalue",""), m.get("dmr_qvalue",""),
            h.get("num_contexts",""), h.get("mean_coverage",""),
            h.get("mean_mod_group_1",""), h.get("mean_mod_group_2",""),
            h.get("mod_fold_change",""), h.get("mod_difference",""),
            h.get("test_statistic",""), h.get("dmr_pvalue",""), h.get("dmr_qvalue",""),
        ]
        w.writerow(row)

print(f"Wrote {len(keys)} CpG islands to {out_file}")
EOF


Wrote 8910 CpG islands to /home/user/mariner_hi-c/biomodal/downstream/modality/exports/cpg_islands/cpg_islands_CG_run-5_mc_hmc.tsv
```