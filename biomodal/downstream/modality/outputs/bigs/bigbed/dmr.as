table dmr
"Differentially Methylated Region from modality XPLR"
    (
    string chrom;        "Chromosome"
    uint chromStart;     "Start position"
    uint chromEnd;       "End position"
    string name;         "Annotation"
    uint score;          "Score (scaled -log10 qvalue, 0-1000)"
    char[1] strand;      "Strand"
    float modDiff;       "Methylation difference (mutant - control)"
    float qvalue;        "Adjusted p-value (FDR)"
    float foldChange;    "Fold change"
    )
