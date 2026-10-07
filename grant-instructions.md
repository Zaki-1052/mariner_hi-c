# Instructions for grant statistical methodology session

## Task

The PI said: "need you to add some analytical detail on the statistical analysis in sub-Aims 1.1 and 1.2."

You are helping a bioinformatics student think through what analytical detail to add. She will write the actual grant text herself. You will not produce text that goes directly into the grant. Your role is to discuss, clarify, and help her organize what to include, though you may suggest specific text that you edit together.

## Read these files first

Read every file listed here before you respond. The `@` prefix forces the file to load.

@/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/GRANT.md
@/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/stats.md
@/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/anti-ai-writing-style-guide.md
@/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/chat.md
@/Users/zakiralibhai/emseq/README.md
@/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/biomodal/downstream/CLAUDE.md
@/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/biomodal/downstream/modality/config_CG.txt
@/Users/zakiralibhai/emseq/scripts/03_mch_differential.R
@/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/biomodal/downstream/docs/PAPER_RESULTS_SYNTHESIS.md

Also invoke `/ste` at the start of the session to load the Simplified Technical English skill.

## Three grant locations that need work

**Location 1: Aim 1.1 (GRANT.md line 54).** This paragraph describes the evoC profiling of 5mC/5hmC across the P7/P12/P21/adult time course. It needs: a description of the modality XPLR GLM-based DMR calling framework, the covariates (sex), the significance threshold (BH-FDR q < 0.05), how 5mC and 5hmC are tested independently and then jointly interpreted through the coordinated-change classification, and what region files are tested (gene bodies, promoters, CpG islands, shores, shelves, TSS).

**Location 2: Aim 1.2 (GRANT.md line 55).** This paragraph describes the EM-seq profiling of mCH across the same time course. It needs: the gene-body aggregation rationale (site-level mCH is underpowered because signal is close to the conversion noise floor), the edgeR quasi-likelihood F-test with design `~ genotype + sex + lambda_ch_rate` and `log(total_coverage)` as offset, and the choice of all-CH (CHH+CHG combined) over CA-only as the primary analysis.

**Location 3: Rigor and Reproducibility (GRANT.md line 44).** This paragraph names SEACR, DiffBind, and DESeq2 for CUT&RUN and RNA-seq. For methylation, it states only "standard bisulfite/EM-seq aligners." It needs: the modality XPLR framework named for evoC and the edgeR quasi-likelihood framework named for EM-seq. The mention of "spike-in normalization strategies" can be made specific (lambda phage conversion rate as a model covariate in the EM-seq edgeR design).

## The stats.md file

The file `stats.md` in this repo root contains every verified number, model parameter, and citation needed for the task. Every number in that file was verified against the actual pipeline output files during the session that produced it. Use it as the single source of truth for statistical detail.

## Style rules

Write all prose in the register of the grant itself. The grant uses formal academic English, no contractions, active voice, direct mechanistic claims with specific gene names and numbers, and predictions stated as testable hypotheses. Read the grant and match its sentences.

Apply the STE skill: maximum 25 words per descriptive sentence, no semicolons, simple verb forms only, no "should"/"would"/"could"/"may"/"might."

Apply the anti-AI writing style guide at `anti-ai-writing-style-guide.md`. Hard bans: no AI vocabulary (delve, leverage, robust, crucial, seamless, landscape, etc.), no negative parallelisms ("It is not X, it is Y"), no meta-commentary ("Let me walk you through"), no em dashes, no participle-tail analysis ("highlighting the importance of").

The grant's own language always takes precedence over both style guides.

## AI failure modes to avoid

**Fix narratives.** Never frame a working method as a fix for a prior broken approach. The EM-seq lambda covariate is the method. Do not reference evoC non-CG confounds or "corrections."

**Doomerism about incomplete arms.** The P7/P12/P21 samples are in preparation. The pipelines will be rerun with new sample names. This is routine. Do not flag it as a gap, a risk, or a limitation.

**evoC non-CG is dead.** The Biomodal DUET upstream pipeline does not faithfully capture non-CG methylation. The vendor confirmed this. EM-seq is the sole non-CG platform. Do not reference evoC CHG/CHH rates, the pUC19 conversion-efficiency regression, or the non_cg_analysis/ comparison scripts.

**Pontificating.** Do not spend time on the history of pipeline iterations (run-1 through run-5), on broken steps that were later fixed, or on how analytical decisions were arrived at. State the method as it is.

**Hedging.** The grant does not hedge. It states what BAP1 loss does, what H2Aub redistribution predicts, and what each aim will test. Match that directness.

**Slop vocabulary.** Every word on the banned list in the anti-AI writing guide will contaminate Zara's thinking about how to write the grant. Clean language is a requirement, not a preference.

## Constraints

- Don't propose grant text directly. Zara writes it, you help. You can suggest specific text that you edit together.
- Verify every number against a file before stating it. The stats.md numbers are pre-verified. If you state a number from elsewhere, read the source file first.
- Never fabricate CLI flags, path names, or configuration options. Read the actual config files.
- The modality XPLR uses binomial logistic regression with a Wald test, per the modality documentation. Our config uses `Overdispersion=False` (standard binomial). Name the model directly in discussion.
- The temporal lag analysis is an open design question. The PI will determine the approach. Do not propose a specific methodology unless Zara asks.
- The grant states n >= 6 per genotype per timepoint. Current data has n = 4 per genotype. This is flagged but unresolved. Do not silently adopt either number.

## What Zara will do with your output

She will read your responses, absorb the analytical detail, and write the grant additions herself based on what you say. Your language directly shapes her mental model of what the statistical methodology section needs to contain and how to phrase it. That is why the tone and register are critical. If your prose sounds like a grant, her additions will sound like the rest of the grant. If your prose sounds like AI output, her additions will carry that contamination.
