# Prompt for loops/stripes methods documentation session

## Context

A collaborator (Challana) asked for the methods section of the Hi-C paper:

> Challana Tea [1:06 PM]
> Hey, finally getting around to methods section
> Can you write a section up for the loop analysis? It'll need details on the loop caller we used, what resolutions we used, our mariner & edgeR parameters, then how we did our loop annotation. Try to have all the software versions included b/c some journals are real sticklers about it
> Challana Tea [4:40 PM]
> Also for stripe detection and differential analysis

## Task

You are helping a bioinformatics student (Zara, she/her) build the reference materials needed to write a paper methods section. She will write the actual methods text herself in a separate session. You will not produce methods text directly. Your role is to explore the codebase, verify every detail against the actual scripts and configs, and produce two reference documents.

This is the same process that was used for a previous grant-writing task. Read these files to understand how that process worked and what failure modes to avoid:

@/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/original-prompts.md
@/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/grant-writing.txt

The grant-writing.txt is a full transcript of the writing session that used the reference materials. Read it to understand the AI failure modes that appeared during writing (negative parallelisms, dramatizing standard methods, substituting numbers for plain facts, shadowboxing against strawmen, "framework" as an AI-ism, etc.). The methods-instructions.md you produce will need to guard against all of those failure modes.

## This is NOT a one-shot task

First get a thorough understanding of both pipelines (loops and stripes) before producing any documents. Read the code, the configs, the READMEs, the CLAUDE.md files. In chat, have a conversation about what you find: what tools are used, what parameters, what the statistical models are, what the annotation pipeline does, what the stripe detection does. Do not start writing until we discuss.

## What to explore

Two pipelines live in this repo. Explore both.

**Chromatin loop analysis** (the `loops/` directory):
- @/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/loops/CLAUDE.md
- @/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/loops/README.md
- The scripts in `loops/scripts/` (the full pipeline: prep, extract, aggregate, edgeR, annotation, APA, downstream)
- The config files in `loops/config/`
- Loop caller: what tool calls loops, at what resolutions (5kb, 10kb, 25kb), with what parameters
- Mariner parameters: the R/Bioconductor mariner package for loop manipulation, merging, counting
- edgeR differential analysis: the model design, covariates, dispersion, offset, FDR threshold
- Loop annotation: how anchors are classified (ChIP-seq marks, genomic features, CTCF, etc.)
- Multi-resolution handling: how results across resolutions are compared or merged
- Software versions for every tool

**Stripe analysis** (the `stripes/` directory):
- @/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/stripes/stripenn/CLAUDE.md
- @/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/stripes/stripenn/README.md
- Stripenn: how stripes are called, at what resolution, with what parameters
- Differential stripe analysis: the statistical model, design, FDR threshold
- Stripe annotation
- The quagga pipeline in `stripes/quagga/` if relevant
- Software versions

For both pipelines, read the actual R scripts and SLURM submission scripts to get the real parameter values. Do not guess from documentation alone.

## Deliverables (produced AFTER discussion)

### 1. `loops-stats.md`

A comprehensive statistical methodology reference, same structure as the methylation `stats.md` at @/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/stats.md (read it for the format and level of detail). Every number must come from a verified script, config, or output file. Covers:

For loops:
- Loop caller name, version, parameters, resolutions
- Mariner package: version, how loops are merged/padded, pixel extraction parameters
- edgeR model: design matrix, covariates, dispersion estimation, offset, FDR threshold, the actual test used
- Multi-resolution comparison method
- Loop annotation: what marks, what annotation tool, how anchors are classified
- APA analysis parameters
- Sample structure (n per condition, replicates, .hic file sources)
- Key verified numbers (total loops called, significant differential, direction split, etc.)
- Software versions for every tool in the pipeline

For stripes:
- Stripenn: version, parameters, resolution
- Differential model: same level of detail as loops
- Key verified numbers
- Software versions

### 2. `methods-instructions.md`

A context handoff document for a future Claude session that will help Zara write the actual methods text. Same structure as @/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/grant-instructions.md (read it for the format). Covers:

- What the task is (write methods section for the Hi-C paper, covering loops and stripes)
- What files to read (with @ paths to force loading)
- The loops-stats.md as single source of truth
- Style rules: paper methods register (not grant), STE skill, anti-AI writing style guide
- AI failure modes to avoid (draw from grant-writing.txt observations AND the general list)
- What Challana specifically asked for (loop caller, resolutions, mariner & edgeR parameters, loop annotation, software versions, stripe detection, differential analysis)
- Constraints: never fabricate versions/flags/parameters, verify against actual files, Zara writes the text
- What Zara will do with the output (read it, absorb it, write the methods herself)

## Read these files

Read every file listed here before you respond:

@/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/chat.md
@/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/anti-ai-writing-style-guide.md
@/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/loops/CLAUDE.md
@/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/loops/README.md
@/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/stripes/stripenn/CLAUDE.md

Also invoke `/ste` at the start of the session to load the Simplified Technical English skill.

## Style rules

Write all prose in the register of a methods section: third person, past tense for completed analyses, formal but plain. Not the same register as a grant (which uses "we will" and testable predictions). A methods section states what was done and how.

Apply the STE skill: maximum 25 words per descriptive sentence, no semicolons, simple verb forms only, no "should"/"would"/"could"/"may"/"might."

Apply the anti-AI writing style guide at @/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/anti-ai-writing-style-guide.md. Hard bans: no AI vocabulary (delve, leverage, robust, crucial, seamless, landscape, etc.), no negative parallelisms, no meta-commentary, no em dashes (unless the paper itself uses them), no participle-tail analysis.

## AI failure modes to avoid

These are drawn from two prior sessions working on the methylation grant. Read @/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/grant-writing.txt for the full context of how each failure appeared.

**Dramatizing standard methods.** Gene-body aggregation for mCH was described with noise-floor numbers and "answering the reviewer question" framing when it is just what people do. Do not over-justify established practices. State what was done.

**Substituting numbers for plain facts.** AI has a tendency to throw in specific numbers (3x more sites, etc.) to make standard choices sound like benchmarking results. A methods section states the method. Numbers belong in results.

**Fix narratives.** Never frame a working method as a fix for a prior broken approach. State the method as it is.

**Negative parallelisms.** "It is not X, it is Y" and all variants. The single most reliable AI tell. Just state Y.

**Shadowboxing against strawmen.** "specified in advance rather than selected post hoc" is AI defending against imagined bad science. Nobody needs to be told you are not doing bad science.

**AI-ism vocabulary.** "Framework" (use "pipeline" or the tool name). "Cascade" (just describe the chain). "Leveraging" (use "using"). "Robust" (banned). Read the full banned list.

**Slop vocabulary contamination.** Your language shapes Zara's mental model of how to write the methods. Clean language is a requirement. If you use slop, she absorbs slop patterns, and those leak into the deliverable.

**Pontificating about pipeline history.** Do not spend time on how analytical decisions were arrived at, what was tried first, or what changed between iterations. State the method as it is now.

**Hedging.** A methods section states what was done. No "may," "could," "might."

## Constraints

- Do not produce methods text directly. Zara writes it with help from the handoff session.
- Verify every number against a file before stating it. Read the actual scripts and configs.
- Never fabricate software versions, CLI flags, path names, or parameter values. Read the source.
- Put both documents in the repo root (not in a docs/ subdirectory).
- Whenever you reference a file in the instructions document, use the @ prefix with the full path so the next session is forced to read it.

## What Zara will do with your output

She will take the two documents to a separate session. That session reads them, discusses the methods content with her, and helps her write the actual methods text. Your language in both documents directly shapes the final methods section through two layers of AI mediation. That is why clean, accurate, plain language matters even in the reference documents.
