# Relevant edgeR Documentation for Mariner Differential Loop Analysis

## Core Concepts

### 2.1 Terminology
- edgeR performs differential abundance analysis for pre-defined genomic features
- Features are referred to as "genes" in documentation but apply to any genomic feature (including chromatin loops)
- Analysis uses read counts (or in your case, aggregated Hi-C contact counts)
- Counts should NOT be transformed before input to edgeR

### 2.6 The DGEList Data Class
edgeR stores data in a `DGEList` object containing:
- `counts`: matrix of integer counts (rows = features/loops, columns = samples)
- `samples`: data.frame with sample information including:
  - `lib.size`: library size (column sums of counts)
  - `norm.factors`: normalization factors (default = 1)
  - `group`: factor identifying sample groups
- `genes`: optional data.frame with feature annotation (genomic coordinates)

Create DGEList:
```r
y <- DGEList(counts=x, group=group)
```

### 2.7 Filtering
Filter lowly expressed features before analysis:
```r
keep <- filterByExpr(y, group=group)
y <- y[keep, , keep.lib.sizes=FALSE]
```

`filterByExpr()`:
- Keeps features with worthwhile counts in minimum number of samples
- Automatically determines thresholds based on library sizes and design
- Filtering is independent of which sample belongs to which group (no bias)
- Recommended approach - avoids manual threshold setting

### 2.8 Normalization

#### 2.8.1 Normalization is Only Necessary for Sample-Specific Effects
- edgeR concerned with relative changes between conditions
- Technical factors that are unrelated to experimental conditions cancel out
- Normalization only needed for sample-specific effects

#### 2.8.2 Sequencing Depth
- Automatically adjusted via library sizes
- Part of basic modeling, no user intervention needed

#### 2.8.3 Effective Library Sizes (TMM Normalization)
`normLibSizes()` normalizes library sizes to minimize log-fold changes between samples for most genes:
```r
y <- normLibSizes(y)
```

- Uses TMM (Trimmed Mean of M-values) by default
- Recommended when majority of features are not differentially expressed
- Accounts for compositional biases
- Creates "effective library size" = original library size × scaling factor
- Normalization factors multiply to unity (geometric mean preserved)

#### 2.8.6 Model-Based Normalization, Not Transformation
- Normalization takes form of correction factors in statistical model
- Original read counts are NOT transformed
- Do NOT input RPKM, FPKM, or other transformed values
- Do NOT add artificial values to counts

### 2.9 Negative Binomial Models

#### 2.9.2 Biological Coefficient of Variation (BCV)
- BCV = coefficient of variation of true abundance between biological replicates
- Represents variation that would remain even with infinite sequencing depth
- Square root of dispersion parameter
- For well-controlled experiments: typical BCV ~0.4 for human, ~0.1 for model organisms

#### 2.9.3 Estimating BCVs
Three approaches:
1. **Common dispersion**: Assumes all features have same dispersion
2. **Trended dispersion**: Dispersion varies with expression level
3. **Tagwise dispersion**: Gene-specific dispersions with empirical Bayes moderation (RECOMMENDED)

### 2.12 The Quasi-Likelihood Pipeline (edgeR v4) - RECOMMENDED

#### 2.12.1 QL Generalization of NB Model
- Extends negative binomial with quasi-likelihood methods
- Accounts for gene-specific variability from biological AND technical sources
- Variance: `var(y) = σ²(μ + φμ²)` where:
  - φ = NB dispersion (global biological variability)
  - σ² = QL dispersion (gene-specific variability)

#### 2.12.2 Testing for DE
```r
fit <- glmQLFit(y, design, robust=TRUE)
qlf <- glmQLFTest(fit, coef=2)  # or contrast=...
topTags(qlf)
```

**Key advantages**:
- More rigorous error rate control than likelihood ratio tests
- Accounts for uncertainty in dispersion estimation
- Uses F-tests instead of chi-square approximations
- `robust=TRUE` strongly recommended - identifies outliers

**Note**: In edgeR v4, `estimateDisp()` step is now optional before `glmQLFit()`

## Experimental Design for Two Groups

### 3.2.1 Introduction
Simplest design: comparing conditions based on independent biological replicates.

### 3.2.3 GLM Approach
For two groups (control vs mutant):

**Design matrix approach 1** (coefficient per group):
```r
design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- levels(y$samples$group)
```

Then test with contrast:
```r
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, contrast=c(-1,1))  # mutant - control
```

**Design matrix approach 2** (intercept + treatment effect):
```r
design <- model.matrix(~group, data=y$samples)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef=2)  # tests treatment effect
```

## No Biological Replicates (Section 2.13)

**Critical for your analysis**: Your data has n=1 per condition (merged replicates).

Options when no replicates available:

1. **Descriptive analysis only** (safest)
   - MDS plots, fold changes
   - No significance testing

2. **Pick reasonable dispersion** based on similar data
   - Typical BCV values: 0.4 (human), 0.1 (model organisms), 0.01 (technical replicates)
   ```r
   bcv <- 0.2
   et <- exactTest(y, dispersion=bcv^2)
   ```
   - Results very sensitive to chosen dispersion
   - Less controlled datasets may have much larger dispersions

3. **Estimate from reduced model** (if multiple factors)
   - Remove least important factors to create residual df
   - Only works if DE features are small proportion

4. **Estimate from control features** (if available)
   - Use housekeeping/non-responding features
   - Requires reasonably large number (dozens to hundreds)

**Important**: None of these are satisfactory alternatives to biological replication. Results should be interpreted cautiously.

## Key Functions

### Creating DGEList
```r
y <- DGEList(counts=counts_matrix, 
             group=factor(c("ctrl","mut")),
             genes=annotation_df)
```

### Filtering
```r
keep <- filterByExpr(y, group=y$samples$group)
y <- y[keep, , keep.lib.sizes=FALSE]
```

### Normalization
```r
y <- normLibSizes(y)  # TMM normalization
```

### Dispersion Estimation (optional in v4)
```r
y <- estimateDisp(y, design, robust=TRUE)
plotBCV(y)  # visualize dispersions
```

### QL Pipeline (RECOMMENDED)
```r
fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)  # visualize QL dispersions

qlf <- glmQLFTest(fit, coef=2)  # or contrast=...
topTags(qlf)
summary(decideTests(qlf))
```

### Results Extraction
```r
# Top differential features
topTags(qlf, n=20)

# Summary of DE features
summary(decideTests(qlf))

# All results
results <- topTags(qlf, n=Inf)$table
```

### Visualization
```r
# MD plot (log-fold-change vs average abundance)
plotMD(qlf)
abline(h=c(-1,1), col="blue")  # 2-fold change lines

# MDS plot (sample relationships)
plotMDS(y, col=c(1,2)[y$samples$group])
```

## Output Interpretation

### topTags() columns
- `logFC`: log2 fold-change (positive = higher in second group)
- `logCPM`: average log2 counts-per-million
- `F`: F-statistic (for QL F-test)
- `PValue`: p-value
- `FDR`: false discovery rate (Benjamini-Hochberg adjusted p-value)

### decideTests() summary
Shows number of features:
- `Down`: significantly decreased (FDR < 0.05 by default)
- `NotSig`: not significant
- `Up`: significantly increased

## Important Notes

1. **Input Requirements**:
   - Integer counts (or aggregated counts as in your case)
   - NOT normalized/transformed values
   - NOT RPKM/FPKM/TPM

2. **Library Sizes**:
   - Automatically calculated as column sums
   - Can be manually set if needed (as you did for equal ctrl/mut pairs)

3. **Normalization**:
   - TMM assumes majority of features are NOT differential
   - If this assumption violated, results may be biased

4. **Multiple Testing**:
   - FDR < 0.05 standard threshold
   - Controls false discovery rate (expected proportion of false positives)

5. **No Replicates**:
   - Severely limits statistical power
   - Cannot reliably estimate biological variation
   - Results should be validated experimentally

6. **Robust Options**:
   - `robust=TRUE` in `estimateDisp()` and `glmQLFit()` strongly recommended
   - Identifies and downweights outlier features

## Relevant Case Studies

### 4.3 Yoruba HapMap (Single Samples)
Shows analysis with multiple samples but demonstrates:
- Quality assessment via correlation
- MDS plotting
- Handling large numbers of features
- Gene set testing approaches

### 4.4 Mouse Mammary Gland
Demonstrates:
- Complete QL pipeline workflow
- Design matrix setup
- Multiple contrasts
- Visualization approaches
- Filtering strategy

## Software Versions
edgeR v4 (latest) includes:
- Optional `estimateDisp()` step
- Improved QL dispersion estimation
- Better handling of small counts
- Enhanced computational efficiency

---


# Additional edgeR Documentation from Reference Manual

## Data Input and Manipulation

### DGEList Construction and Manipulation

#### DGEList() - Creating the Core Object
```r
DGEList(counts, lib.size=NULL, norm.factors=NULL, samples=NULL, 
        group=NULL, genes=NULL, remove.zeros=FALSE)
```

**Key parameters**:
- `counts`: Matrix or data.frame of counts (rows=features, columns=samples)
- `lib.size`: Defaults to `colSums(counts)` if not specified
- `norm.factors`: Normalization factors (default=1 for all samples)
- `samples`: Data.frame with sample information (will be created if NULL)
- `group`: Experimental group for each sample
- `genes`: Data.frame with feature annotation
- `remove.zeros`: Logical, whether to remove rows with 0 total count

**For data.frame input**:
- `annotation.columns`: Specify which columns contain gene IDs vs counts
- Function will auto-detect non-numeric columns as annotation

**Returns**: DGEList object with components:
- `$counts`: Count matrix
- `$samples`: Sample information including lib.size, norm.factors, group
- `$genes`: Feature annotation (if provided)

### Subsetting DGEList Objects
```r
y[i, j, keep.lib.sizes=TRUE]
```

**Behavior**:
- `i`: Subset features (rows)
- `j`: Subset samples (columns)
- `keep.lib.sizes=TRUE`: Preserves original library sizes (recommended for filtering)
- `keep.lib.sizes=FALSE`: Recalculates library sizes from remaining features

**Important**: When filtering features, use `keep.lib.sizes=FALSE` to update library sizes.

### Combining DGEList Objects
```r
cbind(dge1, dge2, ...)  # Combine samples (same features)
rbind(dge1, dge2, ...)  # Combine features (same samples)
```

## Normalization Functions

### normLibSizes() - TMM and Other Methods
```r
normLibSizes(y, method=c("TMM","TMMwsp","RLE","upperquartile","none"),
             refColumn=NULL, logratioTrim=0.3, sumTrim=0.05, 
             doWeighting=TRUE, Acutoff=-1e10, p=0.75)
```

**Methods**:
1. **TMM** (Trimmed Mean of M-values) - DEFAULT, RECOMMENDED
   - Trims extreme log-ratios and A-values
   - `logratioTrim`: Fraction to trim from M-values (default 0.3)
   - `sumTrim`: Fraction to trim from A-values (default 0.05)
   - `doWeighting`: Use precision weights (default TRUE)
   - `refColumn`: Reference sample (auto-selected if NULL)

2. **TMMwsp** (TMM with singleton pairing)
   - More stable with high proportion of zeros
   - Reuses positive counts from zero-in-one-library genes
   - Better for sparse data

3. **RLE** (Relative Log Expression)
   - Median ratio to geometric mean across samples
   - Similar to DESeq2's default

4. **upperquartile**
   - Uses 75th percentile (or quantile specified by `p`)
   - Simple but less robust than TMM

5. **none**
   - Sets all norm.factors to 1
   - Use when no normalization needed

**Output**: Updates `y$samples$norm.factors`

**Effective library size** = `lib.size * norm.factors`

### getNormLibSizes() - Extract Effective Library Sizes
```r
getNormLibSizes(y, log=FALSE)
```
Returns `lib.size * norm.factors` for each sample.

## Filtering Functions

### filterByExpr() - Automatic Filtering
```r
filterByExpr(y, design=NULL, group=NULL, lib.size=NULL,
             min.count=10, min.total.count=15, large.n=10, 
             min.prop=0.7)
```

**Logic**:
- Keeps genes with `>= min.count` reads in a worthwhile number of samples
- "Worthwhile number" = smallest group size (or min inverse leverage from design)
- For large groups (>10 samples), requires expression in at least `min.prop` of samples
- Also requires `>= min.total.count` total reads across all samples

**Parameters**:
- `min.count`: Minimum count threshold (default 10)
- `min.total.count`: Minimum total across samples (default 15)
- `large.n`: What constitutes a "large" group (default 10)
- `min.prop`: Minimum proportion for large groups (default 0.7)

**Returns**: Logical vector indicating which features to keep

**Usage**:
```r
keep <- filterByExpr(y, group=y$samples$group)
y <- y[keep, , keep.lib.sizes=FALSE]
```

**Important**: This is the RECOMMENDED filtering approach. Avoids manual threshold setting.

## Dispersion Estimation

### estimateDisp() - All-in-One Dispersion Estimation
```r
estimateDisp(y, design=NULL, prior.df=NULL, trend.method="locfit",
             tagwise=TRUE, span=NULL, min.row.sum=5, 
             grid.length=21, grid.range=c(-10,10),
             robust=FALSE, winsor.tail.p=c(0.05,0.1))
```

**What it does**:
- Estimates common, trended, AND tagwise dispersions in one call
- Recommended over calling individual estimation functions separately

**Key parameters**:
- `trend.method`: How to estimate trend ("locfit", "none", "movingave", "loess")
- `tagwise=TRUE`: Estimate gene-specific dispersions (RECOMMENDED)
- `robust=TRUE`: Robustify against hypervariable genes (RECOMMENDED)
- `prior.df`: Prior degrees of freedom for EB shrinkage (auto-calculated if NULL)

**Returns**: DGEList with added components:
- `$common.dispersion`: Single dispersion value
- `$trended.dispersion`: Vector of trended dispersions
- `$tagwise.dispersion`: Vector of gene-specific dispersions
- `$prior.df`: Prior degrees of freedom used
- `$span`: Smoothing span used

**Note**: In edgeR v4, this step is OPTIONAL before `glmQLFit()`. The QL fit can estimate dispersions internally.

### Individual Dispersion Functions

#### estimateCommonDisp()
```r
estimateCommonDisp(y, tol=1e-06, rowsum.filter=5, verbose=FALSE)
```
Estimates single dispersion for all genes.

#### estimateTrendedDisp()
```r
estimateTrendedDisp(y, method="bin.spline", df=5, span=2/3)
```
Estimates dispersion trend with abundance.

#### estimateTagwiseDisp()
```r
estimateTagwiseDisp(y, prior.df=10, trend="movingave", span=NULL,
                    method="grid", grid.length=11, grid.range=c(-6,6))
```
Estimates gene-specific dispersions with empirical Bayes shrinkage.

**Parameters**:
- `prior.df`: Weight given to common/trended dispersion (higher = more shrinkage)
- `trend`: Use trended dispersion as prior ("movingave", "loess", "none")
- `method`: "grid" (interpolation) or "optimize" (direct optimization)

## GLM Functions

### glmQLFit() - Quasi-Likelihood GLM Fitting
```r
glmQLFit(y, design=NULL, dispersion=NULL, abundance.trend=TRUE,
         robust=FALSE, winsor.tail.p=c(0.05,0.1),
         legacy=FALSE, top.proportion=NULL)
```

**Key features** (edgeR v4):
- Can estimate NB dispersion internally if `dispersion=NULL`
- Uses adjusted deviances for better small-count handling (`legacy=FALSE`)
- Estimates quasi-dispersions accounting for gene-specific variability

**Parameters**:
- `design`: Design matrix (defaults to `model.matrix(~group)`)
- `dispersion`: NB dispersions (estimated internally if NULL)
- `abundance.trend=TRUE`: Allow QL dispersion to vary with abundance (RECOMMENDED)
- `robust=TRUE`: Robustify against outliers (STRONGLY RECOMMENDED)
- `legacy=FALSE`: Use new v4 method with adjusted deviances (RECOMMENDED)
- `top.proportion`: Proportion of top genes for initial dispersion estimate

**Returns**: DGEGLM object with components:
- `$coefficients`: Fitted coefficients (log scale)
- `$fitted.values`: Fitted values
- `$deviance`: Residual deviances
- `$df.residual.adj`: Adjusted residual df (for legacy=FALSE)
- `$df.prior`: Prior df for QL dispersions
- `$s2.prior`: Prior QL dispersion values
- `$s2.post`: Posterior QL dispersion estimates

### glmQLFTest() - Quasi-Likelihood F-Test
```r
glmQLFTest(glmfit, coef=ncol(glmfit$design), contrast=NULL,
           poisson.bound=TRUE)
```

**Parameters**:
- `coef`: Which coefficient(s) to test (column number or name)
- `contrast`: Numeric contrast vector (takes precedence over coef)
- `poisson.bound=TRUE`: Ensures p-values not less than Poisson LRT

**Returns**: DGELRT object with:
- `$table`: Results table with logFC, logCPM, F, PValue
- `$df.total`: Denominator df for F-test
- `$comparison`: Description of contrast tested

**Usage**:
```r
fit <- glmQLFit(y, design, robust=TRUE)
qlf <- glmQLFTest(fit, coef=2)  # Test coefficient 2
# OR
qlf <- glmQLFTest(fit, contrast=c(-1,1))  # Test contrast
```

### glmFit() - Standard GLM Fitting (Alternative to QL)
```r
glmFit(y, design=NULL, dispersion=NULL, prior.count=0.125,
       start=NULL)
```

**Simpler than QL**:
- Fits NB GLM without quasi-likelihood extension
- Faster but less rigorous error rate control
- Use when speed critical or QL assumptions questionable

**Returns**: DGEGLM object similar to glmQLFit but without QL components

### glmLRT() - Likelihood Ratio Test
```r
glmLRT(glmfit, coef=ncol(glmfit$design), contrast=NULL)
```

**Alternative to glmQLFTest**:
- Uses likelihood ratio test instead of F-test
- Less conservative than QL F-test
- Not recommended for routine use (QL F-test preferred)

## Results Functions

### topTags() - Extract Top Results
```r
topTags(object, n=10, adjust.method="BH", sort.by="PValue", p.value=1)
```

**Parameters**:
- `n`: Number of top features to return (use `n=Inf` for all)
- `adjust.method`: Multiple testing correction ("BH", "BY", "holm", "bonferroni", "none")
- `sort.by`: Sort by "PValue", "logFC", or "none"
- `p.value`: Only return features with adjusted p-value ≤ this threshold

**Returns**: TopTags object containing:
- `$table`: Data.frame with results
  - `logFC`: Log2 fold-change
  - `logCPM`: Average log2 counts-per-million
  - `F` or `LR`: Test statistic
  - `PValue`: Raw p-value
  - `FDR`: Adjusted p-value (if adjust.method="BH")
- `$adjust.method`: Method used
- `$comparison`: Description of test
- `$test`: Type of test performed

**Usage**:
```r
top <- topTags(qlf, n=20)  # Top 20
all_results <- topTags(qlf, n=Inf)$table  # All results as data.frame
```

### decideTests() - Classify Results
```r
decideTests(object, adjust.method="BH", p.value=0.05, lfc=0)
```

**Parameters**:
- `adjust.method`: Multiple testing correction
- `p.value`: Significance threshold (default 0.05)
- `lfc`: Log-fold-change threshold (default 0, NOT RECOMMENDED to change)

**Returns**: TestResults object (matrix with -1/0/1 values)
- `-1`: Significantly down-regulated
- `0`: Not significant
- `1`: Significantly up-regulated

**Usage**:
```r
summary(decideTests(qlf))  # Summary counts
is.de <- decideTests(qlf)
```

**Important**: Using `lfc` threshold not recommended. Use `glmTreat()` instead for fold-change testing.

### glmTreat() - Test Against Fold-Change Threshold
```r
glmTreat(glmfit, coef=ncol(glmfit$design), contrast=NULL, 
         lfc=log2(1.2), null="interval")
```

**Purpose**: Test if fold-change is significantly GREATER than threshold (not just different from zero)

**Parameters**:
- `lfc`: Log2 fold-change threshold (e.g., log2(1.5) for 1.5-fold)
- `null`: "interval" (recommended) or "worst.case"

**Returns**: DGELRT object similar to glmQLFTest

**Usage**:
```r
tr <- glmTreat(fit, coef=2, lfc=log2(1.5))
topTags(tr)
```

**Important**: This is the CORRECT way to incorporate fold-change thresholds, not using `lfc` in `decideTests()`.

## Visualization Functions

### plotMD() - Mean-Difference Plot
```r
plotMD(object, column=1, xlab="Average log CPM", 
       ylab="log-fold-change", status=NULL,
       main=NULL, ...)
```

**For DGEList**: Plots one sample vs average of others
**For DGELRT/DGEGLM**: Plots logFC vs logCPM from test results

**Parameters**:
- `column`: Which sample/coefficient to plot
- `status`: Vector indicating which points to highlight
- Color points by significance:
```r
plotMD(qlf, status=decideTests(qlf))
```

### plotMDS() - Multidimensional Scaling Plot
```r
plotMDS(y, top=500, labels=NULL, pch=NULL, cex=1,
        dim.plot=c(1,2), gene.selection="pairwise",
        method="logFC", prior.count=2)
```

**Purpose**: Visualize sample relationships in 2D

**Parameters**:
- `top`: Number of most variable features to use (default 500)
- `gene.selection`: "pairwise" (default) or "common"
- `method`: "logFC" (default) or "bcv"
- `dim.plot`: Which dimensions to plot (e.g., c(1,2) or c(2,3))

**Usage**:
```r
plotMDS(y, col=c("blue","red")[y$samples$group])
```

### plotBCV() - Plot Biological Coefficient of Variation
```r
plotBCV(y, xlab="Average log CPM", ylab="Biological coefficient of variation")
```

**Shows**:
- Common dispersion (horizontal line)
- Trended dispersion (blue curve)
- Tagwise dispersions (black points)

**Usage**: After `estimateDisp()`
```r
y <- estimateDisp(y, design)
plotBCV(y)
```

### plotQLDisp() - Plot QL Dispersions
```r
plotQLDisp(glmfit, xlab="Average Log2 CPM", 
           ylab="Quarter-Root Mean Deviance")
```

**Shows**:
- Raw QL dispersions (black points)
- Squeezed QL dispersions (red points)
- QL dispersion trend (blue line)

**Usage**: After `glmQLFit()`
```r
fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)
```

## Utility Functions

### cpm() - Counts Per Million
```r
cpm(y, normalized.lib.sizes=TRUE, log=FALSE, prior.count=2)
```

**Parameters**:
- `normalized.lib.sizes=TRUE`: Use effective library sizes (lib.size * norm.factors)
- `log=FALSE`: Return linear scale CPM
- `log=TRUE`: Return log2(CPM + prior.count)
- `prior.count`: Offset added before log transformation

**Returns**: Matrix of CPM values (same dimensions as counts)

**Usage**:
```r
log_cpm <- cpm(y, log=TRUE)
```

### aveLogCPM() - Average Log2 CPM
```r
aveLogCPM(y, normalized.lib.sizes=TRUE, prior.count=2, 
          dispersion=NULL)
```

**Returns**: Vector of average log2-CPM across samples for each gene

**More sophisticated than**:
```r
rowMeans(cpm(y, log=TRUE))  # Simple average
```

**Difference**: `aveLogCPM()` weights larger libraries more heavily and accounts for dispersion.

### addPriorCount() - Add Prior Count for Shrinkage
```r
addPriorCount(y, lib.size=NULL, offset=NULL, prior.count=1)
```

**Purpose**: Add library-size-adjusted prior count to shrink log-fold-changes

**Returns**: List with:
- `$y`: Matrix with added prior counts
- `$offset`: Adjusted offsets

**Note**: This is done automatically in `glmFit()` and `predFC()` with `prior.count` parameter.

## Working Without Replicates

### exactTest() with Fixed Dispersion
```r
exactTest(y, pair=c("control","treatment"), dispersion=0.04)
```

**For two-group comparison without replicates**:
```r
bcv <- 0.2  # Choose based on similar experiments
y <- DGEList(counts=counts, group=group)
y <- normLibSizes(y)
et <- exactTest(y, dispersion=bcv^2)
topTags(et)
```

**Important**: Results VERY sensitive to chosen dispersion. This is NOT a substitute for replication.

### Reasonable BCV Values
- **Technical replicates**: BCV ≈ 0.01
- **Model organisms (controlled)**: BCV ≈ 0.1
- **Human samples (controlled)**: BCV ≈ 0.4
- **Poorly controlled experiments**: BCV may be much larger

**Your case**: With merged biological replicates, reasonable to use BCV ≈ 0.1-0.2 for mouse data.

## Important Technical Notes

### Count Data Requirements
1. **Must be integer or integer-like counts**
   - Your aggregated Hi-C counts are acceptable
   - NOT normalized values (RPKM, FPKM, TPM)
   - NOT log-transformed values

2. **Library sizes**
   - Automatically calculated as column sums
   - Can be manually set if needed:
   ```r
   y$samples$lib.size <- c(1e7, 1e7)  # Equal sizes
   ```

3. **Zero counts**
   - Handled automatically in statistical models
   - Don't add artificial values (0.5, 1, etc.)
   - Filter out features with all zeros using `filterByExpr()`

### Design Matrix Setup

**Two-group comparison** (your case):

**Approach 1** - Group means parameterization:
```r
design <- model.matrix(~0 + group, data=y$samples)
colnames(design) <- levels(y$samples$group)
# Test with contrast
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, contrast=c(-1, 1))  # mut - ctrl
```

**Approach 2** - Treatment effect parameterization:
```r
design <- model.matrix(~group, data=y$samples)
# Test coefficient 2
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef=2)  # Tests treatment effect
```

Both give identical results. Approach 2 is simpler for two-group comparisons.

### Interpretation of Results

**logFC (log fold-change)**:
- Positive: Higher in second group (or treatment)
- Negative: Lower in second group (or treatment)
- log2 scale: logFC=1 means 2-fold change, logFC=2 means 4-fold change

**logCPM**:
- Average expression level across all samples
- Higher values = more highly expressed features
- Used for filtering and visualization

**F-statistic** (from glmQLFTest):
- Larger values = stronger evidence for differential expression
- Converted to p-value using F-distribution

**PValue**:
- Raw p-value from statistical test
- Not adjusted for multiple testing

**FDR** (False Discovery Rate):
- Adjusted p-value controlling expected proportion of false positives
- Standard threshold: FDR < 0.05
- More stringent: FDR < 0.01

### Common Pitfalls to Avoid

1. **Don't transform counts before input**
   - edgeR expects raw counts
   - Normalization handled internally

2. **Don't use logFC cutoffs in decideTests()**
   - Use `glmTreat()` instead for proper fold-change testing

3. **Don't skip filtering**
   - Low-count features reduce power
   - Use `filterByExpr()` before analysis

4. **Don't ignore warnings about dispersions**
   - If dispersions very large (>1), check data quality
   - May indicate outliers or batch effects

5. **Don't use LRT when QL available**
   - QL F-test (`glmQLFTest`) has better error rate control
   - LRT (`glmLRT`) too liberal for small samples

### Recommended Workflow Summary

```r
# 1. Create DGEList
y <- DGEList(counts=counts, group=group, genes=annotation)

# 2. Filter low counts
keep <- filterByExpr(y, group=group)
y <- y[keep, , keep.lib.sizes=FALSE]

# 3. Normalize
y <- normLibSizes(y)

# 4. Design matrix
design <- model.matrix(~group, data=y$samples)

# 5. QL fit (robust=TRUE strongly recommended)
fit <- glmQLFit(y, design, robust=TRUE)

# 6. QL test
qlf <- glmQLFTest(fit, coef=2)

# 7. Results
topTags(qlf, n=20)
summary(decideTests(qlf))

# 8. Visualize
plotMD(qlf, status=decideTests(qlf))
plotQLDisp(fit)
```

This is the RECOMMENDED pipeline for edgeR v4 with your type of data.