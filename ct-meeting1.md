- juicer - mapq score (filtering
- pairix)
- convert containers to singularity #todo
- envs cooler.yml
- dont remove duplicates in cnr
  distance decay not working
- #todo : insulations at tad boundaries - cntrl mut 

```bash
#!/bin/bash
#SBATCH --job-name="nfhic"
#SBATCH --output="nfhic_nolig.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=200G
#SBATCH --account=csd940
#SBATCH --export=ALL
#SBATCH -t 48:00:00

##
module reset
module load singularitypro/4.1.2
source activate env_nf

export NXF_SINGULARITY_CACHEDIR=/expanse/lustre/scratch/$USER/temp_project/.singularity
export NXF_WORK=/expanse/lustre/scratch/$USER/temp_project/job_$SLURM_JOBID
export NXF_SINGULARITY_RUNOPTIONS="--bind /scratch,/expanse"

##Specify directories
SAMPLE_NAME=$1
BASE_DIR='/expanse/lustre/projects/csd940/ctea/homer/trimmed'
FASTQ_INPUT="${BASE_DIR}/trimmed_${SAMPLE_NAME}/*${SAMPLE_NAME}*_{1,2}.fastq.gz"
B2_INPUT="${BASE_DIR}/trimmed_batch2_${SAMPLE_NAME}/*${SAMPLE_NAME}*_{1,2}.fastq.gz"

# Find read 1 and read 2 fastq files
FASTQ_R1=$(ls ${BASE_DIR}/trimmed_${SAMPLE_NAME}/*${SAMPLE_NAME}*_R1.fastq.gz)
FASTQ_R2=$(ls ${BASE_DIR}/trimmed_${SAMPLE_NAME}/*${SAMPLE_NAME}*_R2.fastq.gz)
B2_R1=$(ls ${BASE_DIR}/trimmed_batch2_${SAMPLE_NAME}/*${SAMPLE_NAME}*_R1.fastq.gz)
B2_R2=$(ls ${BASE_DIR}/trimmed_batch2_${SAMPLE_NAME}/*${SAMPLE_NAME}*_R2.fastq.gz)

# Output CSV name
CSV_FILE="${SAMPLE_NAME}_sample_sheet.csv"

# Write to CSV
echo "sample,fastq_1,fastq_2" > "$CSV_FILE"
echo "${SAMPLE_NAME},${FASTQ_R1},${FASTQ_R2}" >> "$CSV_FILE"
echo "${SAMPLE_NAME},${B2_R1},${B2_R2}" >> "$CSV_FILE"

echo "Sample sheet written to: $CSV_FILE"

##

nextflow run nf-core/hic -profile singularity --input "${CSV_FILE}" \
	-r 2.1.0 \
	--outdir "/expanse/lustre/projects/csd940/ctea/nf-hic/trimmed_"${SAMPLE_NAME}"" \
	--genome 'mm10' \
	--digestion 'arima' \
	--bin_size '1000000,500000,250000,100000,40000,25000,10000,5000,2000,1000' \
	--res_dist_decay '250000,100000,50000' \ 
	--res_compartments "1000000,250000,100000" \
	--tads_caller 'insulation,hicexplorer' \
	--res_tads '40000,20000,10000' \
	--hicpro_maps \
	--save_raw_maps \
	-work-dir $NXF_WORK

echo "NF-hic pipeline completed."


```



| bin | count |
| --- | ----- |
| 1   | 2     |
| 2   | 3     |
| 3   | 2     |
| 1-3 | 7     |
| 4-6 | 2     |


- multi-qc 
- 300 mill valid pairs. -- 10 kb 
- 500 to 600: 5
- micro c  - cant confirm peis
- capture hi-c:
	- probes for promoters
	- 
- resolution limits:
	- max 5kb
- singletons: paired ends looped around to one site
- begin 90%, pairing 60%
- 30% duplicate reads from fastqc
- mcool raw - multires
  
- nf-hic
- trim with fastp
- valid_pairs
	- call contact maps
	- pairix failes .pairs
		- juicer
- juicer-pre
	- take pairs
	- same as fastq
- juicer: 1.6 vs 2 vs v9
	- using 2
	- - hictk
- juicebox visualization 
- genova isnst replicate aware
	- need to merge control vs mutant
- no nested tads
- multihicompare - homer- compartments, tads, not loops
	- define consesus tad set
- classify tad changes- if theyve rearranged


- set up analysis:
	- https://www.science.org/doi/10.1126/sciadv.adv2283
- increase or decrease boundaries 
	- define consensus peak/boundary set
- compartment pc - principal component analysis
- give homer a peak file
- multi-hi-compare:
	- contacts including loops
	- more computationally intensive
	- trans interactions no diff
	- differential analysis tool
- tallied signal is ap - signal over region for peak
- manhattan plots
	- gwas studies
	- each dot is a bin
	- at that bin is there a significant contact
- chrom 10
	- syt1/nav3 locus
	- very high insulation score
	- most impacted region
- loop diff analysis - stuck
	- contacts
	- diff stats
	- cLoops?? -- no
- https://github.com/dozmorovlab/multiHiCcompare
	- best at whats a differential contact
- explorer not replicate aware
- hiccups best for loop calling
- merged replicates and control.mutant
- 3-5 billion reads is best vs 600 mill
- t

```bash
#!/bin/bash
#SBATCH --job-name="fastp"
#SBATCH --output="fastp.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=50G
#SBATCH --account=csd940
#SBATCH --export=ALL
#SBATCH -t 8:00:0

module reset
source activate homer1

#Set variables
SAMPLE_NAME=$1
FASTQ_DIR=$2
OUTPUT_DIR="/expanse/lustre/projects/csd940/ctea/homer/trimmed_batch2_${SAMPLE_NAME}"
THREADS=16

# Create output directories
mkdir -p ${OUTPUT_DIR}/logs
mkdir -p ${OUTPUT_DIR}/temp

# Check if required arguments are provided
if [ -z "${SAMPLE_NAME}" ] || [ -z "${FASTQ_DIR}" ]; then
    echo "Usage: sbatch script.sh SAMPLE_NAME FASTQ_DIRECTORY"
    echo "Example: sbatch script.sh sample1 /path/to/fastq/files"
    exit 1
    fi

    # Check if FASTQ directory exists
    if [ ! -d "${FASTQ_DIR}" ]; then
        echo "ERROR: FASTQ directory ${FASTQ_DIR} does not exist"
	exit 1
    fi

# Function to find FASTQ files
find_fastq_files() {
    local sample=$1
    local dir=$2
   # Common patterns for R1/R2 files
    local patterns=(
    "*${sample}*R1_001.fastq.gz *${sample}*R2_001.fastq.gz"
    "${sample}_R1*.fastq.gz ${sample}_R2*.fastq.gz"
    "${sample}_1*.fastq.gz ${sample}_2*.fastq.gz"
    "${sample}*_R1*.fastq.gz ${sample}*_R2*.fastq.gz"
    "${sample}*_1*.fastq.gz ${sample}*_2*.fastq.gz"
    "${sample}.R1*.fastq.gz ${sample}.R2*.fastq.gz"
    "${sample}.1*.fastq.gz ${sample}.2*.fastq.gz")
								        
for pattern in "${patterns[@]}"; do
       set -- $pattern
       local r1_pattern=$1
       local r2_pattern=$2
	       
       local r1_files=($(find ${dir} -name "${r1_pattern}" -type f))
       local r2_files=($(find ${dir} -name "${r2_pattern}" -type f))
																        
       if [ ${#r1_files[@]} -eq 1 ] && [ ${#r2_files[@]} -eq 1 ]; then
    	   	R1_FASTQ="${r1_files[0]}"
        	R2_FASTQ="${r2_files[0]}"
       return 0
																						            fi	
   done

   return 1	
}
# Find R1 and R2 files
echo "Searching for FASTQ files for sample: ${SAMPLE_NAME} in directory: ${FASTQ_DIR}"
if find_fastq_files "${SAMPLE_NAME}" "${FASTQ_DIR}"; then
	echo "Found R1 file: ${R1_FASTQ}"
	echo "Found R2 file: ${R2_FASTQ}"
else
	echo "ERROR: Could not find matching R1 and R2 files for sample ${SAMPLE_NAME}"
        echo "Searched in directory: ${FASTQ_DIR}"
	echo "Expected patterns include:"
        echo "  ${SAMPLE_NAME}_R1*.fastq.gz and ${SAMPLE_NAME}_R2*.fastq.gz"
	echo "  ${SAMPLE_NAME}_1*.fastq.gz and ${SAMPLE_NAME}_2*.fastq.gz"
        echo "  ${SAMPLE_NAME}*_R1*.fastq.gz and ${SAMPLE_NAME}*_R2*.fastq.gz"
	echo "Available files in directory:"
	ls -la ${FASTQ_DIR}/*fastq.gz 2>/dev/null || echo "No .fastq.gz files found"
	exit 1
fi

echo "=== Starting ARIMA Hi-C processing with fastp ==="

# Step 1: Simple 5-base trimming setup
echo "Step 1: Setting up 5-base trimming..."
echo "Will trim first 5 bases from both R1 and R2 reads" 

# Step 2: Run fastp with 5-base trimming from 5' end
echo "Step 1: Trimming first 5 bases with fastp..." 
echo "Start fastp trimming: $(date)"

srun fastp \
       	--in1 ${R1_FASTQ} \
        --in2 ${R2_FASTQ} \
	--out1 ${OUTPUT_DIR}/${SAMPLE_NAME}_trimmed_R1.fastq.gz \
        --out2 ${OUTPUT_DIR}/${SAMPLE_NAME}_trimmed_R2.fastq.gz \
	--trim_front1 5 \
        --trim_front2 5 \
	--trim_tail1 0 \
	--trim_tail2 0 \
    	--cut_front \
    	--cut_tail \
    	--cut_window_size 4 \
    	--cut_mean_quality 20 \
    	--qualified_quality_phred 20 \
    	--unqualified_percent_limit 30 \
    	--length_required 20 \
    	--thread ${THREADS} \
    	--html ${OUTPUT_DIR}/${SAMPLE_NAME}_fastp_report.html \
    	--json ${OUTPUT_DIR}/${SAMPLE_NAME}_fastp_report.json \
    	--report_title "ARIMA Hi-C ${SAMPLE_NAME} - 5bp trimmed" \
    > ${OUTPUT_DIR}/fastp.log 2>&1 || error_exit "fastp trimming failed"

TRIMMED_R1="${OUTPUT_DIR}/${SAMPLE_NAME}_trimmed_R1.fastq.gz"
TRIMMED_R2="${OUTPUT_DIR}/${SAMPLE_NAME}_trimmed_R2.fastq.gz"

echo "fastp 5-base trimming completed: $(date)"
echo "Trimmed 5 bases from 5' end of both reads"

```


- juicer's arrowhead
- insulation plots over tad boundaries #todo 
- Straw - made by juicer ppl
	- dump contacts into counts format
	- most interactions within kb bin
	- input for multihicompare
- 5 KB is the limit
#todo:
- differential loop analysis
	- (multi hi chromosome)
	- pulling straw interactions for loops
	- low numbers
	- look at contacts
- hiccups calling long range loops

- aggregate pixels to call for differential loop analysis
- loops from multi hi compare w known


new tool:
- Mariner: https://academic.oup.com/bioinformatics/article/40/6/btae352/7685386
- https://github.com/EricSDavis/mariner
- loop might be shifted one bin
- #todo : set up mariner
- # main task




- sub compartment analsyis
- dchic
- strong a/b compartment
- less a compartment in mutant
- could try difhic
- juicer has own differential analysis - not replicate aware
	- good for big differences best not the best
  - #todo :
	  - p12 hi-c data
- #todo:
	- https://github.com/Taiji-pipeline/Taiji
		- transcription factors based on atac-seq and rna seq
- https://mail.google.com/mail/u/0/popout?ver=w4viue472nmr&search=inbox&th=%23thread-f:1845457786107077869&cvid=1
	- #todo ask questions about setup?
- TOBIAS:
	- merge into mutant file
	- normalizes somehow??
	- accessibility changes: more or likely diff bound transcription factors
	- works best when not much accessibility change
	- biased by the global accessibility 
	- if globally more is accessible
	- detect more transcription factors
	- more is dysregulated
	- rope in rna seq
	- want to back up tobias with rna seq -- so todo use taiji



# priorities
- #todo mariner
	- if not diff loops 
	- #todo if mariner contradicts hiccups, use cLoops
- #todo taiji
- see whats happening with r loading


- batch1/old batch
	- 600 mill reads/sample
	- loops 
	- do again for a billion reads
		- combined
- 25042 is first batch
- 8 hours default
- counts csd940
- 50GB ram
- max everything on juicer and nextflow
- job name and output name
- partition wont change
- hiccups need gpu
- with R:
- #todo: add logging to https://arc.net/l/quote/ilsnomvs
- #todo debug module loading

--

- https://github.com/aidenlab/JuicerTools
- juicer-pre: format pairs files from nextflow
- format files for juicer pre
- motif-finder
- set up new version of juicer at some point #todo 
	- diff hiccups
	- cpu
	- only detects loops within 3-8 mb
- #todo add readme files to ctea's directories
- maps is for hi-chip
	- singularity container is wrong 
	- more recent setup #todo low priority
- ks27ac and ubiquitinated histone at same site so cant assign reason for change
- hi-c cant tell you which it was, hi-chip will
	- hi-chip has big wig files with a lot of noise
	- methylation levels

- homer: tag directories
- input files: bam
	- no contact maps

- blacklisted files:
	- short read sequencing
	- 100-200 bps
	- blacklist repeated regions
	- cut out of analysis
	- long read sequencing 
		- with repeat sequencing
- blacklist.bed - 254 regions
- http://homer.ucsd.edu/homer/interactions2/index.html:
	- homer_scripts
- #todo: setup some kind of version control?
- superenhancers: modulate multiple promoters
	- maintain enhancer pool 
	- phase separation 
		- sub compartments that are accessible 
		- semi-accessible 
		- local interactions
	- multi-hi-compare:
		- enhancer enhancer contacts
	- cpg islands: repeats
		- methylation repression
		- promoter regions

- #todo : ctea/envs
- nextflow pipeline:
	- #todo: nfhic_pipeline.sb
- add logging to nextflow pipeline
- nfhic
- nf-hic
- trimmer: fastp

- bedpe files - contact files - 2 regions

- merge pairs - multiresolution from hiccups
	- 5, 10, 25 kb
	- long range loops only detectable at 25 kb
- converted hic files from cool are 1.6
- rsync all the trimmed
- start with juicer 2.0, then do 1.6
- strawr takes a lot of time
- goal for mariner:
	- counts list for each loop that is sensitive to shifts between bins
	- to use diffs like edgeR




- mariner - tonight and set that up
- saturday taiji pipeline - and email phd candidate about errors they were getting
- - p12 timepoint new seq data just came in -- next task