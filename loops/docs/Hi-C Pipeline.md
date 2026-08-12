Hi-C Processing  
General linux commands:

- “cd” will change the directory to wherever you want to go  
- “ls” will list all the files/directories in the current directory that you’re in  
- “pwd” will show you the entire file path for the directory you’re currently in  
- “cntrl c” ends any command that is in process  
- “ls \-l” give u more info about a file and file size  
- “rm” deletes files  
- “cp” copy  
- “vim” to go into a script/or view a script  
- “cd ..” is to go backwards (out of 1 directory)  
- “df \-h” is to find the statistics of the server (useful when you think there is no space. Check ‘data’ category)  
- “clear” will clear the screen  
- “mkdir” \- make a new directory  
- cp \<from path\> \<to path\>  
- mv \<from path\> \<to path\>

Script Editor: (press “esc”) first

- “G” to go to bottom  
- “gg” to go to top  
- “ZZ” to exit vim and save and quit  
- Shift “I” how you change/edit the script  
- Esc-button – exits the edit mode  
- “:q\!” exit vim without saving

Screen:

- “screen \-S“ then a screen name to create a new screen  
- ctrl a \+ d to detach from screen  
- ctrl a \+ k → y to kill a screen  
- “screen \-ls” to list a screen  
- “screen \-r \[screen.name\]” to reattach the screen  
  - If multiple screens and want to reattach, then have to type in the full name

SLURM:

- “sacct \-b” for a brief listing of past jobs  
- “sacct \-l \-j \[\#job number\]” for information about your job (time, resources, etc.)  
- “Squeue \-u $USER” for status on your submitted jobs

Sequencing data is available on the FTP server (will be sent out by Cole):   
ftp://igm-storage2.ucsd.edu/20230929\_LH00320\_0051\_B22CYMFLT3/  
Username:  
ferguson  
Password:  
bjHl5OPB3C  
If you cannot access this, you may have to be on UCSD wi-fi or use the VPN.

**If more CPU/GPU hours are needed, ask Cole to exchange credits for more hours.**

#### Before Starting (only need to do the first time): 

1. Ask Cole to make an ACCESS account for you:  
   ![][image1]  
2. Log into the expanse portal at [portal.expanse.sdsc.edu](http://portal.expanse.sdsc.edu) (it may take a few days after Cole gives access to ACCESS to access the expanse portal). Use your ACCESS login through Globus OpenID. You may be prompted to set up Duo.  
3. Once you gain access, the portal should appear like this:  
   ![][image2]  
4. Click on “expanse Shell Access” to open up a shell.  
5. You will be assigned a login node either 01 or 02 – do not run any super computationally intensive processes on these login nodes. These are for organizing, file pushing, and basic processes.  
6. ‘cd /expanse/lustre/projects/csd940’ and check that there is a folder that corresponds to your account  
   ![][image3]

To enter SDSC from your local terminal

1. ssh [username@login.expanse.sdsc.edu](mailto:username@login.expanse.sdsc.edu)  
2. Then enter your password from access  
3. Then go to Duo and put in the 6-digit code from your cilogon

#### Retrieving/uploading files: 

For shallow sequencing files (\<4gb), you can connect to server and drag and drop fastq files to your local directory  
For deep sequencing files (\>20gb), you must offload the files via wget:

1. “wget ‘ftp://ferguson:bjHl5OPB3C@igm-storage.ucsd.edu/250714\_LH00444\_0375\_B22Y2JVLT3/file’ \-P /path/to/local/directory/”  
2. These are usually very big files, it is advisable that you offload to an HD. Cole will usually want this as backup

Importing/exporting from expanse:  
Note: Your home directory (cd \~/) on expanse is limited in space (100gb). I usually keep reference files, bed files, and backup scripts here. Fastq files and hic files should be kept in the lustre filesystem.  
There is also a scratch directory for each user in /expanse/lustre/scratch/$USER/

1. Three options:  
   1. scp \-r path/to/file [username@login.expanse.sdsc.edu](mailto:username@login.expanse.sdsc.edu):/expanse/lustre/projects/csd940/path/to/directory  
   2. rsync \-r path/to/file [username@login.expanse.sdsc.edu](mailto:username@login.expanse.sdsc.edu):/expanse/lustre/projects/csd940/path/to/directory  
      1. \--exclude '\*.tmp' to not transfer anything with .tmp  
      2. –include for vice versa  
      3. Ex: rsync \-av \--exclude='\*.bam' \--exclude='\*.pairs.gz' \--exclude='\*.pairs.gz.px2' \--exclude='\*balanced.txt' \--progress \--stats ctea@login.expanse.sdsc.edu:/expanse/lustre/projects/csd940/ctea/nftest/trimmed\_mut\_M1 /Volumes/Ferguson\_HD/Challana/HiC/  
   3. Pull the files from the FTP server in expanse.  
2. This will take a while for the fastq files – best to upload overnight.

#### Running juicer on shallow sequencing files: 

1. ‘cd /expanse/lustre/projects/csd940/ctea/HiC’  
2. Prepare a folder for your experiment (i.e. 250402\_Bap1)  
   1. For each sample, make a folder corresponding to the sample name (i.e. ctrl\_M1, ctrl\_M2, ctrl\_F3, etc.)  
   2. In each sample name folder, make a directory ‘fastq’ and place the corresponding R1 and R2 fastqs within each folder  
   3. Sample directory should look like 250402\_Bap \> ctrl\_M1 \> fastq \> fastq.R1, fastq.R2  
3. To run juicer, go to the HiC folder and run “sbatch juicer.sb \[sample\_name\] \[sample\_folder\]”  
4. This will submit a job. Check ‘squeue \-u $USER’ for the status.  
   1. Typically, you will either see the run status as “PD” (priority) or “R” (running). If the status is PD, then the run has not started yet and is on the queue. If there are a lot of jobs on queue, then it may be a little while before your job runs.   
      1. You can also click on “active jobs” from the expanse portal page   
   2. Monitor the run maybe for the first 10 minutes after you’ve submitted the job. If the run ends in \<10 minutes, chances are that the run has an error.  
   3. “vim juicer3.\[\#jobnumber\].out” to view the run. You may have to debug.  
   4. Shallow-sequenced samples should take \<6 hours to run through juicer.  
5. Once the run is complete for a sample, the sample name folder should contain 3 folders – fastq, aligned, and splits.  
6. Your summary statistics are in ‘aligned’ as ‘inter.txt’ or ‘inter\_30.txt’  
   1. Inter\_30.txt corresponds to aligned reads with a MAPQ score \> 30\.  
   2. Important things to look for, but consult this [guide](https://arimagenomics.com/wp-content/files/Bioinformatics-User-Guide-Arima-HiC-and-Arima-High-Coverage-HiC.pdf) for more details:  
      1. Aligned reads should be between 80-90%  
      2. PCR duplicates should not be \>20%  
      3. Below MAPQ threshold should ideally not be \>20%

![][image4]  
If QC metrics are all within expected range, then proceed with deep sequencing (600M-1B reads).  
If you ever want to view the contents of a singularity container, run:  
‘singularity exec $JUICER\_CONTAINER ls /’ where JUICER\_CONTAINER=/path/to/container

Some extra juicer notes 260113:

- What doesn’t work – Juicer eigenvector does not work at resolutions higher than 250000 BP, I’ve tried a few things and nothing’s worked to get higher. If you need this, use the cis.vecs.txt files from your nextflow output. Plot these alongside a K27ac bigwig to check that directionality is good.  
- What does work – diffhiccups, arrowhead, APA, compare, validate

#### Running nextflow-hic pipeline on deep sequencing files: 

Note: Keep in mind that we are running the nextflow pipeline instead of the juicer one. It is slightly faster and is based on Hic-pro’s pipeline. (aligns with bowtie2 and produces cool files, validpairs files, and pairix files that are handy for other analysis tools.

Trimming with fastp:

1. As done above, mkdir a folder in the ‘nf-hic’ folder for your fastq directory (i.e. 250402\_Bap1\_deepseq).  
2. Trim your files with fastp before running the pipeline – this will chew back some of the adapter sequence for better alignment.  
3. Run ‘sbatch fastp.sb \[sample\_name\] \[fastq\_dir\]’ (i.e. sbatch fastp.sb ctrl\_M1 250402\_Bap1\_deepseq)  
4. This will create a directory in your fastq directory called ‘trimmed’ containing your trimmed fastq files in your corresponding sample folder

Running nextflow:

1. After trimming, files should be structured as nf-hic/\[fastq\_dir\]/trimmed/\[sample\_name\]/trimmed\_{R1,R2}.fastq.gz  
2. Run ‘sbatch nfhic\_pipeline.sb \[sample\_name\] \[fastq\_dir\]  
   1. If you are running a second batch/lane (where you have multiple sets of fastq files corresponding to the same sample, you can concatenate those files with ‘cat B1fq.gz B2fq.gz \> combined.fq.gz’ or you can specify in your nextflow settings that you are providing multiple runs.  
   2. For the latter, you will need to trim the fastq files from both batches/runs, then specify the directory of the trimmed directory in the nfhic\_pipeline.sb script in ‘B2\_INPUT’. Make sure your B2 sample has the same sample name identifier.  
3. This will take \~25-30 hours to run for a standard 600M read library, maybe \~35-40 hours for a 1B read library.  
   1. Your output folder (Processed\_\[sample\_name\]) should be located within your experiment fastq directory and contain a few folders  
      1. Compartments \- A/B values \- I don’t typically use this  
      2. **Contact\_maps** \- This will contain pre-generated bins, cool files (both raw and cooler-balanced) and txt files for each resolution; mcool files are raw  
      3. Distance\_decay \- distance decay plots; honestly this isn’t working, we can call these on another tool if we want to  
      4. Fastqc \- individual fastqc  
      5. Multiqc \- contains multiqc plots and your report – it is important to look at your **multiqc report** to look at your hi-c metrics or whether you need to tweak your pipeline conditions  
      6. Hicpro \- contains 3 folders  
         1. Mapping \- results of **bowtie2 alignment** – this will contain the bam files needed for Homer tag directories  
         2. Stats \- distribution of reads (these stats are also reported in the multiqc report file)  
         3. **Valid\_pairs** \- this contains your .allValidPairs file that can be used to call contact maps as well as pairix files that can be formatted for juicerpre  
      7. Pipeline\_info \- this is going to have the same info as would be reported in your slurm outs file  
      8. Tads \- hicExplorer folder will contain results from hicExplorer’s [hicFindTADS](https://hicexplorer.readthedocs.io/en/latest/content/tools/hicFindTADs.html) and your insulation folder will contain insulation values for each bin size  
         1. The bin sizes for this can be pre-specified in the nfhic\_pipeline.sb script – currently, it is set to 10kb, 20kb, and 40kb

#### Juicerpre & Hiccups: 

If you want to run hiccups (loops) or arrowhead (TADs), your hic files need to be in juicer v2.0 format. In order to generate this, you need to run juicerpre on your pairs output from the nextflow pipeline.

To do this, you will need to pre-format your pairs files from hicpro/valid\_pairs/pairix 

1. ‘cd /expanse/lustre/projects/csd940/ctea/HiC/juicer\_scripts’  
2. Vim ‘format\_4DN.sb’ and ‘juicer\_pre.sb’ and change ‘Sample\_folder’ to the directory that contains your processed nextflow folders  
3. Run “sbatch format\_4DN.sb \[sample\_name\]”  
   1. This will produce a juicer folder within the processed\_\[sample\_name\] directory with a \[sample\_name\]\_short.pairs to be used for juicerpre  
4. Now run “sbatch juicer\_pre.sb \[sample\_name\]”  
   1. This will produce a juicer 2.0 .hic file in the juicer folder within the processed\_\[sample\_name\] directory  
5. Note: if you want these files to be KR normalized and only KR-normalized (needed for some tools), then add ‘-k KR’ as a flag. Otherwise, it will run the standard VC, VC\_SQRT, GW\_SCALE, SCALE, INTER\_SCALE

Note: if you would like to generate a merged file combining counts from all the replicates together, the best way you can do this is by merging the short.pairs file that is output from format\_4DN.sb

1. ‘cd /expanse/lustre/projects/csd940/ctea/nf-hic’  
2. Vim into the ‘sortmerge\_contacts.sh’ script to change your base directory and output directory as needed  
3. ‘sbatch sortmerge\_contacts.sh \[ctrl\_rep1\_short.pairs\] \[ctrl\_rep2\_short.pairs\] \[ctrl\_rep3\_short.pairs\] \[merged\_ctrl\]’ – change the names accordingly, the output file to be used for juicerpre will then be sorted\_merged\_ctrl\_short.pairs

Now to run hiccups on these files:

1. ‘cd /expanse/lustre/projects/csd940/ctea/hiccups/’  
2. Prepare a exp folder with your hi-c files generated from juicerpre  
   1. Ideally, you should symlink the .hic files with “ln \-s /path/to/hic .”  
3. Run hiccups with “sbatch hiccups.sb path/to/hicfiles”  
4. This will produce a folder inside ‘hiccups/hiccups\_results’ corresponding to your exp folder, with folders corresponding to each sample within the exp folder

#### Cooler functions:

A few functions for manipulating cool files are located within nf-hic/cooler\_scripts  
Challana’s note \- I do not currently have the time to make these batch scripts the most flexible. Vim each script and change the directories as needed.

* [cbalance.sb](https://cooler.readthedocs.io/en/latest/cli.html#cooler-balance) \- If you want to re-balance raw cool files via cooler’s version of ICE  
* [cdump\_only.sb](https://cooler.readthedocs.io/en/latest/cli.html#cooler-dump) \- cooler’s version of juicer straw, never really ended up using this – just use juicer straw on juicerpre files  
* [cload\_validpairs.sb](https://cooler.readthedocs.io/en/latest/cli.html#cooler-cload) \- loads contact maps from validpairs file  
  * If there was a resolution you didn’t specify in your nextflow-hic configuration, then you can load the raw contact map from the validpairs file here  
* [czoomify.sb](https://cooler.readthedocs.io/en/latest/cli.html#cooler-zoomify) \- generates a multi-resolution cool file (mcool) from the resolution of your specified cool file  
  * Specify resolutions with \-r flag  
  * Highest resolution will be from the cool file that you give (ideally your 1kb.cool file)  
* [cinfo.sb](https://cooler.readthedocs.io/en/latest/cli.html#cooler-info) \- provides the info and metadata of the provided cool file

#### Homer \- Differential TAD & PC Analysis:

For the most part, I am following the workflow described [here](http://homer.ucsd.edu/homer/interactions2/index.html)  
Before starting, you will need to follow these [instructions](http://homer.ucsd.edu/homer/introduction/install.html) for installing homer and adding the homer/bin directory to your exec path in your \~/.bash\_profile file

Generating tag directories

1. If two different batches are processed through nextflow, you will need to merge the \_0\_bwt2pairs.bam and \_1\_bwt2pairs.bam files together with “merge\_and\_splitbam.sb”. Change directories as needed. This script could be adjusted to also run homer’s maketagdir script and save a step.  
2. In order to process files through homer, you will need to generate tag directories using the bam/sam files given through juicer/nextflow. These bam files will need to be split into their respective R1 and R2 alignment files.  
3. ‘cd /expanse/lustre/projects/csd940/ctea/homer’  
4. Prepare a folder with your all your bam files for each sample (these can be symlinked)  
5. Run ‘sbatch unmergebam.sb \[folder\_name\]’ → this will create a “split\_bams” directory within this folder with R1 and R2 bam files for each sample.  
6. To make a tag directory for a sample, go to homer\_scripts and edit ‘ALIGNMENT\_DIR’ in maketagdir.sb to match the location of your split bams file  
7. Run ‘sbatch maketagdir.sb \[sample\_name\]’ → this will create a folder ‘tags’ in your ALIGNMENT\_DIR with the sample name as a folder inside of it. Inside your sample folder is the tag directory in a ‘tagdirs’ folder.  
   1. Be sure sample\_name is an indicator in the name of your bam file

TAD analysis – not the most user-friendly atm, you will have to go into each script and change the directories as needed

1. Run ‘sbatch [findtad.sb](http://findtad.sb) \[sample\_name\]’ to find the TADs for an individual sample.  
2. Run ‘sbatch [mergetad2.sb](http://mergetad2.sb)’ to merge the TADs from each tag directory in your ALIGNMENT\_DIR folder from above. This is needed to create a consensus TAD set which you can name in the script.  
   1. Specify the folder with your sample tag directories in TAD\_DIR  
3. Run ‘sbatch [scoretad.sb](http://scoretad.sb)’ to calculate an inclusion ratio (i.e. interaction density) for a TAD from each sample.  
   1. Specify the folder with your sample tag directories in BASE\_DIR  
   2. Specify your consensus TAD set from step 2 as MERGED\_TAD  
   3. This will create a file of your chosen name which is what will be used as the counts table for differential analysis  
4. Run ‘sbatch [getDiffExpression.sb](http://getDiffExpression.sb)’ to perform your TAD differential analysis  
   1. Specify your scores table created step 4 as TAD\_SCORES  
   2. This will produce a txt file with significance values that you can plug into a volcano plot in R 🙂

PC analysis – this will essentially be the same as your TAD analysis but with PC1 values

1. sbatch [PC.sb](http://PC.sb) \[sample\_name\]  
2. sbatch [diffcompartments.sb](http://diffcompartments.sb)  
3. To compare within regions, run ‘sbatch bedgraph\_analysis.sb \[ctrl\_samplename\] \[mut\_samplename\] \[peak\_region\]

Analysis for compaction is in ‘compaction\_scripts’ but I don’t anticipate us using them that much

#### GENOVA:

GENOVA works from cool files and juicer 1.6 hic files. Keep in mind that it is not replicate-aware. It only performs single ctrl-mut comparisons. You can run analyses on separate replicates to more than one instance of a phenotype, but for your final analysis for publication/presentation, I would recommend combining the cool or hic files you plan on using into one replicate.

- The best way to do this merge will be through hictk merge. There is a script for this in ctea/HiC/hictk/[merge.sb](http://merge.sb), but you can make one easily and input your cool/hic files

GENOVA is helpful for:

- Loci visualizations \- insulation, side-by-side comparisons, difference plots  
- Aggregate TAD analysis  
- Aggregate region analysis \- helpful because you can feed in any particular peak file (i.e. dysregulated peaks identified from C\&R)  
- Aggregate loop analysis  
- Compartment (CSS) analysis \- not as set up here though

Currently (10/21/25), GENOVA is set up for the aws instance. If someone would like to set it up on expanse, be my guest. I have also set up a developmental version of GENOVA since I was running into errors with the norm function in GENOVA. The norm function is still kind of funky, so right now, I just feed cooler-balanced cool files and call it a day.

Note to self 260120:

- Currently for some reason there’s an error in load\_contacts using the GENOVA\_hic script, but it loads if I open up an R shell and just run the same commands  
- Okay, I’ve now made a new conda env GENOVA\_2.0 which runs the most updated GENOVA and now accepts balancing=TRUE to use the KR weight column in juicer KR-normalized hic files  
  - Good grief this was such a crazy run-around  
- Rscript load\_contacts\_2.R hic/resorted\_mut.hic \--resolution=10000 \--sample\_name=Mut \--colour=red \--scale\_cis \--balancing

If you are running very deep sequencing files at higher bins (900M-1B reads, 5-10kb res), then the current RAM on the instance is not high enough to load all the files in. You will need to stop the instance, go to the EC2 dashboard, and change the instance type to increase the RAM:   
Actions \> Instance settings \> Change instance type – then in “New instance type”, change from m5.4x large to m5.8x large → then hit “Change instance type” at the bottom

This **doubles the price of instance hours** so be sure to change this instance type back down to m5.4x large when you’re done and also note to the slack aws-instance channel that you are changing the instance type so people are a bit more judicious about leaving the instance on

**Loading contacts**

1. Export your balanced cool files from expanse to your local and then import into the AWS instance in ‘/data2/rs\_256/hic/cool’  
2. ‘cd /data2/rs\_256/hic/GENOVA’  
3. Run ‘conda activate R\_env’  
4. Run ‘Rscript GENOVA.R’  
   1. This will prompt you to enter the file path to your wild type and mutant files as well as your resolution. Currently for flexibility, I have this set up so that you have to enter the full filepath, so ‘cool/filename.cool’, but this could be changed if it gets too annoying  
   2. Enter your resolution as a full number (i.e. 5000, 10000, etc., not 5kb, 10kb)  
5. This will generate a [chromosome matrix](https://rdrr.io/github/robinweide/GENOVA/man/chromosome_matrix.html), [Tad+n jitter plot](https://rdrr.io/github/robinweide/GENOVA/man/intra_inter_TAD.html), and [PE-SCan](https://rdrr.io/github/robinweide/GENOVA/man/PESCAn.html) plot in the ‘initial\_outs’ directory.  
6. This will also generate .RData objects for both your wild type and mutant cool files in your cool directory (i.e. ctrl.cool.RData)

**Loci Visualizations** \- There are two different options here: Loci\_single.R and Loci\_batch.R. I’m almost always just running Loci\_batch.R because why not just test out a few loci?

1. Activate the R\_env conda environment if not already activated.  
2. Run ‘Rscript Loci\_batch.R’ or Loci\_single.R from the GENOVA directory.  
3. List the file paths when prompted for your ctrl and mut **RData** files.  
4. Enter loci as prompted. Be sure you are entering chromosomes in as “chr\#” and not “\#”  
   1. Loci\_single.R will ask for the chromosome, left, and right bounds separately  
   2. Loci\_batch.R will ask for you to input each loci as “chr\# \[left\] \[right\]  
5. Run and this will generate a side-by-side pyramid plot, a differential pyramid plot, and an insulation plot for each loci.  
   1. Pyramid plots will go in the ‘pyramids’ directory. Insulation plots will go in the ‘Insulation’ directory

![][image5]

**Aggregate Peak/Loop Analysis**  
Your peaks/loops should be a bedpe file with two genome coordinates. Most likely you will be getting this from juicer hiccups.

1. Import your .bedpe file into /data2/rs\_256/hic/GENOVA/regions/loops  
2. Only the files of interest for analysis should be in this folder. Otherwise move them from ‘loops’ to ‘regions’.  
3. Activate the R\_env conda environment if not already activated.  
4. Run ‘Rscript APA.R’ from the GENOVA directory.  
5. List the full file paths when prompted for your ctrl and mut **RData** files as well as the corresponding resolution.  
6. In the GENOVA/APA directory, this will output png files showing the pixel enrichment around the loops for each bedpe file in your loops directory.

![][image6]  
Note: If you get a “Error in value\[\[1L\]\] : subscript out of bounds” error, this likely means that the chr column of your cool file contacts are \#, not ‘chr\#’. In this case, you will need to convert the chr column of your bed/bedpe files to \#. There is a line in both this script and the ATA script that you can un\# to do this.

**Aggregate TAD Analysis**  
Your TAD file will be a bed file. This will likely be generated from homer or juicer arrowhead. Usually this is a nice way to validate the results of homer \- are differential TADs identified in homer indeed higher/lower in their contacts

1. Import your bed file into /data2/rs\_256/hic/GENOVA/regions/TAD\_analyze.  
2. Only the files of interest for analysis should be in this folder. Otherwise move them from TAD\_analyze to ‘regions’.  
3. Activate the R\_env conda environment if not already activated.  
4. Run ‘Rscript ATA.R’ from the GENOVA directory.  
5. List the full file paths when prompted for your ctrl and mut **RData** files as well as the corresponding resolution.  
6. In the GENOVA/ATA directory, this will output png files showing the enrichment within and around the TADs for each bed file in your TAD\_analyze directory.

![][image7]

**Aggregate Region Analysis**  
Your region file can be any particular bed file – usually I run this analysis on differential peak regions identified from CUT\&RUN.

1. Import your bed file into ‘/data2/rs\_256/hic/GENOVA/regions/analyze’.  
2. Only the files of interest for analysis should be in this folder. Otherwise move them from ‘analyze’ to ‘regions’.  
3. Run ‘Rscript Region.R’ from the GENOVA directory.  
4. List the full file paths when prompted for your ctrl and mut **RData** files as well as the corresponding resolution.  
5. In the GENOVA/Region\_Analysis directory, this will output png files showing the enrichment and insulation within and around the regions for each bed file in your ‘analyze’ directory.

![][image8]![][image9]

#### MultiHiCompare:

Detailed documentation for multiHiCompare can be found [here](https://www.bioconductor.org/packages/release/bioc/vignettes/multiHiCcompare/inst/doc/multiHiCcompare.html).

Currently, we have our pipeline set up to run Fastlo normalization as it is quite a bit faster and less memory intensive, but we could try out cyclic loess maybe if we move some scripts over to the expanse lustre server.

I also currently have straw set up in the dchic conda environment.

**Straw**  
You have a few different options for your input for straw \- I typically use the .hic files. You can either convert your cool files to hic files via hicexplorer or you can run juicerpre on your pairs files. I would recommend the latter as you will already be running juicerpre for Hiccups.  
Alternatively, you can format your cool files as done in section 2.2.2 from the bioconductor vignette.

1. cd /expanse/lustre/projects/csd940/ctea/HiC/straw\_scripts  
2. Vim ‘straw.sb’ and change BASE\_DIR to the directory of your .hic files and RES to the resolution of interest (usually 5, 10, 25kb)  
3. Now run ‘straw.sb \[sample\_name\]’ where your hic file of interest is \[sample\_name\].hic  
4. This will create a ‘straw’ directory in your base directory and a directory for your resolution within your straw directory. Your upper sparse triangle matrix will be in your resolution directory as ‘contacts\_\[sample\_name\].\[res\].txt’  
5. Export this to your local and import to the AWS instance in ‘/data2/rs\_256/hic/straw’

**Running multiHiCompare on EC2:**  
The instance type we typically have set is m5.4x large. You will need to up the instance type to m5.8x large as noted above for GENOVA in order to have enough RAM to load contacts.   
(11/03/25 CT: We should probably set this up in the SDSC server, but I originally had this set up on AWS as it was easier to finneagle around and troubleshoot. Now that it’s working though, we can consider migrating these steps over to Expanse.)

1. cd /data2/rs\_256/hic/straw  
2. These contact files are very, very large and multiHiCompare calls its significant contacts only within chromosomes. Thus, the higher the resolution, the more you will want to split up your contact files into chunks of chromosomes.  
   1. Make sure only the contact files for your samples of interest are in this folder. If someone else’s contacts are present, move them into a separate folder or delete them  
   2. For 10kb straw matrices, run ./split\_contacts.sh → this will split contacts by 5 chromosomes  
   3. For 5kb straw matrics, run ./split\_contacts\_5kb.sh → this will split contacts into 3 chromosomes each  
   4. For 25kb straw matrices, run you can edit either of the following above scripts to split into 10 chromosomes, or you can just split contacts by 5 chromosomes as in 10kb to make your life easier  
   5. This will produce folders for each chromosome subset in the straw directory.  
3. Conda activate R\_env  
4. In ‘data2/rs\_256/hic’, run ‘Rscript multiHiCompare.R \[one folder of subdivided contacts (i.e. chr1\_3)’. Repeat for every chromosome subset.  
   1. Currently the script is set up to operate with 3 controls and 3 mutants, but this can be made more flexible in the future  
   2. Once the run is complete, this will create an output folder named after your chromosome subset with few visualizations (I usually only look at the manhattan plot) and a ‘multihicompare\_results.txt’ file which can be used for downstream analysis of significant contacts.

#### dchic:

For running on the Expanse server:

1. cd /expanse/lustre/projects/csd940/zalibhai/dchic/DCHIC  
2. Make 3 folders: 5kb, 10kb, 25kb, mkdir straw in each; symlink your contact txt maps to “.”  
   1. \`ln \-s /path/to/contacts\_\[sample\_name\].\[res\].txt /expanse/lustre/projects/csd940/zalibhai/dchic/DCHIC/{res}kb/\`  
3. Run split contacts like before, then: mv chr\* straw/  
4. There are three existing multihiCompare.R scripts in each res folder, their wd is set to /{res}kb  
   1. Make sure each working directory is vimmed to /expanse/lustre/projects/csd940/zalibhai/dchic/DCHIC/5kb etc  
   2. This can be edited to take an argument but I think it’s simpler to separate them  
5. Then just run ./[submit.sh](http://submit.sh) {resolution folder}  
   1. Eg: ./[submit.sh](http://submit.sh) 5kb  
6. It will submit \`run\_multihic.sb\` for you for each chrom in that resolution  
   1. This can be edited to run all three resolutions at once but it might overwhelm slurm

To document:

- Hictk  
- hicexplorer

**Hi-chIP (In progress, not fully set up)**  
![][image10]

singularity exec \--bind /scratch,/expanse /cm/shared/apps/containers/singularity/arimamaps/arima\_maps\_hichip\_v0.1.sif /MAPS/bin/Arima-MAPS\_v2.0.sh \-C 0 \-m HichIP\_referencepeak/P51\_K119ub\_ctrl\_intersect.bed \-I 250402\_HiChIP/$1 \-O $SLURM\_SUBMIT\_DIR/$1/output \-o mm10 \-b references/Mus\_musculus\_assembly10.fa  
sta \-t 8 \-f 1

singularity exec \--bind /scratch,/expanse /cm/shared/apps/containers/singularity/arimamaps/arima\_maps\_hichip\_v0.1.sif /MAPS/bin/Arima-MAPS\_v2.0.sh \-C 0 \-I /M  
APS/test\_data/fastq/Arima-MAPS-test \-O $SLURM\_SUBMIT\_DIR/output \-m /MAPS/test\_data/ENCFF247YHM.UW.bed \-o hg19 \-b /MAPS/test\_data/hg19/hg19.fa \-t 8 \-f 1

![][image11]











