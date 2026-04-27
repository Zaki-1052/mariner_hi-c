
import subprocess
import os
from pathlib import Path
from deeptools_plotting import heatmap_plot


GENOME_SIZES = {'mm10': 2494787188, 'hg38': 2913022398}

def bam_coverage(bam_file,out_file,blacklisted_regions,genome='mm10'):
    # test if the data is PE
    samtools_func = 'samtools view -c -f 1 ' + bam_file
    out = subprocess.check_output(samtools_func,shell=True)
    out = int(out.decode('utf-8').split('\n')[0])
    if out == 0: extend_append = ' -e 100'
    else: extend_append = ' -e'
    genome_size = GENOME_SIZES[genome]
    func = 'bamCoverage --bam ' + bam_file + ' -o ' + out_file  + ' -bl ' + blacklisted_regions + ' -p 4 --normalizeUsing RPKM --effectiveGenomeSize ' + str(genome_size) + extend_append
    subprocess.run(func, shell=True)


# must provide either a bam_dict or bigWig_dict
    # bigWig_dict will be used directly, whereas files in a bam_dict will first be RPKM normalized using bamCoverage
# bed_dict is a dictionary containing paths to bed files of interest
def bed_pileup(bed_dict,out_dir,blacklisted_regions,bigWig_dict=None,bam_dict=None,out_name='bed_pileup_matrix',up_down=500,color_dict=None,vmax_groups=None,line_measure='mean',use_height=0.000075,pileup_type='referencePoint'):
    if (bigWig_dict is None) and (bam_dict is not None):
        
        bigWig_dict = {}
        for name,path in bam_dict.items():
            bigWig_out_name = os.path.dirname(path) + '/' + Path(path).stem + '_RPKM.bw'
            if not os.path.exists(os.path.dirname(path) + '/' + Path(path).stem + '_RPKM.bw'):
                print('performing RPKM normalization of bam files')
                bam_coverage(bam_file=path,out_file=bigWig_out_name,blacklisted_regions=blacklisted_regions)
            bigWig_dict[name] = bigWig_out_name

    bigWig_name_list = ''
    bigWig_file_list = ''
    for name,path in bigWig_dict.items():
        bigWig_name_list = bigWig_name_list + ' ' + name
        bigWig_file_list = bigWig_file_list + ' ' + path

    bed_name_list = ''
    bed_list = ''
    for bed_name,bed in bed_dict.items():
        bed_name_list = bed_name_list + ' ' + bed_name
        bed_list = bed_list + ' ' + bed


    # sort order is not preserved with different bigwigs i think - maybe i should fix that
    # maybe need to use --sortUsingSamples, but im not sure

    # not 100% sure on skipping zeros. by default they are included
    if pileup_type == 'referencePoint':
        out_path = out_dir + '/' + out_name + '_values'
        body_length = None
        func = 'computeMatrix reference-point --referencePoint center --missingDataAsZero -S ' + bigWig_file_list + ' -R' + bed_list + ' -o ' + out_dir + '/' + out_name + ' -bl ' + blacklisted_regions + ' -p 4 --outFileNameMatrix ' + out_dir + '/' + out_name + '_values -b ' + str(up_down) + ' -a ' + str(up_down) + ' --sortRegions keep --samplesLabel ' + bigWig_name_list
    elif pileup_type == 'scale-regions':
        out_path = out_dir + '/' + out_name + '_values_TSS-TES'
        body_length = 5000
        func = 'computeMatrix scale-regions --missingDataAsZero -S ' + bigWig_file_list + ' -R' + bed_list + ' -o ' + out_dir + '/' + out_name + ' -bl ' + blacklisted_regions + ' -p 4 --outFileNameMatrix ' + out_dir + '/' + out_name + '_values_TSS-TES -b ' + str(up_down) + ' -a ' + str(up_down) + ' --regionBodyLength ' + str(body_length) + ' --sortRegions keep --samplesLabel ' + bigWig_name_list
    subprocess.run(func, shell=True)

    heatmap_plot(values_path=out_path,
                 color_dict=color_dict,vmax_groups=vmax_groups,line_measure=line_measure,use_height=use_height,pileup_type=pileup_type,up_down=up_down,body_length=body_length)
