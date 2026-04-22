import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import bioframe
import cooler
import cooltools
import cooltools.lib.plotting
import json
from mpl_toolkits.axes_grid1 import Divider, Size

from plotting import factor_int,ax_selection,box
from statistics_functions import kruskal_wilcoxon


plt.style.use('seaborn-poster')

with open('/Users/tessapopay/example_data/custom_params.json') as f:
    custom_params = json.loads(f.read())

font_dict = {'big_title':10,'title':8,'labels':6,'axlab':7,'plt-text':5,'legend':6,'fonttype':'arial'}


# bed dataframes need to have the columns 'chrom','start','end' and if they have strand, it needs to be 'strand' not 'Strand'
# use smaller flank for bedpe than bed
# split_diagonal can only be used with a bed_dict
def mcool_pileup(mcool_dict,out_dir,out_name,flank=500000,bed_dict=None,bedpe_dict=None,over_expected=True,resolution=10000,split_diagonal=False,palette=None,v_range=None,genome='hg38'):
    arms = make_arms(genome=genome)

    if v_range is None: v_range = [-1,1]

    if bed_dict is not None:
        use_dict = bed_dict
        type = 'bed'
    elif bedpe_dict is not None:
        use_dict = bedpe_dict
        type = 'bedpe'
        split_diagonal = False

    if split_diagonal:
        for pair in range(0,len(mcool_dict),2):
            first = list(mcool_dict.keys())[pair]
            second = list(mcool_dict.keys())[pair + 1]

            bed_mtx_dict1,_ = mcool_get_matrix(mcool_file=mcool_dict[first],resolution=resolution,type='bed',flank=flank,bed_dict=use_dict,over_expected=over_expected,arms=arms)
            bed_mtx_dict2,_ = mcool_get_matrix(mcool_file=mcool_dict[second],resolution=resolution,type='bed',flank=flank,bed_dict=use_dict,over_expected=over_expected,arms=arms)
            
            new_mtx_dict = {}
            for bed_name in bed_mtx_dict1:
                matrix1 = bed_mtx_dict1[bed_name]
                matrix2 = bed_mtx_dict2[bed_name]
                
                bin_count = len(matrix1)
                out_matrix = matrix1.copy()
                for row in range(0,bin_count):
                    new_row = matrix2[row][0:row].tolist() + matrix1[row][row:bin_count].tolist()
                    out_matrix[row] = new_row
                new_mtx_dict[bed_name] = out_matrix

            pileup_plotting(bed_mtx_dict=new_mtx_dict,resolution=resolution,flank=flank,out_dir=out_dir,out_name=out_name + '_' + first + '_vs_' + second)
    else:
        for mcool_name,mcool_file in mcool_dict.items():
            if type == 'bed': bed_mtx_dict = mcool_get_matrix(mcool_file=mcool_file,resolution=resolution,type=type,flank=flank,bed_dict=use_dict,over_expected=over_expected,arms=arms)
            if type == 'bedpe': bed_mtx_dict = mcool_get_matrix(mcool_file=mcool_file,resolution=resolution,type=type,flank=flank,bed_dict=use_dict,over_expected=over_expected,arms=arms)
            pileup_plotting(bed_mtx_dict=bed_mtx_dict,flank=flank,resolution=resolution,out_dir=out_dir,out_name=out_name + '_' + mcool_name,over_expected=over_expected,v_range=v_range)


def mcool_get_matrix(mcool_file,resolution,flank,bed_dict,type,over_expected,arms,mean_array=True):
    clr = cooler.Cooler(mcool_file + '::/resolutions/' + str(resolution))
    if over_expected: expected = cooltools.expected_cis(clr, view_df=arms, nproc=2, chunksize=1_000_000)
    else: expected = None

    # if type == 'bed': use_cols = ['chrom','start','end']
    # if type == 'bedpe': use_cols = ['chrom','start','end']

    # worth only plotting the max intensity for a given genes
    # i just realized that the way I run this means it plots not the max intensity contact for a given gene across the multiple conditions, but whichever is the max intensity contact for the given condition
    # which is kinda useless. Basically I might be quantifying different contacts for the different conditions
    # will try mean instead i guess?
    # otherwise I will have to do it a different way
    bed_mtx_dict = {}
    for name,sites in bed_dict.items():
        stack = cooltools.pileup(clr, sites, view_df=arms, expected_df=expected, flank=flank)
        if 'strand' in sites.columns:
            mask = np.array(sites.strand == '-', dtype=bool)
            # the order of the stack being returned from pileup seems to be different than expected, so I have adjusted the code to account for this
            stack[mask, :, :] = stack[mask, ::-1, ::-1]

        if mean_array:
            mtx = np.nanmean(stack, axis=0)
            bed_mtx_dict[name] = mtx
        else: bed_mtx_dict[name] = stack.copy()

    if type == 'bedpe': return bed_mtx_dict
    if type == 'bed': return bed_mtx_dict

def pileup_plotting(bed_mtx_dict,flank,resolution,out_dir,out_name,v_range,over_expected=True,rescale_flank=None):
    # need to make sure cmap range is

    nrows,ncols,fig,axs,divider = pileup_plot_setup(subplot_count=len(bed_mtx_dict))
    
    
    i = 0
    order = list(bed_mtx_dict.keys())
    repeats = len(order)
    for row in range(nrows):
        for col in range(ncols):
            fig_ax = ax_selection(nrows=nrows,ncols=ncols + 1,axs=axs,row=row,col=col,divider=divider)
            if i + 1 > repeats:
                fig.delaxes(fig_ax)
                continue

            dataset = order[i]
            mtx = bed_mtx_dict[dataset].copy()

            i += 1
            if over_expected:
                im = fig_ax.imshow(
                    np.log2(mtx),
                    vmax = v_range[1],
                    vmin = v_range[0],
                    cmap='coolwarm',
                    interpolation='none')
            else:
                im = fig_ax.imshow(
                    np.log2(mtx),
                    cmap='coolwarm',
                    interpolation='none')

            fig_ax.set(title=dataset,xticks=[],yticks=[])

    fig_ax = ax_selection(nrows=nrows,ncols=ncols + 1,axs=axs,row=nrows-1,col=0,divider=divider)
    if flank is None and rescale_flank is not None:
        fig_ax.set(xticks=[mtx.shape[0]/(3/rescale_flank),mtx.shape[0]-mtx.shape[0]/(3/rescale_flank)],
            xticklabels=['start','end'],
            yticks=[],
            yticklabels=[],
            xlabel='relative position')

    if flank is not None:
        fig_ax.set(xticks=[0,flank/resolution,2 * flank/resolution],
            xticklabels=[str(int(-flank/1000)) + 'kb','',str(int(flank/1000)) + 'kb'],
            yticks=[0,flank/resolution,2 * flank/resolution],
            yticklabels=[],
            xlabel='relative position')

    for row in range(nrows):
        cb_ax = ax_selection(nrows=nrows,ncols=ncols + 1,axs=axs,row=row,col=ncols,divider=divider)
        if row == 0:
            cbar = fig.colorbar(im,orientation='vertical',cax=cb_ax)
            cbar.set_label('log2obs/exp',labelpad=-30)
            cbar.outline.set_color('black')
            cbar.outline.set_linewidth(0.25)
            cbar.ax.tick_params(width=0.25,color='black',labelcolor='black')
        else: fig.delaxes(cb_ax)

    plt.savefig(out_dir + '/' + out_name + '.png',dpi=300)
    plt.savefig(out_dir + '/' + out_name + '.pdf',dpi=300)
    plt.show()
    #plt.figure()

def pileup_plot_setup(subplot_count):
    sns.set_theme(font='arial',style='ticks',rc=custom_params)

    nrows,ncols = factor_int(subplot_count)

    ax_height,ax_width = 1,1

    fig_height = (ax_height + 0.4) * nrows
    fig_width = (ax_width + 0.25) * (ncols + 0.5)

    fig, axs = plt.subplots(nrows=nrows, ncols=ncols + 1,figsize=(fig_width,fig_height))

    sc = Size.Scaled(2)
    fh = Size.Fixed(ax_width)
    fv = Size.Fixed(ax_height)

    h = [sc, fh] * ncols + [Size.Scaled(1),Size.Fixed(0.1),sc]
    v = [sc, fv] * nrows + [sc]
    divider = Divider(fig, (0.0, 0.0, 1., 1.), h, v)

    return nrows,ncols,fig,axs,divider

def make_arms(genome='hg38'):
    # Use bioframe to fetch the genomic features from the UCSC.
    chromsizes = bioframe.fetch_chromsizes(genome)
    cens = bioframe.fetch_centromeres(genome)
    # # create a view with chromosome arms using chromosome sizes and definition of centromeres
    arms = bioframe.make_chromarms(chromsizes,  cens)
    return arms