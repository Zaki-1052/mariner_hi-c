
import os
import pandas as pd
import seaborn as sns
import json
import numpy as np
from plotting import ax_selection,randomize_hex

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'custom_params.json')) as f:
    custom_params = json.loads(f.read())

sns.set_theme(font='arial',style='ticks',rc=custom_params)

# only accepts 'mean' for line measure
def heatmap_plot(values_path,use_height,pileup_type,color_dict=None,vmax_groups=None,line_measure='mean',bed_color_dict=None,up_down=None,body_length=None,xticklabels=None):
    df = pd.read_csv(values_path,sep='\t',header=None,nrows=1)
    datasets = [i.replace('#','') for i in df.loc[0,:].tolist()]
    dataset_names = [i.split(':')[0] for i in datasets]
    df = pd.read_csv(values_path,sep='\t',header=None,nrows=1,skiprows=1)
    down,up,bin_size = int(df.loc[0,0].replace('#downstream:','')),int(df.loc[0,1].replace('upstream:','')),int(df.loc[0,3].replace('bin size:',''))
    df = pd.read_csv(values_path,sep='\t',header=None,nrows=1,skiprows=2)
    data_cols = [i for i in df.loc[0].tolist() if i not in datasets]  
    bam_order = [i for i in df.loc[0].unique() if i not in datasets]  
    df = pd.read_csv(values_path,sep='\t',header=None,skiprows=3)

    n_bins = df.shape[1]/len(list(set(data_cols)))

    df['mean'] = df.astype(float).mean(axis=1)
    used_rows = 0
    for dataset in datasets:
        set_name,nrows = dataset.split(':')
        df.loc[used_rows:used_rows + int(nrows) - 1,'group'] = set_name
        used_rows = used_rows + int(nrows)

    df.sort_values('mean',ascending=False,ignore_index=True,inplace=True)
    df.drop(columns='mean',inplace=True)

    df.loc[-1] = data_cols + [np.nan]
    df.index = df.index + 1  # shifting index
    df.sort_index(inplace=True)


    all_dict = {}
    vmax_dict = {}
    line_dict = {}
    line_max_dict = {}
    for bam in bam_order:
        vmax_dict[bam] = np.percentile(df.loc[1:,df.loc[0] == bam].to_numpy(),98.0,axis=(0,1))
        all_dict[bam] = {}
        line_dict[bam] = {}
        used_rows = 1
        for dataset in dataset_names:
            all_dict[bam][dataset] = df.loc[df['group'] == dataset,df.loc[0] == bam]
            if line_measure == 'mean': line_dict[bam][dataset] = all_dict[bam][dataset].mean().tolist()
        line_dict[bam] = pd.DataFrame(line_dict[bam])
        value_vars = list(line_dict[bam].columns)
        line_dict[bam]['index'] = range(int(-n_bins/2 * bin_size + bin_size/2),int(n_bins/2 * bin_size + bin_size/2),bin_size)
        line_dict[bam] = line_dict[bam].melt(id_vars='index',value_vars=value_vars,var_name='bed',value_name='values')
        line_max_dict[bam] = line_dict[bam]['values'].max() * 1.05

    if vmax_groups is not None:
        for group in vmax_groups:
            group_max = []
            group_max_line = []
            for bam,vmax in vmax_dict.items():
                if bam in group: group_max.append(vmax)
            for bam,vmax in line_max_dict.items():
                if bam in group: group_max_line.append(vmax)
            group_max = max(group_max)
            group_max_line = max(group_max_line)
            for bam in group:
                vmax_dict[bam] = group_max
                line_max_dict[bam] = group_max_line

    fig,axs,divider = plot_setup(dataset_nrows=[int(i.split(':')[1]) for i in datasets],bam_ncols=len(list(set(data_cols))),use_height=use_height)

    if color_dict is None:
        color_dict = {}
        for bam in bam_order: color_dict[bam] = 'Blues'
    
    if bed_color_dict is None:
        bed_color_dict = {}
        for dataset in dataset_names:
            color_pick = randomize_hex()
            while color_pick in bed_color_dict.values(): color_pick = randomize_hex()
            
            bed_color_dict[dataset] = color_pick


    col_rep = 0
    for bam in bam_order:
        #if col_rep == 2: break
        if (col_rep != len(list(set(data_cols))) - 1) or (len(datasets) == 1): legend = False
        else: legend = True
        ax = ax_selection(nrows=len(datasets) + 2,ncols=len(list(set(data_cols))),axs=axs,row=0,col=col_rep,divider=divider)
        line = sns.lineplot(line_dict[bam],x='index',y='values',ax=ax,hue='bed',legend=legend,palette=bed_color_dict)
        line.set_xlim(-n_bins/2 * bin_size,n_bins/2 * bin_size)
        if pileup_type == 'referencePoint':
            line.set_xticks([-n_bins/2 * bin_size,0,n_bins/2 * bin_size])
            if xticklabels is not None:
                line.set_xticklabels(xticklabels)
            else:
                line.set_xticklabels([str((-n_bins/2 * bin_size)/1000) + 'kb',0,'+' + str((n_bins/2 * bin_size)/1000) + 'kb'])
        elif pileup_type == 'scale-regions':
            print(line.get_xticks())
            set_val = -n_bins*(up_down*2 + body_length)
            line.set_xticks([-up_down,body_length/2])
            line.set_xticklabels(['TSS','TES'])

        line.set_title(bam)
        line.set_ylim([0,line_max_dict[bam]])

        if (col_rep == len(list(set(data_cols))) - 1) and (len(datasets) != 1):
            ax.legend(title=None,loc='center left', bbox_to_anchor=(1, 0.5),labelspacing=0.2,fontsize=5)

        if col_rep != 0:
            ax.set_ylabel('')
            ax.set_xticklabels([])
            line.set_xlabel('')

        row_rep = 1
        for dataset in datasets:
            use_set = dataset.split(':')[0]
            use_df = all_dict[bam][use_set].copy().astype(float).reset_index(drop=True)
            #use_df.columns = range(int(-n_bins/2 * bin_size + bin_size/2),int(n_bins/2 * bin_size + bin_size/2),bin_size)

            cmap = sns.color_palette(color_dict[bam],as_cmap=True)

            ax = ax_selection(nrows=len(datasets) + 2,ncols=len(list(set(data_cols))),axs=axs,row=row_rep,col=col_rep,divider=divider)

            if use_df.shape[0] < 1000: interpolation = 'nearest'
            else: interpolation = 'bilinear'

            heat = ax.imshow(use_df, interpolation=interpolation, cmap=cmap,vmin=0,vmax=vmax_dict[bam],aspect='auto')
            ax.tick_params('y',which='both',left=False,labelleft=False)
            ax.tick_params('x',which='both',bottom=False,labelbottom=False)
            if col_rep == 0: ax.set_ylabel(dataset.split('.bed')[0],rotation=90,ha='center',va='center_baseline')

            row_rep += 1

        cbar_ax = ax_selection(nrows=len(datasets) + 2,ncols=len(list(set(data_cols))),axs=axs,row=len(datasets) + 1,col=col_rep,divider=divider)
        cbar = fig.colorbar(heat,orientation='horizontal',cax=cbar_ax)
        cbar.outline.set_color('black')
        cbar.outline.set_linewidth(0.25)
        cbar.ax.tick_params(width=0.25,color='black',labelcolor='black')

        col_rep += 1

    fig.savefig(values_path + '.png',dpi=300)
    fig.savefig(values_path + '.pdf',dpi=300)

def plot_setup(dataset_nrows,bam_ncols,use_height):
    sns.set_theme(font='arial',style='ticks',rc=custom_params)

    from mpl_toolkits.axes_grid1.axes_size import Fixed, Scaled
    from mpl_toolkits.axes_grid1 import Divider, Size
    import matplotlib.pyplot as plt

    ax_height_list = [use_height * i for i in dataset_nrows]
    ax_width = 0.5

    fig_height = sum(ax_height_list) + 2
    fig_width = (ax_width + 0.75) * bam_ncols

    fig, axs = plt.subplots(nrows=2 + len(dataset_nrows),ncols=bam_ncols,figsize=(fig_width,fig_height))

    sc = Size.Scaled(1)
    fh = Size.Fixed(ax_width)

    h = [Size.Fixed(0.5)] + [Size.Fixed(ax_width),Size.Fixed(0.2)] * (bam_ncols) + [Size.Fixed(1)]
    v = [Size.Fixed(0.2),Size.Fixed(0.05)]
    for ax_height in reversed(ax_height_list): v = v + [Size.Fixed(0.05),Size.Fixed(ax_height)]
    v = v + [Size.Fixed(0.15),Size.Fixed(ax_width),Size.Fixed(0.1)]

    divider = Divider(fig, (0.0, 0.0, 1., 1.), h, v)

    return fig,axs,divider