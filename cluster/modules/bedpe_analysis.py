#!/usr/bin/env python

import pandas as pd
import subprocess
import copy
import os

import plotting
from cluster_tools import comparison_type,sort_by_strings
from statistics_functions import kruskal_wilcoxon


def loop_size(out_dir,bedpe_df=None,bedpe_dict=None,title=None,logY=False):

    bedpe_df_cp = pd.DataFrame()
    data_names,data_list = [],[]
    for sample in bedpe_dict:
        bedpe = bedpe_dict[sample].copy()
        bedpe[['x1','x2','y1','y2']] =  bedpe[['x1','x2','y1','y2']].astype(int)
        bedpe['loop size (kb)'] = (bedpe[['x1','y1']].max(axis=1) - bedpe[['x1','y1']].min(axis=1))/1000
        bedpe['sample'] = sample
        bedpe_df_cp = pd.concat([bedpe_df_cp,bedpe])
        data_names.append(sample)
        data_list.append(bedpe['loop size (kb)'].tolist())
    stats_df = kruskal_wilcoxon(data_names=data_names,data_list=data_list)
    stats_df.to_csv(out_dir + '/loop_size.stats.txt',sep='\t',header=True,index=False)

    comparison = comparison_type(data_cols=list(bedpe_dict.keys()))[0]
    if comparison == 'multiple': hue_col,palette = None,None
    if comparison == 'pairwise':
        bedpe_df_cp['treatment'] = bedpe_df_cp['sample'].str.split('_',expand=True)[0]
        bedpe_df_cp['sample'] = bedpe_df_cp['sample'].str.split('_',expand=True)[1]
        hue_col = 'treatment'
        palette = ['darkgrey', 'forestgreen']

    plotting.box(melted_df=bedpe_df_cp,
          xcol='sample',
          ycol='loop size (kb)',
          out_dir=out_dir,
          measure='loop size (kb)',
          title=title,
          hue_col=hue_col,
          palette=palette,
          out_name='loop_size',
          logY=logY
          )

    plotting.strip(melted_df=bedpe_df_cp,
          xcol='sample',
          ycol='loop size (kb)',
          out_dir=out_dir,
          measure='loop size (kb)',
          title=title,
          hue_col=hue_col,
          palette=palette,
          out_name='loop_size_strip',
          logY=logY
          )

# type can be 'anchor_intersect','within_intersect', or 'anchor_closest'
# if wanted, provide a FPKM_df containing columns at leat of 'gene_id','gene_name', and 'FPKM', and this will return
    # the FPKM of the intersecting genes
def bedtools_annotation(out_dir=None,type='anchor_intersect',bedpe_dict=None,FPKM_df=None,temp_dir='~/temp'):
    os.makedirs(temp_dir,exist_ok=True)

    default_bed_name = 'gene promoters'
    default_gene_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'mm10_knownGene_pp.bed')

    gene_file_path = default_gene_path
    B_bed_cols = ['chr','start','end','gene_id','idk','strand','gene_name']

    bedpe_df_cp = pd.DataFrame()
    analysis_type = 'dict'
    annotated_bedpe_dict = copy.deepcopy(bedpe_dict)
    for sample in annotated_bedpe_dict:
        annotated_bedpe_dict[sample]['sample'] = sample
        bedpe_df_cp = pd.concat([bedpe_df_cp,annotated_bedpe_dict[sample]])

    A_bedpe_cols = bedpe_df_cp.columns
    A_bedpe_cols = [col for col in A_bedpe_cols if col not in ['chr1','x1','x2','chr2','y1','y2','sample']]

    bedpe_figures_df = pd.DataFrame()
    bedpe_df_cp[['x1','x2','y1','y2']] = bedpe_df_cp[['x1','x2','y1','y2']].astype(int)

    if (type == 'anchor_intersect') or (type == 'anchor_closest'): anchor_list = [['chr1','x1','x2','chr2','y1','y2','sample'],['chr2','y1','y2','chr1','x1','x2','sample']]
    if type == 'within_intersect': anchor_list = [['chr1','x1','y2','chr2','x2','y1','sample']]

    for anchor in anchor_list:
        B_bed_cols_mod = B_bed_cols.copy()
        # bed file has to be sorted by chr and start for closest analysis, so I'll just sort everything
        bedpe_df_cp[anchor].sort_values([anchor[0],anchor[1]]).to_csv(temp_dir + '/all_loops.tsv',
                                                                                   sep='\t',
                                                                                   header=False,
                                                                                   index=False
                                                                                   )

        if (type == 'anchor_intersect') or (type == 'within_intersect'): func = 'bedtools intersect -wao -a ' + temp_dir + '/all_loops.tsv -b ' + gene_file_path + ' > ' + temp_dir + '/loop_genes.tsv'
        if type == 'anchor_closest': func = 'bedtools closest -d -a ' + temp_dir + '/all_loops.tsv -b ' + gene_file_path + ' > ' + temp_dir + '/loop_genes.tsv'
        subprocess.run(func, shell=True)

        bedpe_out = pd.read_csv(temp_dir + '/loop_genes.tsv',
                                sep='\t',
                                header=None
                                )

        if (type == 'anchor_intersect') or (type == 'within_intersect'): extra_col = ['overlap']
        if (type == 'anchor_closest'): extra_col = ['distance']

        bedpe_out.columns = list(anchor) + list(A_bedpe_cols) + list(B_bed_cols_mod) + list(extra_col)
        if (type == 'anchor_closest'): B_bed_cols_mod.append('distance')
        if FPKM_df is not None:
            bedpe_out = bedpe_out.merge(FPKM_df[['gene_id','gene_name','FPKM']],
                                        on=['gene_id','gene_name'],
                                        how='left'
                                        )
            B_bed_cols_mod.append('FPKM')
        bedpe_out[[anchor[1],anchor[2],anchor[4],anchor[5],B_bed_cols_mod[1],B_bed_cols_mod[2]]] = bedpe_out[[anchor[1],anchor[2],anchor[4],anchor[5],B_bed_cols_mod[1],B_bed_cols_mod[2]]].astype('Int64')
        bedpe_out.loc[bedpe_out[B_bed_cols_mod[0]] == '.',B_bed_cols_mod] = pd.NA
        bedpe_out = bedpe_out[list(anchor) + list(A_bedpe_cols) + list(B_bed_cols_mod)]
        bedpe_chr = bedpe_out.copy()
        if type != 'within_intersect': bedpe_chr.loc[bedpe_chr[B_bed_cols_mod[0]].notnull(),'anchor_intersect'] = anchor[0] # mark which anchor the intesect is be drawn from
        bedpe_figures_df = pd.concat([bedpe_figures_df,bedpe_chr])
        
        if (type == 'anchor_intersect') or (type == 'anchor_closest'): gene_cols = [anchor[0] + '_' + i for i in B_bed_cols_mod]
        if (type == 'within_intersect'): gene_cols = B_bed_cols_mod
        bedpe_out.columns = anchor + gene_cols

        if analysis_type == 'dict':
            for sample in annotated_bedpe_dict:
                annotated_bedpe_dict[sample] = annotated_bedpe_dict[sample].merge(bedpe_out,
                                                                                  on=anchor,
                                                                                  how='left'
                                                                                  ).reset_index(drop=True)
                if (anchor == ['chr2','y1','y2','chr1','x1','x2','sample']) or (anchor == ['chr1','x1','y2','chr2','x2','y1','sample']): annotated_bedpe_dict[sample].drop('sample',axis=1,inplace=True)
    

    bedpe_figures_df = sort_by_strings(df=bedpe_figures_df,
                                    order=list(annotated_bedpe_dict.keys()),
                                    sort_col='sample'
                                    )

    if FPKM_df is not None:
        bedpe_figures_grouped = bedpe_figures_df[['chr1','x1','x2','chr2','y1','y2','sample','FPKM']].groupby(by=['chr1','x1','x2','chr2','y1','y2','sample'],as_index=False,sort=False).mean()
        bedpe_figures_df[['chr1','x1','x2','chr2','y1','y2','sample','FPKM','gene_name','gene_id']].to_csv(out_dir + '/loop_FPKM_ungrouped.txt',sep='\t',header=True,index=False)
        bedpe_figures_grouped.to_csv(out_dir + '/loop_FPKM.txt',sep='\t',header=True,index=False)

        data_names,data_list = [],[]
        for sample in bedpe_figures_grouped['sample'].unique():
            data_names.append(sample)
            data_list.append(bedpe_figures_grouped.loc[bedpe_figures_grouped['sample'] == sample,'FPKM'].tolist())
        stats_df = kruskal_wilcoxon(data_names=data_names,data_list=data_list)
        stats_df.to_csv(out_dir + '/loop_FPKM.stats.txt',sep='\t',header=True,index=False)

        plotting.box(melted_df=bedpe_figures_grouped,
            xcol='sample',
            ycol='FPKM',
            out_dir=out_dir,
            measure='chromatin loop mean FPKM',
            title='',
            out_name='FPKM_' + type)

    if analysis_type == 'dict': return annotated_bedpe_dict
