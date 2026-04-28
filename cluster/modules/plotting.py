#!/usr/bin/env python

# this doesn't know how to deal with uneven number of clusters and I should probably fix that at some point
# I want the column widths to be the same regardless of the number of samples and regardless of whether I am subplotting, but I am yet to achieve this

# possible to extract title??

import pandas as pd
import numpy as np
from copy import deepcopy
import seaborn as sns
import matplotlib.pyplot as plt
import math
from matplotlib.ticker import (FormatStrFormatter)
import json
import random
from mpl_toolkits.axes_grid1 import Divider, Size
from PIL import ImageFont, ImageDraw, Image
from matplotlib import font_manager
import os

colors = sns.color_palette()
font_dict = {'big_title':10,'title':8,'labels':6,'axlab':7,'plt-text':5,'legend':6,'fonttype':'arial'}

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'custom_params.json')) as f:
    custom_params = json.loads(f.read())

sns.set_theme(font='arial',style='ticks',rc=custom_params)

def line(melted_df,xcol,ycol,xcol_measure,ycol_measure,out_dir,out_name,title='',hue_col=None,palette=None,subplots=False,order=None,subplot_col=None,ncols=None,matchxy=False,sharex=False,sharey=False,logY=False,add_n=None,Y_range=None,estimator='mean'):
    plot_type = 'sns.lineplot'

    # if you want to use matchxy and sharex/sharey, then both sharex and sharey must be true. I will make this a default
    if matchxy and (sharex or sharey):
        match_share = True
        sharex,sharey = True,True
    else: match_share = False

    #nrows,fig,axs,divider,repeats,order,xmin,xmax,ymin,ymax = plot_setup(plot_type=plot_type,count_matrix=count_matrix,ycol=ycol,xcol=xcol,hue_col=hue_col,order=order,subplot_col=subplot_col,subplots=subplots,ncols=ncols,palette=palette,sharex=False,sharey=sharey)
    nrows,ncols,fig,axs,divider,repeats,order,palette = plot_setup(plot_type=plot_type,
                                                                    melted_df=melted_df,
                                                                    xcol=xcol,
                                                                    ycol=ycol,
                                                                    xcol_measure=xcol_measure,
                                                                    ycol_measure=ycol_measure,
                                                                    hue_col=hue_col,
                                                                    order=order,
                                                                    subplot_col=subplot_col,
                                                                    subplots=subplots,
                                                                    ncols=ncols,
                                                                    palette=palette,
                                                                    sharex=sharex,
                                                                    sharey=sharey)

    rep = 0
    for row in range(nrows):
        for col in range(ncols):
             # set which axes we will use. This has to be done this way because of the arrays are set up.
            fig_ax = ax_selection(nrows=nrows,ncols=ncols,axs=axs,row=row,col=col,divider=divider)
            if rep + 1 > repeats:
                fig.delaxes(fig_ax)
                continue

            if subplots:
                dataset = order[rep]
                data = melted_df[melted_df[subplot_col] == dataset].copy(deep=True)
            else:
                data = melted_df.copy()
                dataset = title

            palette_cp = deepcopy(palette)
            if add_n == 'hue_col':
                if (hue_col is not None):
                    for hue in data[hue_col].unique():
                        if hue == 'non-sig.': continue
                        n_count = data[data[hue_col] == hue].drop([xcol,ycol],axis=1).drop_duplicates().shape[0]
                        #n_count = data[data[hue_col] == hue].shape[0]
                        data.loc[data[hue_col] == hue,hue_col] = hue + ' (n=' + str(n_count) + ')'
                        palette_cp[hue + ' (n=' + str(n_count) + ')'] = palette_cp[hue]
            if (add_n == 'xcol'):
                if data[xcol].dtype == 'object':
                    for x in data[xcol].unique():
                        if x == 'non-sig.': continue
                        n_count = data[data[xcol] == x].shape[0]
                        data.loc[data[xcol] == x,xcol] = x + ' (n=' + str(n_count) + ')'
            if (add_n == 'subplot') or (add_n == 'title'):
                if (hue_col is not None): n_count = data.drop([xcol,ycol,hue_col],axis=1).drop_duplicates().shape[0]
                else: n_count = data.drop([xcol,ycol],axis=1).drop_duplicates().shape[0]
                dataset = dataset + ' (n=' + str(n_count) + ')'
                        
                # n_count = data.shape[0]
                # dataset = dataset + ' (n=' + str(n_count) + ')'

            rep += 1

            if logY:
                if data[data[ycol] == 0].shape[0] > 0:
                    print('data contains zero values in ' + ycol + ' that cannot be converted to log scale')
                    print('please adjust input data before continuing')
                    return
                data[ycol] = np.log10(data[ycol])
                use_ycol_measure = 'log$_{10}$(' + ycol_measure + ')'
            else: use_ycol_measure = ycol_measure
            if estimator == 'median': errorbar = ('ci')
            if estimator == 'mean': errorbar = 'se'
            sns.lineplot(data=data, x=xcol, y=ycol, hue=hue_col, ax=fig_ax, palette=palette_cp,errorbar=errorbar,estimator=estimator)

            if match_share:
                axis_min = min(fig_ax.get_xlim()[0],fig_ax.get_ylim()[0])
                axis_max = max(fig_ax.get_xlim()[1],fig_ax.get_ylim()[1])

                fig_ax.set_xlim(axis_min,axis_max)
                fig_ax.set_ylim(axis_min,axis_max)

            plot_defaults_a(xlabel=xcol_measure,ylabel=use_ycol_measure,fig_ax=fig_ax,title=dataset,hue_col=hue_col)
            if Y_range is None: plot_defaults_b(fig_ax=fig_ax,axis='y')
            else: fig_ax.set_ylim(Y_range)
            if xcol_measure is not None: plot_defaults_b(fig_ax=fig_ax,axis='x')

            if matchxy and (not match_share):
                max_axis = max(fig_ax.get_xlim()[1] - fig_ax.get_xlim()[0],fig_ax.get_ylim()[1] - fig_ax.get_ylim()[0])
                if max_axis == (fig_ax.get_xlim()[1] - fig_ax.get_xlim()[0]):
                    axis_lim = fig_ax.get_xlim()
                    axis_ticks = fig_ax.get_xticks()
                if max_axis ==(fig_ax.get_ylim()[1] - fig_ax.get_ylim()[0]):
                    axis_lim = fig_ax.get_ylim()
                    axis_ticks = fig_ax.get_yticks()
                fig_ax.set_xlim(axis_lim)
                fig_ax.set_ylim(axis_lim)
                #fig_ax.set_xticks(axis_ticks)
                fig_ax.set_yticks(axis_ticks)



    sns.despine()
    #plt.tight_layout()
    plt.savefig(out_dir + '/' + out_name + '.png',dpi=300)
    plt.savefig(out_dir + '/' + out_name + '.pdf',dpi=300)
    plt.figure()
    plt.close(fig)


# if you are not using a hue_col, but want to pass a color to plot with,
    # then palette should be a single color (e.g. palette='#e13661')
def strip(melted_df,xcol,ycol,out_dir,measure,title,out_name,hue_col=None,palette=None,subplots=False,order=None,subplot_col=None,ncols=None,sharey=False,logY=False,label_median=True,add_n=False,Y_range=None):
    plot_type = 'sns.stripplot'

    #nrows,fig,axs,divider,repeats,order,xmin,xmax,ymin,ymax = plot_setup(plot_type=plot_type,count_matrix=count_matrix,ycol=ycol,xcol=xcol,hue_col=hue_col,order=order,subplot_col=subplot_col,subplots=subplots,ncols=ncols,palette=palette,sharex=False,sharey=sharey)
    nrows,ncols,fig,axs,divider,repeats,order,palette = plot_setup(plot_type=plot_type,melted_df=melted_df,ycol=ycol,xcol=xcol,ycol_measure=measure,xcol_measure=None,hue_col=hue_col,order=order,subplot_col=subplot_col,subplots=subplots,ncols=ncols,palette=palette,sharex=False,sharey=sharey)
    rep = 0
    overlay_done =[]
    for row in range(nrows):
        for col in range(ncols):
            # set which axes we will use. This has to be done this way because of the arrays are set up.
            fig_ax = ax_selection(nrows=nrows,ncols=ncols,axs=axs,row=row,col=col,divider=divider)
            if rep + 1 > repeats:
                fig.delaxes(fig_ax)
                continue
            
            if subplots:
                dataset = order[rep]
                data = melted_df[melted_df[subplot_col] == dataset].copy()
            else:
                data = melted_df.copy()
                dataset = title

            palette_cp = deepcopy(palette)
            if add_n == 'hue_col':
                # i'm not super sure this is the best way to do this
                if (hue_col is not None):
                    for hue in data[hue_col].unique():
                        if hue == 'non-sig.': continue
                        n_count = data[data[hue_col] == hue].drop([xcol,ycol],axis=1).drop_duplicates().shape[0]
                        #selse: n_count = data[data[hue_col] == hue].shape[0]
                        data.loc[data[hue_col] == hue,hue_col] = hue + ' (n=' + str(n_count) + ')'
                        palette_cp[hue + ' (n=' + str(n_count) + ')'] = palette_cp[hue]
            if (add_n == 'xcol'):
                for x in data[xcol].unique():
                    if x == 'non-sig.': continue
                    if hue_col is not None: n_count = data[data[xcol] == x].drop([ycol,hue_col],axis=1).drop_duplicates().shape[0]
                    else: n_count = data[data[xcol] == x].shape[0]
                    data.loc[data[xcol] == x,xcol] = x + ' (n=' + str(n_count) + ')'
            if (add_n == 'subplot') or (add_n == 'title'):
                if (hue_col is not None): n_count = data.drop([xcol,ycol,hue_col],axis=1).drop_duplicates().shape[0]
                else: n_count = data.drop([xcol,ycol],axis=1).drop_duplicates().shape[0]
                dataset = dataset + ' (n=' + str(n_count) + ')'

            medians = []
            for cat in data[xcol].unique():
                if hue_col is not None:
                    for hue in data[hue_col].unique():
                        if data.loc[(data[xcol] == cat) & (data[hue_col] == hue)].shape[0] == 0: medians.append(np.nan)
                        else: medians.append(round(data.loc[(data[xcol] == cat) & (data[hue_col] == hue),ycol].median(),2))
                else: medians.append(round(data.loc[data[xcol] == cat,ycol].median(),2))

            if logY:
                if data[data[ycol] == 0].shape[0] > 0:
                    print('data contains zero values in ' + ycol + ' that cannot be converted to log scale')
                    print('please adjust input data before continuing')
                    return
                data[ycol] = np.log10(data[ycol])
                ycol_measure = 'log$_{10}$(' + measure + ')'
            else: ycol_measure = measure

            ypos = []
            for cat in data[xcol].unique():
                if hue_col is not None:
                    for hue in data[hue_col].unique():
                        if data.loc[(data[xcol] == cat) & (data[hue_col] == hue)].shape[0] == 0: ypos.append(np.nan)
                        else: ypos.append((data.loc[(data[xcol] == cat) & (data[hue_col] == hue),ycol].median()))
                else: ypos.append(round(data.loc[data[xcol] == cat,ycol].median(),4))

            rep += 1

            if (hue_col is None) and (palette_cp is not None): color,palette_cp = palette_cp,None
            else: color = None

            sns.stripplot(data=data,
                          x=xcol,
                          y=ycol,
                          hue=hue_col,
                          size=1,
                          linewidth=0,
                          ax=fig_ax,
                          dodge=True,
                          palette=palette_cp,
                          color=color
                          )
            if hue_col is not None: h, l = fig_ax.get_legend_handles_labels()

            # plot the median line and add the median values
            box_plot = sns.boxplot(medianprops={'visible': True,'color': 'black', 'ls': '-', 'lw': 0.5},
                                whiskerprops={'visible': False},
                                zorder=10,
                                hue=hue_col,
                                data=data,
                                x=xcol,
                                y=ycol,
                                showfliers=False,
                                showbox=False,
                                showcaps=False,
                                ax=fig_ax)

            vertical_offset = int(fig_ax.get_ylim()[1])/40
            if subplots or (hue_col is not None): horizontal_offset = 0.05
            else: horizontal_offset = 0.25

            if label_median: 
                x_coords = []
                for c in box_plot.collections:
                    offsets = c.get_offsets()
                    if len(offsets) != 0:
                            xs = [x[0] for x in offsets]
                            x_coords.append(np.nanmean(xs))
                xtick_rep = 0
                for xtick in x_coords:
                    if math.isnan(ypos[xtick_rep]):
                        xtick_rep += 1
                        continue
                    fig_ax.text(xtick - horizontal_offset,ypos[xtick_rep] + vertical_offset,medians[xtick_rep],
                                ha='center',va='bottom',size=font_dict['plt-text'],color='black',rotation='vertical') 
                    xtick_rep += 1

            if Y_range is None: plot_defaults_b(fig_ax=fig_ax,axis='y')
            else: fig_ax.set_ylim(Y_range)

            plot_defaults_a(xlabel=None,ylabel=ycol_measure,fig_ax=fig_ax,title=dataset,hue_col=hue_col)
            if hue_col is not None:
                fig_ax.get_legend().remove()
                legend = fig_ax.legend(h,l,
                                    title='',
                                    markerscale=2,
                                    handletextpad=-0.4,
                                    scatteryoffsets=[0.5]
                                    )

    if sharey: plot_defaults_b(fig_ax=fig_ax,axis='y')
    
    sns.despine()
    plt.savefig(out_dir + '/' + out_name + '.png',dpi=300)
    plt.savefig(out_dir + '/' + out_name + '.pdf',dpi=300)
    plt.figure()
    plt.close('all')


def box(melted_df,xcol,ycol,out_dir,measure,title,out_name,hue_col=None,palette=None,subplots=False,order=None,subplot_col=None,ncols=None,sharey=False,logY=False,add_n=None,Y_range=None,showmeans=False):
    plot_type = 'sns.boxplot'

    nrows,ncols,fig,axs,divider,repeats,order,palette = plot_setup(plot_type=plot_type,melted_df=melted_df,ycol=ycol,xcol=xcol,ycol_measure=measure,xcol_measure=None,hue_col=hue_col,order=order,subplot_col=subplot_col,subplots=subplots,ncols=ncols,palette=palette,sharex=False,sharey=sharey)

    rep = 0
    for row in range(nrows):
        for col in range(ncols):
            # set which axes we will use. This has to be done this way because of the arrays are set up.
            fig_ax = ax_selection(nrows=nrows,ncols=ncols,axs=axs,row=row,col=col,divider=divider)

            if rep + 1 > repeats:
                fig.delaxes(fig_ax)
                continue

            if subplots:
                dataset = order[rep]
                data = melted_df[melted_df[subplot_col]==dataset].copy(deep=True)
            else:
                data = melted_df.copy(deep=True)
                dataset = title
            
            palette_cp = deepcopy(palette)
            if add_n == 'hue_col':
                # i'm not super sure this is the best way to do this
                if (hue_col is not None):
                    for hue in data[hue_col].unique():
                        if hue == 'non-sig.': continue
                        if data[xcol].dtype == 'object':
                            n_count = data[data[hue_col] == hue].drop([xcol,ycol],axis=1).drop_duplicates().shape[0]
                        else: n_count = data[data[hue_col] == hue].shape[0]
                        data.loc[data[hue_col] == hue,hue_col] = hue + ' (n=' + str(n_count) + ')'
                        palette_cp[hue + ' (n=' + str(n_count) + ')'] = palette_cp[hue]
            if (add_n == 'xcol'):
                if data[xcol].dtype == 'object':
                    for x in data[xcol].unique():
                        if x == 'non-sig.': continue
                        n_count = data[data[xcol] == x].shape[0]
                        data.loc[data[xcol] == x,xcol] = x + ' (n=' + str(n_count) + ')'
            if (add_n == 'subplot') or (add_n == 'title'):
                if (hue_col is not None): n_count = data.drop([xcol,ycol,hue_col],axis=1).drop_duplicates().shape[0]
                else: n_count = data.drop([xcol,ycol],axis=1).drop_duplicates().shape[0]
                dataset = dataset + ' (n=' + str(n_count) + ')'

            rep += 1

            if logY:
                if data[data[ycol] == 0].shape[0] > 0:
                    print('data contains zero values in ' + ycol + ' that cannot be converted to log scale')
                    print('please adjust input data before continuing')
                    return
                data[ycol] = np.log10(data[ycol])
                ycol_measure = 'log$_{10}$(' + measure + ')'
            else: ycol_measure = measure

            try: sns.boxplot(data=data,
                             x=xcol,
                             y=ycol,
                             hue=hue_col,
                             ax=fig_ax,
                             dodge=True,
                             palette=palette_cp,
                             flierprops={'markersize': 1},
                             showfliers=False,
                             linewidth=0.5,
                             linecolor='black',
                             showmeans=showmeans
                             )
            except: sns.boxplot(data=data,
                                x=xcol,
                                y=ycol,
                                hue=hue_col,
                                ax=fig_ax,
                                dodge=True,
                                palette=palette_cp,
                                flierprops={'markersize': 1},
                                showfliers=False,
                                linewidth=0.5,
                                showmeans=showmeans
                                )
            

            plot_defaults_a(xlabel=None,ylabel=ycol_measure,fig_ax=fig_ax,title=dataset,hue_col=hue_col)
            
            if Y_range is None: plot_defaults_b(fig_ax=fig_ax,axis='y')
            else: fig_ax.set_ylim(Y_range)



    sns.despine()
    plt.savefig(out_dir + '/' + out_name + '.png',dpi=300)
    plt.savefig(out_dir + '/' + out_name + '.pdf',dpi=300)
    plt.figure()
    plt.close(fig)

def heat(count_matrix,data_cols,out_dir,measure,data_type,title,out_name,subplot_col=None,subplots=False,proportional=False,order=None,add_n=False,v_range=None,cmap=None):
    # determine the color palette to use based on the range of the data
    if cmap is None:
        min_out = min(count_matrix[data_cols].min())
        if min_out < 0: cmap = sns.color_palette('vlag', as_cmap=True)
        else: cmap = sns.light_palette('tomato', as_cmap=True)

    if v_range is not None: vmin,vmax = v_range[0],v_range[1]
    else: vmin,vmax = None,None

    ax = sns.heatmap(count_matrix[data_cols], cmap=cmap)
    im = ax.imshow(count_matrix[data_cols],cmap=cmap,vmin=vmin,vmax=vmax)
    plt.figure()

    height_ratio = None
    if subplots and proportional:
        if order == None: order = count_matrix[subplot_col].unique()
        prop_info = pd.DataFrame()
        for group in order:
            prop_info.loc[group,subplot_col] = group
            prop_info.loc[group,'proportion'] = count_matrix[count_matrix[subplot_col] == group].shape[0]/count_matrix.shape[0]
        prop_info.reset_index(drop=True,inplace=True)
        ratio = prop_info['proportion'].tolist()
        height_ratio = ratio + [np.mean(ratio),np.mean(ratio)]
        height_ratio = [i/max(height_ratio) for i in height_ratio]
    nrows,fig,axs,divider,repeats,order = heatmap_plot_setup(count_matrix=count_matrix.copy(),data_cols=data_cols,measure=measure,order=order,subplot_col=subplot_col,subplots=subplots,sharex=True,add_n=add_n,height_ratio=height_ratio)

    rep = 0
    ylabel_list = []
    for row in range(nrows):
        if (subplots):
            if row >= nrows - 2:
                fig_ax = axs[row,0]
                fig.delaxes(fig_ax)
                continue

            fig_ax = axs[row,0]
            fig_ax.set_axes_locator(divider.new_locator(nx=1, ny=2 * nrows - 2 * row - 1))

            cluster = order[rep]
            data = count_matrix[count_matrix[subplot_col]==cluster].copy()
            title = cluster

        else:
            fig_ax = axs[0]
            fig_ax.set_axes_locator(divider.new_locator(nx=1, ny=1))  
            data = count_matrix

        rep += 1

        sns.heatmap(data[data_cols], ax=fig_ax, cmap=cmap, yticklabels = False, cbar=False,vmin=vmin,vmax=vmax)
        
        if not subplots: fig_ax.set_ylabel(data_type + ' (n=' + str(data.shape[0]) + ')', va='center',fontsize = font_dict['axlab'],color='black',labelpad=15)
        if subplots:
            if add_n: y_label = cluster + ' (n=' + str(data.shape[0]) + ')'
            else: y_label = cluster
            ylabel_list.append(y_label)
            #add cluster information to each subplot as a ylab
            fig_ax.set_ylabel(y_label, labelpad=-1, rotation=0,va='center',ha='right',fontsize = font_dict['plt-text'],color='black')
            # this removes the ticks from the individual heatmaps except the bottom one
            if rep != len(order): fig_ax.tick_params(axis='x',bottom=False,labelbottom=False)

        # add a border around each individual heatmap
        for _, spine in fig_ax.spines.items():
            spine.set_visible(True)
            spine.set_linewidth(0.25)
            spine.set_color('black')

        fig_ax.tick_params(axis='x',labelrotation=90)
        fig_ax.set_xlabel('')

    for row in range(nrows):
        #if row != nrows - 1: continue
        cb_ax = ax_selection(nrows=nrows,ncols=2,axs=axs,row=row,col=1,divider=divider)
        if nrows > 2: cb_row = int(math.ceil((nrows - 2)/2))
        else: cb_row = 0
        if row == cb_row:
            cbar = fig.colorbar(im,orientation='vertical',cax=cb_ax)
            cbar.set_label(measure,labelpad=-30)
            cbar.outline.set_color('black')
            cbar.outline.set_linewidth(0.25)
            cbar.ax.tick_params(width=0.25,color='black',labelcolor='black')

            cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            cbar_ticks = []
            for i in cbar.get_ticks(): cbar_ticks.append(str(round(i,1)))
        else: fig.delaxes(cb_ax)


    plt.savefig(out_dir + '/' + out_name + '.png',dpi=300)
    plt.savefig(out_dir + '/' + out_name + '.pdf')
    plt.figure()
    plt.close(fig)


def joint(melted_df,ycol,xcol,ycol_measure,xcol_measure,out_dir,out_name,title,hue_col=None,palette=None,ncols=1,logX=False,logY=False,modify_labels=False,add_n=False,X_range=None,Y_range=None):
    plot_type = 'sns.scatterplot'

    nrows,ncols,fig,axs,divider,repeats,order,palette = plot_setup(plot_type=plot_type,melted_df=melted_df,ycol=ycol,xcol=xcol,ycol_measure=ycol_measure,xcol_measure=xcol_measure,hue_col=hue_col,order=None,subplot_col=None,subplots=False,ncols=1,palette=palette,sharex=True,sharey=True)

    data = melted_df.copy()

    palette_cp = deepcopy(palette)
    if add_n == 'hue_col':
        if (hue_col is not None):
            for hue in data[hue_col].unique():
                #if hue == 'non-sig.': continue
                n_count = data[data[hue_col] == hue].shape[0]
                if n_count == 0: continue
                data.loc[data[hue_col] == hue,hue_col] = hue + ' (n=' + str(n_count) + ')'
                palette_cp[hue + ' (n=' + str(n_count) + ')'] = palette_cp[hue]


    if logX:
        if data[data[xcol] == 0].shape[0] > 0:
            print('data contains zero values in ' + xcol + ' that cannot be converted to log scale')
            print('please adjust input data before continuing')
            return
        data[xcol] = np.log10(data[xcol])
        use_xcol_measure = 'log$_{10}$(' + xcol_measure + ')'
    else: use_xcol_measure = xcol_measure

    if logY:
        if data[data[ycol] == 0].shape[0] > 0:
            print('data contains zero values in ' + ycol + ' that cannot be converted to log scale')
            print('please adjust input data before continuing')
            return
        data[ycol] = np.log10(data[ycol])
        use_ycol_measure = 'log$_{10}$(' + ycol_measure + ')'
    else: use_ycol_measure = ycol_measure

    fig_ax = sns.jointplot(data=data,
                            x=xcol,
                            y=ycol,
                            hue=hue_col,
                            height=2.5,
                            space=0.2,
                            ratio=3,
                            palette=palette_cp,
                            joint_kws={'linewidth':0,'s':3},
                            marginal_kws={'linewidths':0.5,'common_norm':False})

    if modify_labels: xlabel,ylabel = title.split('_vs_')[0] + '\n' + use_xcol_measure,title.split('_vs_')[1] + '\n' + use_ycol_measure
    else: xlabel,ylabel = use_xcol_measure,use_ycol_measure
    if ('_vs_' in title) and (len(title) > 20): subplot_title = title.split('_vs_')[0] + ' vs\n' + title.split('_vs_')[1]
    else: subplot_title = title

    #plot_defaults_a(xlabel=None,ylabel=None,fig_ax=fig_ax,title=None,hue_col=hue_col)
    fig_ax.set_axis_labels(xlabel, ylabel)

    axis_min = min(fig_ax.ax_joint.get_xlim()[0],fig_ax.ax_joint.get_ylim()[0])
    axis_max = max(fig_ax.ax_joint.get_xlim()[1],fig_ax.ax_joint.get_ylim()[1])

    plot_defaults_b(fig_ax=fig_ax.ax_joint,axis='x')
    plot_defaults_b(fig_ax=fig_ax.ax_joint,axis='y')

    xticks = fig_ax.ax_joint.get_xticks()
    yticks = fig_ax.ax_joint.get_yticks()

    if list(xticks) != list(yticks):
        if xticks[len(xticks) - 1] == fig_ax.ax_joint.get_xlim()[1]: fig_ax.ax_joint.set_yticks(xticks)
    elif yticks[len(yticks) - 1] == fig_ax.ax_joint.get_ylim()[1]: fig_ax.ax_joint.set_xticks(yticks)

    fig_ax.ax_joint.legend(title='',loc='best',markerscale=1,handletextpad=0.1)
    fig_ax.ax_marg_x.set_title(subplot_title)

    if X_range is None: fig_ax.ax_joint.set_xlim(-4,0)
    else: fig_ax.ax_joint.set_xlim(X_range[0],X_range[1])
    
    if Y_range is None: fig_ax.ax_joint.set_ylim(-4,0)
    else: fig_ax.ax_joint.set_ylim(Y_range[0],Y_range[1])

    sns.despine()
    plt.savefig(out_dir + '/' + out_name + '.png',dpi=300)
    plt.savefig(out_dir + '/' + out_name + '.pdf',dpi=300)
    plt.figure()

def stacked(measure,title,out_dir,out_name,count_table=None,count_dict=None,labels=False,palette=None,order=None,sem_table=None,sem_dict=None,Y_range=None):
    if (count_table is None) and (count_dict is None): print('must provide either count_table or count_dict')

    if count_dict is not None:
        subplots = True
        subplot_col = 'subplot_col'
        count_table_submit = pd.DataFrame()
        for dataset in count_dict:
            df = pd.melt(count_dict[dataset].reset_index(drop=False).rename(columns={'index':'xcol'}),id_vars='xcol',value_name='ycol',var_name='stack')
            df['subplot_col'] = dataset
            count_table_submit = pd.concat([count_table_submit,df])
        
    elif (count_table is not None):
        if (count_dict is not None): print('you provided both a count_table and count_dict, will default to using only the count_table')
        subplots = False
        subplot_col = None
        count_table_submit = pd.melt(count_table.reset_index(drop=False).rename(columns={'index':'xcol'}),id_vars='xcol',value_name='ycol',var_name='stack')

    plot_type = 'pd.plot'
    nrows,ncols,fig,axs,divider,repeats,order,palette = plot_setup(plot_type=plot_type,
                                                                   melted_df=count_table_submit,
                                                                   ycol='ycol',
                                                                   xcol='xcol',
                                                                   ycol_measure=measure,
                                                                   xcol_measure=None,
                                                                   hue_col='stack',
                                                                   order=order,
                                                                   subplot_col=subplot_col,
                                                                   subplots=subplots,
                                                                   palette=palette,
                                                                   sharex=False,
                                                                   sharey=False
                                                                   )

    rep = 0
    for row in range(nrows):
        for col in range(ncols):
            fig_ax = ax_selection(nrows=nrows,ncols=ncols,axs=axs,row=row,col=col,divider=divider)
            if rep + 1 > repeats:
                fig.delaxes(fig_ax)
                continue
            
            sem = None
            if subplots:
                dataset = order[rep]
                data = count_dict[dataset]
                if sem_dict is not None: sem = sem_dict[dataset]
            else:
                dataset = title
                data = count_table
                if sem_table is not None: sem = sem_table

            rep += 1

            data.plot(kind='bar', stacked=True,color=palette,ax=fig_ax,yerr=sem)
            plot_defaults_a(xlabel='',ylabel=measure,fig_ax=fig_ax,title=dataset,hue_col=True)
            fig_ax.legend(bbox_to_anchor = (0.1,-0.05))
            plot_defaults_b(fig_ax=fig_ax,axis='y')
            
            if Y_range is not None: fig_ax.set_ylim([Y_range[0],Y_range[1]])
            else: fig_ax.set_ylim([0,100])

            if labels:
                done = []
                for enhancer in list(data.columns):
                    done.append(enhancer)
                    for xtick in fig_ax.get_xticks():
                        if pd.isnull(data.loc[data.index[xtick],enhancer]) or round(data.loc[data.index[xtick],enhancer],1) <= 0.0: continue
                        vo = data.loc[data.index[xtick],done].sum() - 0.5 * data.loc[data.index[xtick],enhancer]
                        fig_ax.text(xtick,vo,round(data.loc[data.index[xtick],enhancer],1),
                                    ha='center',va='center',size=font_dict['plt-text'],color='black',rotation='horizontal')
    sns.despine()
    fig.savefig(out_dir + '/' + out_name + '.png',dpi=300)
    fig.savefig(out_dir + '/' + out_name + '.pdf',dpi=300)
    plt.figure()

def bar(melted_df,xcol,ycol,out_dir,measure,title,out_name,hue_col=None,palette=None,subplots=False,order=None,subplot_col=None,ncols=None,sharey=False,overlay_col=None,Y_range=None):
    plot_type = 'sns.barplot'

    nrows,ncols,fig,axs,divider,repeats,order,palette = plot_setup(plot_type=plot_type,melted_df=melted_df,ycol=ycol,xcol=xcol,ycol_measure=measure,xcol_measure=None,hue_col=hue_col,order=order,subplot_col=subplot_col,subplots=subplots,ncols=ncols,palette=palette,sharex=False,sharey=sharey)

    rep = 0
    for row in range(nrows):
        for col in range(ncols):
            # set which axes we will use. This has to be done this way because of the arrays are set up.
            fig_ax = ax_selection(nrows=nrows,ncols=ncols,axs=axs,row=row,col=col,divider=divider)

            if rep + 1 > repeats:
                fig.delaxes(fig_ax)
                continue

            if subplots:
                dataset = order[rep]
                data = melted_df[melted_df[subplot_col]==dataset]
            else:             
                data = melted_df
                dataset = title

            rep += 1
            sns.barplot(data=data, x=xcol, y=ycol, hue=hue_col, ax=fig_ax, palette=palette,edgecolor='black',linewidth=0.25,errorbar='se')
            h, l = fig_ax.get_legend_handles_labels()
            
            if overlay_col:
                sns.boxplot(showmeans=True,
                            meanline=True,
                            meanprops={'color': 'black', 'ls': 'dotted', 'lw': 1},
                            medianprops={'visible': False},
                            whiskerprops={'visible': False},
                            zorder=10,
                            data=data,
                            x=xcol,
                            y=overlay_col,
                            hue=hue_col,
                            showfliers=False,
                            showbox=False,
                            showcaps=False,
                            ax=fig_ax)

            
            
            plot_defaults_a(xlabel=None,ylabel=measure,fig_ax=fig_ax,title=dataset,hue_col=hue_col)
            try:
                fig_ax.get_legend().remove()
                fig_ax.legend(h,l,title='')
            except AttributeError:
                pass
            if (not sharey) and (Y_range is None): plot_defaults_b(fig_ax=fig_ax,axis='y')

            if Y_range is not None: fig_ax.set_ylim(Y_range[0],Y_range[1])
    if sharey: plot_defaults_b(fig_ax=fig_ax,axis='y')

    sns.despine()
    plt.savefig(out_dir + '/' + out_name + '.png',dpi=300)
    plt.savefig(out_dir + '/' + out_name + '.pdf',dpi=300)
    plt.figure()

def plot_setup(plot_type,melted_df,ycol,xcol,ycol_measure,xcol_measure,hue_col,order,subplot_col,subplots,palette,sharex,sharey,ncols=None):
    sns.set_theme(font='arial',style='ticks',rc=custom_params)

    if (hue_col is not None) and (palette is None):
        palette = {}
        for i in melted_df[hue_col].unique(): palette[i] = randomize_hex()
    if (melted_df[xcol].dtypes == 'object'):
        categories = melted_df[xcol].unique()
        xcat_width = []
        for cat in categories: xcat_width.append(text_inches(string=cat)[0])
        xcat_width = max(xcat_width) * 3
        
        if subplots and (subplot_col is not None): 
            xcat_length = []
            for subplot in melted_df[subplot_col].unique(): xcat_length.append(len(melted_df.loc[(melted_df[subplot_col] == subplot),xcol].unique()))
            xcat_length = max(xcat_length)
            if xcat_length > 1: xcat_length = xcat_length
        else: xcat_length = len(categories)
    else: 
        xcat_width = 0.5
        xcat_length = 10

    if (ycol is not None) and (melted_df[ycol].dtypes == 'object'):
        categories = melted_df[ycol].unique()
        ycat_width = []
        for cat in categories: ycat_width.append(text_inches(string=cat)[0])
        ycat_width = max(ycat_width) * 1.2
        ycat_length = len(categories)/3
    else: 
        ycat_width = 0.5
        ycat_length = 10

    if subplots:
        subplot_titles = melted_df[subplot_col].unique()
        if order is None: order = melted_df[subplot_col].unique()
        if ncols is None: nrows,ncols = factor_int(len(order))
        else: nrows = math.ceil(len(order)/ncols)
        subplot_titles = [i.split('_vs_')[0] if (('_vs_' in i) & (len(i) > 20)) else i for i in subplot_titles]

        repeats = len(order)
        letters_title_width = (max([text_inches(string=str(i))[0] for i in subplot_titles]))
        letters_title_height = (max([text_inches(string=str(i))[1] for i in subplot_titles])) * 2
    else:
        ncols,nrows = 1,1
        repeats = 1
        letters_title_width = 0
        letters_title_height = 0

    if ycol_measure is not None:
        ylab_height = (text_inches(string=ycol_measure)[1])
        if '\n' in ycol_measure:
            ylab_height = ylab_height * len(ycol_measure.split('\n'))
            ycol_measure = ycol_measure.split('\n')[0]
        ylab_width = (text_inches(string=ycol_measure)[0])
        
    else: ylab_width,ylab_height = 0,0

    ax_height = 1

    hue_multiplier = 1
    if ((plot_type == 'sns.stripplot') or (plot_type == 'sns.boxplot') or (plot_type == 'sns.barplot')):
        ax_height,width_multiplier = 1.2,1.4
        if (hue_col is not None): 
            if subplots and (subplot_col is not None): 
                hue_count = []
                for subplot in melted_df[subplot_col].unique(): hue_count.append(len(melted_df.loc[(melted_df[subplot_col] == subplot),hue_col].unique()))
                hue_count = max(hue_count)
                if hue_count > 1: hue_multiplier = hue_count * 2
            else: hue_multiplier = len(melted_df[hue_col].unique()) * 2
    if plot_type == 'sns.lineplot': width_multiplier = 1.5

    ax_width = 0.1 * xcat_length * hue_multiplier
    if (ycol is not None) and (melted_df[ycol].dtypes == 'object') and (not subplots): ax_height = 0.5 * ycat_length * hue_multiplier
    if (ycol is not None) and (melted_df[ycol].dtypes == 'object') and subplots: ax_height = 0.12 * ycat_length * hue_multiplier

    if ax_width < 1: ax_width = 1
    if  (plot_type == 'sns.scatterplot') and (melted_df[xcol].dtypes != 'object') and (melted_df[ycol].dtypes != 'object'): ax_height = ax_width

    fig_height = (max(ax_height + xcat_width,ylab_width) * 2 + letters_title_height) * nrows
    fig_width = (max(ax_width + ycat_width,letters_title_width) + ylab_height) * (ncols + 1)

    fig, axs = plt.subplots(nrows=nrows, ncols=ncols,figsize=(fig_width,fig_height),sharex=sharex,sharey=sharey)

    sc = Size.Scaled(2)
    fh = Size.Fixed(ax_width)
    fv = Size.Fixed(ax_height)

    h = [sc, fh] * ncols + [sc]
    v = [sc, fv] * nrows + [sc]
    divider = Divider(fig, (0.0, 0.0, 1., 1.), h, v)

    return nrows,ncols,fig,axs,divider,repeats,order,palette

def heatmap_plot_setup(count_matrix,data_cols,measure,order,subplot_col,subplots,sharex,add_n,height_ratio=None):
    sns.set_theme(font='arial',style='ticks',rc=custom_params)

    if subplots:
        subplot_titles = count_matrix[subplot_col].unique()
        if add_n: [i + ' (n' + str(count_matrix[count_matrix[subplot_col] == i].shape[0]) + ')' for i in subplot_titles]
        if order is None: order = count_matrix[subplot_col].unique()
        nrows = len(order) + 2
        repeats = len(order)
    else:
        nrows = 1
        repeats = 1

    ylab_width = (text_inches(string=measure)[0] * font_dict['axlab']/1000)/72
    ylab_height = (text_inches(string=measure)[1] * font_dict['axlab']/1000)/72

    if subplots:
        letters_title_width = (max([text_inches(string=str(i))[0] for i in subplot_titles]))
        letters_title_height = (max([text_inches(string=str(i))[1] for i in subplot_titles]))
    else:
        letters_title_width = 0
        letters_title_height = 0

    ax_height = 0.25

    width_multiplier = 0.9
    height_multiplier = 1.5

    if not subplots: height_multiplier,ax_height = 4,1

    ax_width = 0.2 * len(data_cols)

    #fig_height = (max(ax_height,ylab_width) + letters_title_height) * height_multiplier * nrows * 2
    fig_height = (ax_height + letters_title_height) * nrows * height_multiplier
    fig_width = (max(ax_width,letters_title_width) + ylab_height) * width_multiplier  * 2

    fig_width = fig_width + (1.25/fig_width)**3

    fig, axs = plt.subplots(nrows=nrows, ncols=2,figsize=(fig_width,fig_height),sharex=False)

    sc = Size.Scaled(2)
    fh = Size.Fixed(ax_width)
    fv = Size.Fixed(ax_height)

    h = [Size.Fixed(0.5), Size.Fixed(ax_width), Size.Fixed(0.2), Size.Fixed(0.1),Size.Fixed(0.1)]

    if height_ratio is not None:
        v = []
        ax_height = ax_height * 2
        for height in height_ratio[::-1]:
            v = v + [sc,Size.Fixed(ax_height*height)]
        v = v + [sc]
    else: v = [sc, fv] * nrows + [sc]

    divider = Divider(fig, (0.0, 0.0, 1., 1.), h, v)

    return nrows,fig,axs,divider,repeats,order


def plot_defaults_a(xlabel,ylabel,fig_ax,title,hue_col=None):
    fig_ax.set_title(title)

    fig_ax.xaxis.set_tick_params(labelbottom=True)
    fig_ax.yaxis.set_tick_params(labelleft=True)
    fig_ax.tick_params(axis='x',labelrotation=90)
    if hue_col is not None: fig_ax.legend(title='',loc='best',markerscale=1,handletextpad=0.1)

    fig_ax.set_xlabel(xlabel)
    fig_ax.set_ylabel(ylabel)
    fig_ax.xaxis.get_label().set_visible(True)
    fig_ax.yaxis.get_label().set_visible(True)

def plot_defaults_b(fig_ax,axis):
    if axis == 'x':
        axis_min,axis_max = fig_ax.get_xticks()[1],fig_ax.get_xticks()[len(fig_ax.get_xticks()) - 2]
        if fig_ax.get_xticks()[len(fig_ax.get_xticks()) - 2] < fig_ax.get_xlim()[1]:
            axis_max = fig_ax.get_xticks()[len(fig_ax.get_xticks()) - 2] + (fig_ax.get_xticks()[len(fig_ax.get_xticks()) - 2] - fig_ax.get_xticks()[len(fig_ax.get_xticks()) - 3])
        if axis_min < 0: min_val,max_val = (-max([abs(math.floor(axis_min)),abs(math.ceil(axis_max))])),(max([abs(math.floor(axis_min)),abs(math.ceil(axis_max))]))
        elif axis_max <= 1: min_val,max_val = 0,axis_max
        else: min_val,max_val = 0,math.ceil(axis_max)
        fig_ax.set_xlim(min_val,max_val)
        fig_ax.set_xticks(fig_ax.get_xticks())

    if axis == 'y':
        axis_min,axis_max = fig_ax.get_yticks()[1],fig_ax.get_yticks()[len(fig_ax.get_yticks()) - 2]
        if fig_ax.get_yticks()[len(fig_ax.get_yticks()) - 2] < fig_ax.get_ylim()[1]:
            axis_max = fig_ax.get_yticks()[len(fig_ax.get_yticks()) - 2] + (fig_ax.get_yticks()[len(fig_ax.get_yticks()) - 2] - fig_ax.get_yticks()[len(fig_ax.get_yticks()) - 3])
        if axis_min < 0: min_val,max_val = (-max([abs(math.floor(axis_min)),abs(math.ceil(axis_max))])),(max([abs(math.floor(axis_min)),abs(math.ceil(axis_max))]))
        elif axis_max <= 1: min_val,max_val = 0,axis_max
        else: min_val,max_val = 0,math.ceil(axis_max)
        fig_ax.set_ylim(min_val,max_val)
        fig_ax.set_yticks(fig_ax.get_yticks())

def ax_selection(nrows,ncols,axs,row,col,divider,subplots=True):
    # set which axes we will use. This has to be done this way because of the arrays are set up.
    if (nrows == 1) and (ncols == 1): subplots = False
    if (subplots) and (nrows > 1) and (ncols > 1):
        fig_ax = axs[row,col]
        fig_ax.set_axes_locator(divider.new_locator(nx=col * 2 + 1, ny=2 * nrows - 2 * row - 1))    
    elif (subplots) and ((nrows == 1) or (ncols == 1)):
        if nrows == 1: rep = col
        if ncols == 1: rep = row
        fig_ax = axs[rep]
        fig_ax.set_axes_locator(divider.new_locator(nx=col * 2 + 1, ny=2 * nrows - 2 * row - 1))  
    elif not subplots:
        fig_ax = axs
        #fig_ax.set_axes_locator(divider.new_locator(nx=0, ny=0))  
        fig_ax.set_axes_locator(divider.new_locator(nx=1, ny=1))  

    return fig_ax

def text_inches(string, font_size_points=font_dict['labels'], dpi=72):

    file = font_manager.findfont('arial')

    # Load the font
    font = ImageFont.truetype(file, font_size_points)

    # Create a dummy image and draw object to measure text
    # The actual image content doesn't matter for measurement
    dummy_image = Image.new("RGB", (1, 1))
    draw = ImageDraw.Draw(dummy_image)

    # Get the bounding box of the text (left, top, right, bottom)
    # Using anchor='lt' (left, top) for consistent measurement
    bbox = draw.textbbox((0, 0), string, font=font, anchor='lt')

    # Calculate width in pixels
    width_pixels = bbox[2] - bbox[0]

    # Convert pixels to inches
    width_inches = width_pixels / dpi

    text_height_pixels = bbox[3] - bbox[1]

    # Convert pixels to inches
    text_height_inches = text_height_pixels / dpi

    return width_inches,text_height_inches

def randomize_hex():
    choice = random. randint(0, 9)
    pal = sns.color_palette('muted')
    color_pick = pal.as_hex()[choice]

    return color_pick

def factor_int(n):
    val = math.ceil(math.sqrt(n))
    val2 = int(n/val)
    while val2 * val != float(n):
        val -= 1
        val2 = int(n/val)
    return val, val2
