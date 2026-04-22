import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import json
from mpl_toolkits.axes_grid1 import Divider, Size
import os
from pathlib import Path

with open('/Users/tessapopay/example_data/custom_params.json') as f:
    custom_params = json.loads(f.read())

sns.set_theme(font='arial',style='ticks',rc=custom_params)


def heatmap_plot(path,normalize=False,cmap_color='Blues',vmin=None,vmax=None):

    df = pd.read_csv(path,sep='\t',header=0,index_col=0)
    print(df)
    if 'Base' in df.index: df.drop(index='Base',inplace=True)
    if 'Genome %' in df.columns:
        cols = list(df.columns)
        cols.remove('Genome %')
        df.drop(columns=cols[len(cols) - 1],inplace=True)
        print(df)
        df.columns = cols
        
    #print(df)

    if normalize: df = 100 * df/df.sum()

    fig,axs,divider = plot_setup(df=df)

    cmap = sns.color_palette(cmap_color,as_cmap=True)

    ax = axs[0]
    ax.set_axes_locator(divider.new_locator(nx=1, ny=1))  

    sns.heatmap(df,square=True, cmap=cmap,ax=ax,cbar=False,vmin=vmin,vmax=vmax)

    im = ax.imshow(df,cmap=cmap,vmin=vmin,vmax=vmax)

    cb_ax = axs[1]
    cb_ax.set_axes_locator(divider.new_locator(nx=3, ny=1))  

    cbar = fig.colorbar(im,orientation='vertical',cax=cb_ax)
    #cbar.set_label('log10obs/exp',labelpad=-30)
    cbar.outline.set_color('black')
    cbar.outline.set_linewidth(0.25)
    cbar.ax.tick_params(width=0.25,color='black',labelcolor='black')

    out_path = os.path.dirname(path) + '/' + Path(path).stem
    plt.savefig(out_path + '.png',dpi=300)
    plt.savefig(out_path + '.pdf',dpi=300)
    plt.close()

def plot_setup(df):
    sns.set_theme(font='arial',style='ticks',rc=custom_params)

    ax_height = 0.1 * df.shape[0]
    ax_width = 0.1 * df.shape[1]

    #fig_height = (max(ax_height,ylab_width) + letters_title_height) * height_multiplier * nrows * 2
    fig_height = ax_height + 1.5
    fig_width = ax_width + 3

    #if not subplots: fig_height = fig_height * 2

    fig, axs = plt.subplots(nrows=1,ncols=2,figsize=(fig_width,fig_height))

    sc = Size.Scaled(2)
    fh = Size.Fixed(ax_width)
    fv = Size.Fixed(ax_height)

    h = [sc, fh, Size.Fixed(0.05),Size.Fixed(0.1),Size.Scaled(1)]
    v = [sc, fv, sc]
    divider = Divider(fig, (0.0, 0.0, 1., 1.), h, v)

    return fig,axs,divider