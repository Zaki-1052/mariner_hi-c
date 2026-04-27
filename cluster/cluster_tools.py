
import pandas as pd

from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

def elbow(out_dir,count_matrix,data_cols):
    # create elbox plot using range of kmeans (1-15)
    sse = []
    for k in range(1, 15):
        km = KMeans(n_clusters=k,n_init=10)
        km.fit(count_matrix[data_cols])
        sse.append(km.inertia_)
    plt.figure(figsize=(6, 6))
    plt.plot(range(1, 15), sse, '-o', c = 'maroon')
    plt.xlabel('Count of Clusters')
    plt.ylabel('SSE')
    plt.savefig(out_dir + '/elbow_plot.png',dpi=300)
    plt.show()
    plt.close()

def sort_clusters(cluster_sort_df,data_cols):
    clusters = pd.DataFrame()
    for cluster in cluster_sort_df['GROUP'].astype(int).unique():
        clusters.loc[cluster,'GROUP'] = int(cluster)
        clusters.loc[cluster,'average'] = cluster_sort_df.loc[cluster_sort_df["GROUP"] == cluster,data_cols].to_numpy().mean()

    clusters.sort_values('average',axis=0,inplace=True,ascending=False,ignore_index=True)
    clusters['new_order'] = 'clust' + (clusters.index + 1).astype(str)
    return dict(zip(clusters['GROUP'], clusters['new_order'])),list(clusters['new_order'])

# the following piece of code is designed to determine if we are making pairwise comparisons (i.e. if for each of the treatments, there is matching timepoints or whatever)
def comparison_type(data_cols):
    treatments = []
    for i in data_cols: treatments.append(i.split('_')[0])
    treatments = list(set(treatments))

    treatment_numbers = []
    for treat in treatments: treatment_numbers.append(len([k for k in data_cols if treat in k]))
    treatment_numbers = list(set(treatment_numbers))
    if len(treatment_numbers) > 1: comparison = 'multiple'
    elif len(treatments) == len(data_cols): comparison = 'multiple'
    else: comparison = 'pairwise'

    return comparison,treatments

# sort col must be just one column
# in the future, might be nice to sort by multiple user-provided columns
# must provide either the sort_col and the order you want to sort that column,
    # or a sort_dict, if you want to sort multiple columns by strings
    # dict format is {'sort_columnA':'orderA','sort_columnB':'orderB'}
def sort_by_strings(df,sort_col=None,order=None,sort_dict=None,pre_sort_cols=[],post_sort_cols=[]):
    out_df = df.copy()
    
    if sort_col is not None: sort_col_list,sort_order = [sort_col],order
    elif sort_dict is not None: sort_col_list = list(sort_dict.keys())

    for col in sort_col_list:
        if sort_dict is not None: sort_order = sort_dict[col]
        rep = 0
        for i in sort_order:
            out_df[col] = out_df[col].replace(i,rep)
            rep += 1

    out_df.sort_values(pre_sort_cols + sort_col_list + post_sort_cols,ascending=True,inplace=True)

    for col in sort_col_list:
        if sort_dict is not None: sort_order = sort_dict[col]
        rep = 0
        for i in sort_order:
            out_df[col] = out_df[col].replace(rep,i)
            rep += 1
    return out_df