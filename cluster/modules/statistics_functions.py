
import itertools
import pandas as pd
from scipy.stats import kruskal,ranksums
from statsmodels.stats.multitest import multipletests


def kruskal_wilcoxon(data_names,data_list):
    statistics_df = pd.DataFrame()
    row = 0
    if len(data_names) > 2:
        statistics_df.loc[row,'test_comparison'] = ','.join(data_names)
        statistics_df.loc[row,'test_type'] = 'Kruskal-Wallis'
        statistic,pvalue = kruskal(*data_list,nan_policy='omit')
        statistics_df.loc[row,'statistic'],statistics_df.loc[row,'pvalue'] = statistic,pvalue

        if pvalue < 0.05:
            row_list,pval_list = [],[]
            for combo in itertools.combinations(range(0,len(data_list)),r=2):
                row += 1
                statistics_df.loc[row,'test_comparison'] = data_names[combo[0]] + ',' + data_names[combo[1]]
                statistics_df.loc[row,'test_type'] = 'Wilcoxon rank-sum'
                wil_stat,wil_pval = ranksums(data_list[combo[0]],data_list[combo[1]],alternative='two-sided')
                statistics_df.loc[row,'statistic'],statistics_df.loc[row,'pvalue'] = wil_stat,wil_pval
                row_list.append(row)
                pval_list.append(wil_pval)
            _,corrected_pval,_,_ = multipletests(pval_list,method='bonferroni')
            statistics_df.loc[row_list,'multiple_testing_method'] = 'bonferroni'
            statistics_df.loc[row_list,'padj'] = corrected_pval
    elif len(data_names) == 2: # not super sure if this is running as intended
        statistics_df.loc[row,'test_comparison'] = ','.join(data_names)
        statistics_df.loc[row,'test_type'] = 'Wilcoxon rank-sum'
        wil_stat,wil_pval = ranksums(data_list[0],data_list[1],alternative='two-sided',nan_policy='omit')
        statistics_df.loc[row,'statistic'],statistics_df.loc[row,'pvalue'] = wil_stat,wil_pval

    return statistics_df
