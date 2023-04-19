# Functions used to parse plate reader data
def get_auc(g):
    import numpy as np
    return np.trapz(g.OD)

# Functions to perform statistical tests
def parse_pval(pvalue):
    if pvalue < 0.0001:
        return '****'
    elif pvalue < 0.001:
        return '***'
    elif pvalue < 0.01:
        return '**'
    elif pvalue < 0.05:
        return '*'
    else:
        return 'ns'
    
def get_Tukey(df, factor1, factor2, var):
    import pandas as pd
    import statsmodels.api as sm
    from statsmodels.formula.api import ols
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
    
    anova_df = df[[factor1, factor2, var]]
    
    # Performing two-way ANOVA
    model = ols(f'{var} ~ C({factor1}) + C({factor2}) +\
    C({factor1}):C({factor2})',
            data=anova_df).fit()
    anova_res = sm.stats.anova_lm(model, type=2) # Should be checked to make sure the factors explain a significant amount of variation
  
    anova_df['combination'] = anova_df.apply(lambda row: str(row[factor1])+' / '+str(row[factor2]), axis=1)
    
    # perform multiple pairwise comparison (Tukey HSD)
    m_comp = pairwise_tukeyhsd(endog=anova_df[var], groups=anova_df['combination'], alpha=0.05)
    
    # convert tukeyhsd table to a DataFrame
    tukey_data = pd.DataFrame(data=m_comp._results_table.data[1:], columns = m_comp._results_table.data[0])

    # convert p-adj to ***
    tukey_data['stat_significance'] = tukey_data['p-adj'].apply(parse_pval)
    
    # filter comparison of respiration factor for each condition
    tukey_data['condition1'] = tukey_data.group1.apply(lambda x: x.split('/ ')[1])
    tukey_data['condition2'] = tukey_data.group2.apply(lambda x: x.split('/ ')[1])
    tukey_res = tukey_data[tukey_data.condition1 == tukey_data.condition2][['condition1', 'group1', 'group2', 'p-adj', 'reject', 'stat_significance']]
    
    return anova_res, tukey_res 