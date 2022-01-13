from __future__ import print_function
import sys 
import os
import argparse 
import pandas as pd
import numpy as np
from collections import defaultdict
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import plotly.express as px
import plotly.graph_objects as go
import plotly.offline as py
import gzip
import io
import seaborn as sns
from scipy import stats
import matplotlib.animation as animation
import matplotlib
from cyvcf2 import VCF
from interlap import InterLap
from sklearn import metrics
from rpy2.robjects.packages import STAP
from sklearn.ensemble import RandomForestClassifier
from math import e # Euler's constant 

## Plotting Colors
AR_color = sns.color_palette("colorblind")[0]
AD_color = sns.color_palette("colorblind")[2]
HI_color = sns.color_palette("colorblind")[3]
Other_color = sns.color_palette("colorblind")[1]
OR_color = sns.color_palette("colorblind")[4]
NonEss_color = sns.color_palette("colorblind")[9]


deciles = ["0.0-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5","0.5-0.6","0.6-0.7","0.7-0.8","0.8-0.9","0.9-1.0"]

def get_decile_score(score):
    
    if score == "Missing":
        return("None")
    
    for d in deciles:
        
        mind = float(d.strip().split("-")[0])
        maxd = float(d.strip().split("-")[1])
        
        if maxd == 1.0:
            maxd += 0.1
            
        maxd += 0.1 if maxd == 1.0 else 0.0
        
        if score >= mind and score < maxd:
            return(d)

def get_score_bin(row, score_column):
    """
    add_score_bin
    =============
    Method to add the a score/decile bin for a particular row in a dataframe. (Usually used for an apply function)
     it will iterate over deciles (10 bins) from 0.0 to 1.0. It will assign the bin/decile where the
     score is >= the lower limit and the score is < the upper limit: [lowerlimit - upper limit)
     The method will return the assigned decile/bin [lowerlimit - upper limit)

    Parameters:
    -----------
    1) row:           (pandas Series) A row in a pandas dataframe
    2) score_column: (str) The name of the column with the score to get bin for

    Returns:
    ++++++++
    1) (str) A decile/bin where the score is between the lower and upper limit.
            example: 0.3-0.4, where 0.3 is the lower limit and 0.4 is the upper limit
    NOTE: If the row does not have a score a "-1.0" bin is returned
    """

    decile_bins = ["0.0-0.1",
                   "0.1-0.2",
                   "0.2-0.3",
                   "0.3-0.4",
                   "0.4-0.5",
                   "0.5-0.6",
                   "0.6-0.7",
                   "0.7-0.8",
                   "0.8-0.9",
                   "0.9-1.0"]


    for b in decile_bins:
        min_score = float(b.strip().split("-")[0])
        max_score = float(b.strip().split("-")[1])

        max_score += 0.1 if max_score == 1.0 else 0.0

        if float(row[score_column]) >= min_score and float(row[score_column]) < max_score:
            return(b)

    return("-1.0")


def get_alternative_gene_symbols(gene_info_file):
    """
    get_alternative_gene_symbols
    ============================
    This method is used to parse a file with accepted and alternative gene symbols and create a mapping between the pair.
     This function is set up to work with the HGNC Database mapping file. http://www.genenames.org/. The HGNC has put together
     A large data file that contains mapping between genes across the different data providers and each of the unique ids that
     each provider makes. Here, the function will use the accepted symbols and the alternative symbols from this file
     to create a mapping between the symbols. File url: ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt

    Parameters:
    ----------
    1) gene_info_file: (str) The file path to the HGNC symbol mapping file. NOTE: The file needs a "#" at the beginning of the 1st line in order to work.

    Returns:
    +++++++
    1) (dictionary) A dict where keys represent a single gene symbol, and values represent a set of synonymous symbols. The keys will include the accepted and
                    alternative gene symbols.
    2) (dictionary) A dict where keys are the main gene symbol and values are sets of alternative symbols
    3) (dictionary) A dict where keys are alt symbols and values are the main symbol it maps to

    NOTE: Any gene with no alternative symbols will be included as a key with the value as an empty set
    """

    ## Open the alt gene file
    gene_info_fh = gzip.open(gene_info_file, "rt", encoding = "utf-8") if gene_info_file.endswith(".gz") else io.open(gene_info_file, "rt", encoding = "utf-8")

    ## Main dictionary for alternative symbols
    alternative_gene_symbols = defaultdict(set)
    main_to_alt_gene_symbols = defaultdict(set)
    alt_to_main_gene_symbols = defaultdict(str)

    ## Get an index from the header
    header_line = gene_info_fh.readline()
    assert "#" in header_line, "The alternative gene symbol file does not have a header staring with '#'. A header is required. Please add a header and try again"
    header_index = header_line.strip().replace("#","").split("\t")

    ## Iterate over the file
    for line in gene_info_fh:

        ## Get a dictionary for the current line
        line_dict = dict(zip(header_index, line.strip().split("\t")))

        ## Get the official gene symbol and the synonymous symbols
        gene_symbol = line_dict["symbol"]

        ## Get that alternative gene symbols
        synonymous_symbols = set(x for x in line_dict["alias_symbol"].strip().replace("\"","").split("|")) if line_dict["alias_symbol"] else set()
        ## Add any previous symbols
        synonymous_symbols.update(set(x for x in line_dict["prev_symbol"].strip().replace("\"","").split("|")) if line_dict["prev_symbol"] else set())
        ## Add the current symbol to the set
        synonymous_symbols.add(gene_symbol)

        ## Add to dict
        alternative_gene_symbols[gene_symbol].update(synonymous_symbols)
        main_to_alt_gene_symbols[gene_symbol].update(synonymous_symbols)
        alt_to_main_gene_symbols[gene_symbol] = gene_symbol

        ## iterate over each synonymous symbol
        for synonymous_symbol in synonymous_symbols:

            ## add synonymous symbol
            alternative_gene_symbols[synonymous_symbol].update(set([x for x in synonymous_symbols if x != synonymous_symbol] + [gene_symbol]))

            ## map alt symbol to main
            alt_to_main_gene_symbols[synonymous_symbol] = gene_symbol


    ## Close fh
    gene_info_fh.close()

    return(alternative_gene_symbols, main_to_alt_gene_symbols, alt_to_main_gene_symbols)


def get_normalized_matrix_per_label(groupby_object, mutation_rate_dict):
    """
    get_normalized_matrix_per_label
    ===============================
    This method is used to get the marginal distribution for each reference allele at each bin (exp: 0.8-1.0), 
     and normalize the distribution using Ref to alt count divided by total from the marginal distribution.
    Example Matrix:
      
      D+4
              ALT
              ---
        | A | C | G | T |
       -|---|---|---|---|
       A| # | # | # | # |
    R| -|---|---|---|---|
    E| C| # | # | # | # |
    F| -|---|---|---|---|
       G| # | # | # | # |
       -|---|---|---|---|
       T| # | # | # | # | 
       ------------------
    Marginal Distribution is based on Ref allele (Row)
    Normalized marginal distribution counts will add up to 1 
    Parameters:
    -----------
    1) groupby_object:     (Pandas df) A subseted pandas dataframe that represents ref to alt counts for a specific SpliceAI delta score bin (exp: "0.0-0.2")
    2) mutation_rate_dict: (dict)      A dictionary that contains the ref allele specific mutation rates. (Updated by this function)
    Returns:
    ++++++++
    Nothing is returned. Instead, the mutation_rate_dict is dynamically updated. 
    """
    
    label = groupby_object.delta_score_bin.unique()[0]

    ## Create an empty matrix
    zeroton_matrix = np.zeros((len(groupby_object.ref.unique()),len(groupby_object.alt.unique())), dtype= np.float)
    zeroton_plus_singleton_matrix = np.zeros((len(groupby_object.ref.unique()),len(groupby_object.alt.unique())), dtype= np.float)
    row_index = dict(zip(groupby_object.ref.unique(),[0,1,2,3]))
    column_index = dict(zip(groupby_object.alt.unique(),[0,1,2,3]))

    ## create a matrix for ref to alt counts
    for row in groupby_object.itertuples():
        ## If ref is the same as alt (no mutation)
        if row.ref == row.alt:
            
            ## Zeroton Counts
            zeroton_matrix[row_index[row.ref], column_index[row.alt]] += row.zerotons

            ## Zeroton plus Singleton counts
            zeroton_plus_singleton_matrix[row_index[row.ref], column_index[row.alt]] += row.zerotons

        else:
            ## Zeroton Counts
            ### If ref to alt, counts = non zertons
            zeroton_matrix[row_index[row.ref], column_index[row.alt]] = row.non_zerotons

            ## Zeroton plus Singleton counts
            ### Add counts to alt for all non zerotons and non singletons
            zeroton_plus_singleton_matrix[row_index[row.ref], column_index[row.alt]] = row.non_zerotons_plus_singletons
            ### Add singleton counts to the ref to ref cell
            zeroton_plus_singleton_matrix[row_index[row.ref], column_index[row.ref]] += row.singletons

    ## Get the normalized row counts
    zeroton_row_normalized_matrix = zeroton_matrix / zeroton_matrix.sum(axis = 1)[:, np.newaxis]
    zeroton_plus_singleton_row_normalized_matrix = zeroton_plus_singleton_matrix / zeroton_plus_singleton_matrix.sum(axis = 1)[:, np.newaxis]

    ## get mutation rates based on 1 - ref to ref normalized count
    zeroton_mutation_rates = [1.0 - zeroton_row_normalized_matrix[row_index[x[0]],column_index[x[1]]] for x in zip(row_index,column_index)]
    zeroton_plus_singletons_mutation_rates = [1.0 - zeroton_plus_singleton_row_normalized_matrix[row_index[x[0]],column_index[x[1]]] for x in zip(row_index,column_index)]

    ## Mutation rate to dict
    mutation_rate_dict[label]["zeroton"] = dict(zip(groupby_object.ref.unique(),zeroton_mutation_rates))
    mutation_rate_dict[label]["zeroton_plus_singleton"] = dict(zip(groupby_object.ref.unique(),zeroton_plus_singletons_mutation_rates))


def get_per_bin_rate_ratio(df,
                           bin_column,
                           case_group,
                           control_group,
                           cohort_column,
                           score_column,
                           low_score_cutoff,
                           high_score_cutoff,
                           exclude_bins = [],
                           offset_zero = False):
    """


    1) df: (pandas DataFrame)
    2) bin_column: (str) The column that represents different score bins to get the rate ratio for. (Example: constraint_percentile_bin)
    3)
    ) offset_zero: (bool) whether or not to offset scores by 0.5 to avoid divide by zero errors
    """

    r_get_rateratio = STAP(r_rateratio_func, "r_pkg")

    rr_list = []

    for xbin, group in df.groupby(bin_column):

        ## Skip an excluded specific bins
        if xbin in exclude_bins:
            continue


        print("Bin {} Rate Ratio Info".format(xbin))


        ## Get the group totals for the current bin
        group_totals = group.groupby(cohort_column)[score_column].count().reset_index()

        if offset_zero:
            ## Offset all scores with 0.1 to avoid dividing by 0
            group_totals[score_column] += 2

        ## Get the per group count of variants that are <= to the low score cutoff
        group_low_counts = group.loc[(group[score_column] < low_score_cutoff) & (group[score_column] >= 0)].groupby(cohort_column)[score_column].count().reset_index()

        ## Check if case and control group exists
        if case_group not in group_low_counts[cohort_column].tolist():
            group_low_counts.loc[len(group_low_counts.index)] = [case_group, 0]
        if control_group not in group_low_counts[cohort_column].tolist():
            group_low_counts.loc[len(group_low_counts.index)] = [control_group, 0]

        if offset_zero:
            ## Offset all scores with 0.1 to avoid dividing by 0
            group_low_counts[score_column] += 1

        ## Get the per group count of variants that are >= to the high score cutoff
        group_high_counts = group.loc[group[score_column] >= high_score_cutoff].groupby(cohort_column)[score_column].count().reset_index()

        ## Check if case and control group exists
        if case_group not in group_high_counts[cohort_column].tolist():
            group_high_counts.loc[len(group_high_counts.index)] = [case_group, 0]
        if control_group not in group_high_counts[cohort_column].tolist():
            group_high_counts.loc[len(group_high_counts.index)] = [control_group, 0]

        if offset_zero:
            ## Offset all scores with 0.1 to avoid dividing by 0
            group_high_counts[score_column] += 1


        ## Get rate ratios
        (low_rate_ratio,
         low_pvalue,
         low_ci_low,
         low_ci_high) = r_get_rateratio.get_rateratio(float(group_low_counts.loc[group_low_counts[cohort_column] == case_group][score_column].values[0]),
                                                      float(group_low_counts.loc[group_low_counts[cohort_column] == control_group][score_column].values[0]),
                                                      float(group_totals.loc[group_totals[cohort_column] == case_group][score_column].values[0]),
                                                      float(group_totals.loc[group_totals[cohort_column] == control_group][score_column].values[0]))

        (high_rate_ratio,
         high_pvalue,
        high_ci_low,
        high_ci_high) = r_get_rateratio.get_rateratio(float(group_high_counts.loc[group_high_counts[cohort_column] == case_group][score_column].values[0]),
                                                      float(group_high_counts.loc[group_high_counts[cohort_column] == control_group][score_column].values[0]),
                                                      float(group_totals.loc[group_totals[cohort_column] == case_group][score_column].values[0]),
                                                      float(group_totals.loc[group_totals[cohort_column] == control_group][score_column].values[0]),)

        ## Add rate ratio info to list
        rr_list.append([xbin,
                        low_rate_ratio,
                        low_pvalue,
                        low_ci_low,
                        low_ci_high,
                        high_rate_ratio,
                        high_pvalue,
                        high_ci_low,
                        high_ci_high])


        print(group_totals)
        print(group_low_counts)
        print(group_high_counts)
        print("\n\n")


    ## Create rate ratio dataframe
    rr_df = pd.DataFrame(rr_list,
                         columns = ["score_bin",
                                   "unconstrained_rr",
                                   "unconstrained_p_value",
                                   "unconstrained_ci_bottom",
                                   "unconstrained_ci_top",
                                   "constrained_rr",
                                   "constrained_p_value",
                                   "constrained_ci_bottom",
                                   "constrained_ci_top"])
    return(rr_df)


def plot_rate_ratio(rr_df, 
                    bin_column, 
                    constrained_rr_column, 
                    unconstrained_rr_column,
                    constrained_label,
                    unconstrained_label,
                    constrained_color = "firebrick",
                    unconstrained_color = "dodgerblue",
                    legend_title = "Variant Classification",
                    xlabel = "Decile/Bin",
                    ylabel = "Rate Ratio",
                    title = "Rate Ratio plot",
                    plot_ymax = -1,
                    constrained_error = [],
                    unconstrained_error = []):
    
    rr_df["cons_x_axis"] = np.arange(0,len(rr_df[bin_column])) + 0.25
    rr_df["uncons_x_axis"] = np.arange(0,len(rr_df[bin_column])) - 0.25
    
    fig, ax = plt.subplots(figsize=(20,15))


    sns.pointplot(x = cons_x_axis, 
                  y = constrained_rr_column, 
                  data = rr_df,
                  join=False, 
                  color = constrained_color, 
                  scale = 2)
    
    if constrained_error:
        plt.errorbar(x = rr_df["cons_x_axis"], 
             y = rr_df[constrained_rr_column], 
             yerr = constrained_error, 
             fmt='none', 
             color = constrained_color,
             capsize = 4)


    sns.pointplot(x = uncons_x_axis, 
                  y = unconstrained_rr_column, 
                  data = rr_df,
                  join=False, 
                  color = unconstrained_color,
                  scale = 2)
    
    if unconstrained_error:
        plt.errorbar(x = rr_df["uncons_x_axis"], 
             y = rr_df[unconstrained_rr_column], 
             yerr = unconstrained_error, 
             fmt='none', 
             color = unconstrained_color,
             capsize = 4)

    ax.axhline(1, ls='--', color = "black")

    import matplotlib.patches as mpatches 
    from matplotlib.lines import Line2D


    high_legend = Line2D([0], [0], color = "white", marker='o', markersize=15, markerfacecolor=constrained_color, label = constrained_label)
    low_legend = Line2D([0], [0], color = "white",  marker='o', markersize=15, markerfacecolor=unconstrained_color, label = unconstrained_label)
    plt.legend( handles = [high_legend, low_legend],bbox_to_anchor=(1.01,1), loc="upper left", fontsize = 18, title = legend_title, title_fontsize = 18)
    
    
    plt.ylabel(ylabel, fontsize = 20)
    plt.xlabel(xlabel, fontsize = 20)
    plt.title(title, fontsize = 22)
    
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    
    if plot_ymax > 0:
        plt.ylim(0,plot_ymax)



r_odds_ratio_func = """

get_odds_ratio = function(groupA_event_count, groupA_no_event_count, groupB_event_count, groupB_no_event_count){

#    get_odds_ratio
#    ==============
#    Method to get the odds ratio using a fisher's exact test. 
#    
#    Parameters:
#    -----------
#    1) groupA_event_count:    (int) The number values for group A for the event
#    2) groupA_no_event_count: (int) The number values for group A not for the event
#    3) groupB_event_count:    (int) The number values for group B for the event
#    4) groupB_no_event_count: (int) The number values for group B not for the event
#    
#    Returns:
#    ++++++++
#    1) (float) odds ratio
#    2) (float) p_value 
#    3) (float) lower bound of the 95% confidence interval
#    4) (float) upper bound of the 95% confidence interval


    m = matrix(c(groupA_event_count,
                 groupB_event_count,
                 groupA_no_event_count,
                 groupB_no_event_count), 
               2,2, 
               dimnames=list(c("Group A","Group B"),c("Event","No Event")))
    
    print("Contingency Table")
    print(m)
                                                                                                                        
    or_object = fisher.test(m)
    #print(or_boject)
    
    or_value = or_object$estimate
    p_value = or_object$p.value
    ci_lower = or_object$conf.int[1] 
    ci_upper = or_object$conf.int[2] 

    return(c(or_value, p_value, ci_lower, ci_upper))   
}

"""


r_rateratio_func = """


get_rateratio = function(value_one, value_two, total_one, total_two){

#    get_rateratio
#    =============
#    Method to get the rate ratio using a Poisson exact test.
#
#    Parameters:
#    -----------
#    1) value_one: (int) The number of interesting values in group one
#    2) value_two: (int) The number of interesting values in group two
#    3) total_one: (int) The total number of values in group one
#    4) total_two: (int) The total number of values in group two
#
#    Returns:
#    ++++++++
#    1) (float) rate ratio for the fold change in rate between group one and group two
#    2) (float) bottom of the 95% confidence interval
#    3) (float) top of the 95% confidence interval


    library(rateratio.test)
    rr = rateratio.test(c(value_one,value_two),c(total_one,total_two))
    print(rr)

    conf_int = c(rr$conf.int[1],rr$conf.int[2])
    rate_ratio = rr$estimate["Rate Ratio"]
    p_value = rr$p.value

    return(c(rate_ratio, p_value, conf_int))

}
"""


def get_roc_vars(query_df, label_column, score_column,):
    """
    get_roc_var
    ==========
    Method to get the variables for a roc curve. Variables include:
     - False positive rate
     - True positive rate
     - Thresholds
     - Area under the curve 
    This method uses sklearn's metrics.roc_curve module to get the roc variables.
    
    Parameters:
    -----------
    1) query_df: The dataframe to get the roc variables for.
    2) label_column: The name of the column in the query dataframe that represents the true positive and true 
                     negative labels.
    3) score_column: The name of the column in the query dataframe that represents the scores associated with 
                     each label. (An example would be the Observed_over_Expected column)
    
    Returns:
    ++++++++
    1) A numpy array of false positive rates
    2) A numpy array of true positive rates
    3) A numpy array of thresholds used to generate the false positive and true positive rates
    4) The area under the curve (AUC). 
    """
    
    
    ## Get a copy of the dataframe, with only the labels and the scores
    query_roc_df = query_df[[label_column,score_column]].copy()

    ## Remove all NA or NaN values from the dataframe
    query_roc_df = query_roc_df.dropna()
    ## Cast all values of the subseted dataframe to float
    query_roc_df = query_roc_df.astype("float")
    ## Get the false positive rate and true positive rate using sklearn roc_curve module 
    fp_rate, tp_rate, roc_thresholds = metrics.roc_curve(query_roc_df[label_column],query_roc_df[score_column])
    ## Get the area under the curve 
    area_under_curve = metrics.roc_auc_score(query_roc_df[label_column],query_roc_df[score_column])
    
    ##Yoden's J Statistic
    jindex = np.argmax(tp_rate - fp_rate)
    jscore = tp_rate[jindex] - fp_rate[jindex]
    j_thresh_index = roc_thresholds[jindex]

    print('ROC AUC=%.3f' % (area_under_curve))

    return(fp_rate, tp_rate, roc_thresholds, area_under_curve, jscore, j_thresh_index, jindex)




def plot_roc_curve(list_of_fpr, 
                   list_of_tpr, 
                   area_under_curve, 
                   plot_title="title", 
                   line_color = "dodgerblue", 
                   return_plot=False,
                   save_plot = False, 
                   output_name = "output.png"):
    """
    plot_roc_curve
    ==============
    Method to plot a roc curve from a list of false positive rates and true positive rates.
    
    Parameters:
    ----------
    1) list_of_fpr: A list of false positive rates
                    NOTE: This list should corresponds to a list of true positive rates create from the same thresholds
    2) list_of_tpr: A list of true positive rates
                    NOTE: This list should corresponds to a list of false positive rates create from the same thresholds
    3) area_under_curve: The Area under the curve for the ROC. 
    4) plot_title: The title of the plot. (Default = "title")
    5) line_color: Color of the ROC Curve line. This needs to be a color from the matplotlib color palette. 
                   (Default = "dodgerblue")
    6) return_plot: Whether to return the plot or not. If True, the plot will be returned. If False, the plot will be plotted. 
                (Default = False)
    7) save_plot: Whether or not to save the plot rather than plot it. (Default = False)
    8) output_name: If save_plot is set to True, the name of the output plot. (Default = output.png)
    """
    
    fig, ax = plt.subplots(1, figsize = (15,10))
    
    sns.set(font='Arial', style = "white")
    sns.set_context("paper", font_scale = 2)

    ax.plot([0, 1], [0, 1], 'k--', )
    ax.plot(list_of_fpr,list_of_tpr, color=line_color, linewidth=2, label = "AUC = %.3f" % (area_under_curve))

    ax.set_title(plot_title, fontsize = 30)
    ax.set_xlabel('False positive rate', fontsize = 28)
    ax.set_ylabel('True positive rate', fontsize = 28)
    ax.legend(loc="lower right",prop={'size': 22})
    
    ax.tick_params(axis="x", labelsize=22)
    ax.tick_params(axis="y", labelsize=22)
    
    
    
def plot_combined_roc_curve(list_of_kargs,
                            plot_title="title",
                            save_plot = False, 
                            output_name = "output.png"):
    """
    plot_combined_roc_curve
    =======================
    Method to plot multiple roc curves on the same plot.
    
    Parameters:
    ----------
    1) list_of_kargs: A list of keyword arguments, where the kargs is a dictionary, with each karg dict 
                       representing information for a single ROC curve to plot. The number of ROC curves 
                       to plot on the same plot will be directly linked to the number of karg dicts in this
                       list.
        Required keys in each karg dict:
          - fpr_list: The list of false positive rates
          - tpr_list: The list of true positive rates
          - auc: The area under the curve 
          - line_label: The label to associate with the ROC curve line
          - line_color: The color of the ROC curve line
        If any of these keys are missing for any of the karg dicts, the plotting will fail.
        
    2) plot_title: The title of the plot. (Default = "title")
    3) save_plot: Whether or not to save the plot rather than plot it. (Default = False)
    4) output_name: If save_plot is set to True, the name of the output plot. (Default = output.png)
    """
    
    ## Required keyword argument keys
    required_karg_set = set(["fpr_list",
                             "tpr_list",
                             "auc",
                             "jscore",
                             "j_threshold",
                             "jindex",
                             "line_label",
                             "line_color"])
    
    fig, ax = plt.subplots(1, figsize = (15,10))
    
    sns.set(font='Arial', style = "white")
    sns.set_context("paper", font_scale = 2)
    sns.despine()
                           
    for karg_dict in list_of_kargs:
        ## Check for required kargs
        assert all(True if x in required_karg_set else False for x in karg_dict.keys()), "!!ERROR!! Mislabeled keys in kargs dictionary. Required kargs = {}. Correct for the mislabeled keyword arguments and try again.".format(", ".join(list(required_karg_set)))
        assert all(True if x in karg_dict.keys() else False for x in  required_karg_set), "!!ERROR!! Missing keys in kargs dictionary. Required kargs = {}. Correct for the missing keyword arguments and try again.".format(", ".join(list(required_karg_set)))
        
        
        ## Add plot
        ax.plot(karg_dict["fpr_list"],
                karg_dict["tpr_list"], 
                color = karg_dict["line_color"], 
                linewidth = 1.5,
                label = "{} = {:.3f}, max J = {:.3f}".format(karg_dict["line_label"],karg_dict["auc"],karg_dict["jscore"]))
        
        ax.scatter(karg_dict["fpr_list"][karg_dict["jindex"]],
                karg_dict["tpr_list"][karg_dict["jindex"]], 
                s = 100, 
                color = karg_dict["line_color"] )
    

    ## Add diagonal line
    ax.plot([0, 1], [0, 1], 'k--')
    
    ## Set plot labels
    ax.set_title(plot_title, fontsize = 25)
    ax.set_xlabel('False positive rate', fontsize = 30)
    ax.set_ylabel('True positive rate', fontsize = 30)
    ax.legend(title = "ROC AUC, Yoden's J",loc="lower right",prop={'size': 20})
    
    ax.tick_params(axis="x", labelsize=25)
    ax.tick_params(axis="y", labelsize=25)
    
    if save_plot:
        plt.savefig(output_name)
        plt.close()
    else:
        plt.show()
        
        
        
        
def get_pr_vars(query_df, label_column, score_column):
    """
    get_pr_var
    ==========
    Method to get the variables for a precision recall curve. Variables include:
     - precision
     - recall
     - Thresholds
    This method uses sklearn's metrics.precision_recall_curve module to get the pr variables.
    
    Parameters:
    -----------
    1) query_df: The dataframe to get the pr variables for.
    2) label_column: The name of the column in the query dataframe that represents the true positive and true 
                     negative labels.
    3) score_column: The name of the column in the query dataframe that represents the scores associated with 
                     each label. (An example would be the Observed_over_Expected column)
    
    Returns:
    ++++++++
    1) A numpy array of precisions
    2) A numpy array of recall
    3) A numpy array of thresholds used to generate the precision and recall lists
    4) The PR AUC
    5) The Average Precision Score
    """
    
    
    ## Get a copy of the dataframe, with only the labels and the scores
    query_pr_df = query_df[[label_column,score_column]].copy()
    ## Remove all NA or NaN values from the dataframe
    query_pr_df = query_pr_df.dropna()
    ## Cast all values of the subseted dataframe to float
    query_pr_df = query_pr_df.astype("float")
    ## Get the false positive rate and true positive rate using sklearn roc_curve module 
    precision_list, recall_list, pr_thresholds = metrics.precision_recall_curve(query_pr_df[label_column],query_pr_df[score_column])
    ## Get AUC
    pr_auc = metrics.auc(recall_list, precision_list)
    ## Get average precision score 
    average_per_score = metrics.average_precision_score(query_pr_df[label_column],query_pr_df[score_column])
    
    
    print('PR AUC=%.3f, Avg. Precision Score=%.3f' % (pr_auc, average_per_score))
    
    return(precision_list, recall_list, pr_thresholds, pr_auc, average_per_score)


def plot_combined_pr_curve(list_of_kargs,
                            base_line = {"x": [0,1], "y": [0.5,0.5]},
                            plot_title="title",
                            legend_loc = "lower center",
                            save_plot = False, 
                            output_name = "output.png"):
    """
    plot_combined_pr_curve
    =======================
    Method to plot multiple roc curves on the same plot.
    
    Parameters:
    ----------
    1) list_of_kargs: A list of keyword arguments, where the kargs is a dictionary, with each karg dict 
                       representing information for a single ROC curve to plot. The number of ROC curves 
                       to plot on the same plot will be directly linked to the number of karg dicts in this
                       list.
        Required keys in each karg dict:
          - fpr_list: The list of false positive rates
          - tpr_list: The list of true positive rates
          - auc: The area under the curve 
          - line_label: The label to associate with the ROC curve line
          - line_color: The color of the ROC curve line
        If any of these keys are missing for any of the karg dicts, the plotting will fail.
        
    2) plot_title: The title of the plot. (Default = "title")
    3) save_plot: Whether or not to save the plot rather than plot it. (Default = False)
    4) output_name: If save_plot is set to True, the name of the output plot. (Default = output.png)
    """
    
    ## Required keyword argument keys
    required_karg_set = set(["per_list","rec_list",
                             "pr_auc","av_pr_score",
                             "line_label","line_color"])

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    plt.rcParams['svg.fonttype'] = 'none'

    fig, ax = plt.subplots(1, figsize = (15,10))

    sns.set(font='Arial', style = "white")
    sns.set_context("paper", font_scale = 2)
    sns.despine()
                           
    for karg_dict in list_of_kargs:
        ## Check for required kargs
        assert all(True if x in required_karg_set else False for x in karg_dict.keys()), "!!ERROR!! Mislabeled keys in kargs dictionary. Required kargs = {}. Correct for the mislabeled keyword arguments and try again.".format(", ".join(list(required_karg_set)))
        assert all(True if x in karg_dict.keys() else False for x in  required_karg_set), "!!ERROR!! Missing keys in kargs dictionary. Required kargs = {}. Correct for the missing keyword arguments and try again.".format(", ".join(list(required_karg_set)))
        
        
        ## Add plot
        ax.plot(karg_dict["rec_list"],
                karg_dict["per_list"], 
                color = karg_dict["line_color"], 
                linewidth = 3,
                label = "{} = {:.3f}, $\\bf{:.3f}$".format(karg_dict["line_label"],karg_dict["pr_auc"], karg_dict["av_pr_score"]))
    
        ax.plot()
    
    ## Add baseline
    ax.plot(base_line["x"], base_line["y"], 'k--', lw = 3, dashes=(4, 4))
    
    ## Set plot labels
    ax.set_title(plot_title, fontsize = 55)
    ax.set_xlabel('Recall', fontsize = 50)
    ax.set_ylabel('Precision', fontsize = 50)


    ax.legend(title = "PR AUC, $\\bf{Ave. PR Score}$",loc=legend_loc,prop={'size': 24}, title_fontsize = 26)
    
    ax.tick_params(axis="x", labelsize=40)
    ax.tick_params(axis="y", labelsize=40)
    
    if save_plot:
        plt.savefig((output_name.replace(".svg","") + ".svg"), dpi=300)
        plt.close()
    else:
        plt.show()
        

        
def ConSpliceML_train(training_df,
                      feature_col_names = [],
                      label_col = "label",
                      n_estimators = 100,
                      random_state = 123):

    ## RF models
    rf = RandomForestClassifier(random_state = random_state, n_estimators = n_estimators)

    rf.fit(training_df[feature_col_names],training_df[label_col].to_numpy())
    
    
    return(rf)



def ConSpliceML_score(rf,
                      var_df, 
                      feature_col_names, 
                      new_col_name = "RF_Score", 
                      sort_by_cols = []):
        

         
    var_df[new_col_name] = rf.predict_proba(var_df[feature_col_names])[:,1]
    scored_df = var_df

    return(scored_df)


def enrichment_by_bin(df, score_col, plabel_col):
    """
    enrichment_by_bin
    =================
    Method to calculate an enrichment score of pathogenic vs benign variants for 10 bins
     a score of 0.0 and 1.0.
    
    Enrichment based on odds ratio:
        
        (Patho vars in bin [i] / benign vars in in bin [i]) 
        /
        (Patho vars not in bin [i] / benign vars not in bin [i[]])
     
                      
    NOTE: Pathogenic labels are expected to be 1
          Benign labeles are expected to be 0
    
    Parameters:
    -----------
    1) df:        (pandas DF) A dataframe with a score column and a pathogenic/benign label column
    2) score_col:       (str) The name of the score column in the df
    3) plabel_col:      (str) The name of the pathogenic/benign column in df
    
    Returns:
    ++++++++
    1) (pandas DF) A pandas dataframe with an "Enrichment" column and a "bins" column. 
                    The "Enrichment" column holds the pathogenic enrichment score for 
                    the associated bin.
    """
    
    bins =["0.0-0.1",
           "0.1-0.2",
           "0.2-0.3",
           "0.3-0.4",
           "0.4-0.5",
           "0.5-0.6",
           "0.6-0.7",
           "0.7-0.8",
           "0.8-0.9",
           "0.9-1.1"]

    
    
    ## Get the total number of patho and benign variants
    total_benign = df.loc[df[plabel_col] == 0].shape[0]
    total_patho = df.loc[df[plabel_col] == 1].shape[0]
    ## Offset to deal with zero values
    offset = 0.5
    enrichment_list = []
    se_list = []
    ci_lower_list = []
    ci_upper_list = []
    for b in bins:
        
        ## Get bin min and max values
        min_b, max_b = map(float,b.strip().split("-"))
        
        ## Get variant counts for current bin
        bin_df = df.loc[(df[score_col] >= min_b) & (df[score_col] < max_b)]
        total_bin_count = bin_df.shape[0] + (offset * 2)
        bin_patho_count = bin_df.loc[bin_df[plabel_col] == 1].shape[0] + offset
        bin_benign_count = bin_df.loc[bin_df[plabel_col] == 0].shape[0] + offset
        
        ## Make sure patho and benign bin combine to equal the total bin count
        assert bin_patho_count + bin_benign_count == total_bin_count
        
        
        
        ## Get variant counts for not in current bin
        ## subset df with min and max bin values
        notbin_df = df.loc[(df[score_col] < min_b) | (df[score_col] >= max_b)] 
        total_notbin_count = notbin_df.shape[0] + (offset * 2)
        notbin_patho_count = notbin_df.loc[notbin_df[plabel_col] == 1].shape[0] + offset
        notbin_benign_count = notbin_df.loc[notbin_df[plabel_col] == 0].shape[0] + offset
        
        ## Make sure patho and benign bin combine to equal the total bin count
        assert bin_patho_count + bin_benign_count == total_bin_count
        
        ## Make sure the in and outside bin counts match
        assert (bin_patho_count - offset) + (notbin_patho_count - offset) == total_patho, "{} + {} != {}".format((bin_patho_count - ofsset), (notbin_patho_count - ofsset), total_patho)
        assert (bin_benign_count - offset) + (notbin_benign_count - offset) == total_benign, "{} + {} != {}".format((bin_benign_count - ofsset), (notbin_benign_count - ofsset), total_benign)
        
        ## Get the proportion of patho to benign variants in the bin
        bin_patho_over_benign = bin_patho_count / bin_benign_count
        
        ## Get the proportion of patho to benign variants not in the bin
        notbin_patho_over_benign = notbin_patho_count / notbin_benign_count
        
        ## Calculate OR
        odds_ratio = bin_patho_over_benign / notbin_patho_over_benign 
        
        
        ## Calculate 95% CI from standard error 
        se = np.sqrt(((1/bin_patho_count) + (1/bin_benign_count) + (1/notbin_patho_count) + (1/notbin_benign_count)))


        
        #from math import e ## Euler's constant 
        ### np.log == ln
        ci_lower_bound = pow(e, ((np.log(odds_ratio) ) - (1.96 * se ) ) )
        ci_upper_bound = pow(e, ((np.log(odds_ratio) ) + (1.96 * se ) ) )
    
        
        
        ## Add enrichment score to the list
        enrichment_list.append(odds_ratio)
        se_list.append(se)
        ci_lower_list.append(ci_lower_bound)
        ci_upper_list.append(ci_upper_bound)
        
    
    ## Update last bin label
    bins[-1] = "0.9-1.0"

    ## Create enrichment df
    enrichment_df = pd.DataFrame(np.array([enrichment_list, 
                                           bins, 
                                           se_list, 
                                           ci_lower_list, 
                                           ci_upper_list]).T,
                                 columns = ["Enrichment","bins","se", "ci_lower", "ci_upper"])
    
    enrichment_df.Enrichment = enrichment_df.Enrichment.astype(float)
    enrichment_df.bins = enrichment_df.bins.astype(str)
    enrichment_df.se = enrichment_df.se.astype(float)
    enrichment_df.ci_lower = enrichment_df.ci_lower.astype(float)
    enrichment_df.ci_upper = enrichment_df.ci_upper.astype(float)
    
    ## Return enrichment df
    return(enrichment_df)


def plot_conspliceml_enrichment(enrich_df, 
                                score_col, 
                                bin_col, 
                                ci_lower_bound_col, 
                                ci_upper_bound_col, 
                                title = "",
                                plot_value_offset = 2):
    
    """
    plot_conspliceml_enrichment
    ===========================
    This method takes an enrichment df created by the 'enrichment_by_bin' method and plots the 
     pathogenic fold enrichment by bin.
     
    Parameters:
    -----------
    1) enrich_df:     (pandas DF) An enrichment df from the "enrichment_by_bin" method
    2) score_col:           (str) The name of the score column (Typically = "Enrichment")
    3) bin_col:             (str) The name of the bin column (Typically = "bins")
    4) ci_lower_bound_col:  (str) The name of the column with the 95% CI lower bound
    5) ci_upper_bound_col:  (str) The name of the column with the 95% CI upper bound
    6) title:               (str) The title of the plot
    7) plot_value_offset: (Float) The offset value to plot the OR enrichment score on the bar from the top
    
    
    """
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    plt.rcParams['svg.fonttype'] = 'none'

    fig, ax = plt.subplots( figsize = (20,10))
    sns.set(font='Arial', style = "white")
    sns.set_context("paper", font_scale = 2)
    sns.despine()
    fig.tight_layout(pad=2.0)

    ## Plot the shifted Enrichment scores with a score 
    bars = sns.barplot(x = bin_col, 
                       y = score_col, 
                       data = enrich_df,
                       color = sns.color_palette()[9],
                       order = enrich_df[bin_col].tolist())
    
    for i, row in enumerate(enrich_df.itertuples()):
        
        #bars.text(i, getattr(row, score_col) - ((getattr(row, score_col)) / plot_value_offset), "%.2f" % getattr(row, score_col),  color = "white", ha = "center", fontsize = 32)
        bars.text(i, getattr(row, score_col) + ((getattr(row, score_col)) / plot_value_offset), "%.3f" % getattr(row, score_col),  color = "black", ha = "center", fontsize = 40)

    
    yerr = [enrich_df[score_col] - enrich_df[ci_lower_bound_col], 
            enrich_df[ci_upper_bound_col] - enrich_df[score_col]]
    
    ax.errorbar(x = enrich_df[bin_col], 
                y = enrich_df[score_col], 
                yerr = yerr, 
                fmt='none', 
                color = "black",
                capsize = 12,
                capthick = 4,
                elinewidth = 4,
                #lw = 6,
                label = "95% CI")
    
    


    ## Plot a line at 0.0 (1.0 after fixing offset) to represent the enrichment 
    ax.plot([-0.5,9.5], [1,1], 'k--', lw = 3, dashes=(5, 7.5))

    ## Update label info
    ax.set_ylabel("Odds ratio of pathogenic versus benign", fontsize = 50)
    ax.set_xlabel("ConSpliceML score bin", fontsize = 50)
    ax.set_title(title, fontsize = 55)

    ax.tick_params(axis="x", labelsize=40, rotation = 90)
    ax.tick_params(axis="y", labelsize=40)

    ax.set_yscale("log")
    
    ax.legend(fontsize = 30)


