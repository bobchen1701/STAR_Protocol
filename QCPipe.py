import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pylab as plt
from scipy import stats
from skimage.filters import threshold_minimum, threshold_mean, threshold_otsu
from sklearn.preprocessing import normalize

def anndata_cleanup(adata_in):
    adata_in.X = adata_in.raw.X.copy()
    for key in adata_in.obs_keys():
        del(adata_in.obs[key])
    for key in adata_in.uns_keys():
        del(adata_in.uns[key])
    for key in adata_in.obsm_keys():
        del(adata_in.obsm[key])
    for key in adata_in.varm_keys():
        del(adata_in.varm[key])
    for key in list(adata_in.layers):
        del(adata_in.layers[key])
    return(adata_in)

def norm_dr(adata,dr='umap',normalization = 'arcsinh', scale=False,target_sum=1,cofactor=1e6,neighborhood=True, n_jobs=16):
    if(normalization == 'arcsinh'):
        print('arcsinh transforming with',cofactor,'cofactor')
        sc.pp.normalize_total(adata,target_sum=target_sum)
        arcsinh_transform(adata,cofactor=cofactor)
    elif(normalization == 'gf-icf'):
        print('gf-icf transforming')
        gf_icf(adata)
        adata.X = adata.layers['gf-icf']
    elif(normalization == 'norm_log1p'):
        print('median normalization and log1p transformation')
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
    elif(normalization == 'none'):
        print('no normalization')
    if(scale):
        print('zscore transform')
        sc.pp.scale(adata)
    sc.pp.pca(adata)
    if(neighborhood):
        print('calculating neighborhood based on sqrt(n_obs) neighbors')
        sc.pp.neighbors(adata,n_neighbors=np.sqrt(adata.n_obs).astype(int))
    if(dr=='umap'):
        print('running umap')
        sc.tl.umap(adata,min_dist=0.25)
    if(dr=='tsne'):
        print('running tsne with',n_jobs,'jobs')
        sc.tl.tsne(adata,perplexity=np.sqrt(adata.n_obs).astype(int),n_jobs=n_jobs)
    
def relative_diversity(adata_in,zscore_minimum = 0, manual_diversity = False, mito_tag = 'MT-'):
    adata_in = adata_in.copy()
    cell_qc_metrics = sc.pp.calculate_qc_metrics(adata_in)[0]
    adata_in.obs['total_counts'] = cell_qc_metrics['total_counts']
    
    adata_in.obs['mito_counts'] = adata_in.X[:,np.where(adata_in.var.index.str.startswith(mito_tag))[0]].sum(axis=1)
    adata_in.obs['mito_percent'] = adata_in.obs['mito_counts']/adata_in.obs['total_counts']

    adata_in.obs['log1p_n_genes_by_counts'] = cell_qc_metrics['log1p_n_genes_by_counts']
    adata_in.obs['log1p_n_genes_by_counts_zscore'] = stats.zscore(adata_in.obs['log1p_n_genes_by_counts'])
    diversity_cutoff = threshold_minimum(adata_in.obs['log1p_n_genes_by_counts_zscore'][adata_in.obs['log1p_n_genes_by_counts_zscore'] >= zscore_minimum])
    
    if(manual_diversity):
        diversity_cutoff = manual_diversity
        adata_in.obs['high_relative_diversity'] = adata_in.obs['log1p_n_genes_by_counts_zscore']>diversity_cutoff
    else:
        adata_in.obs['high_relative_diversity'] = adata_in.obs['log1p_n_genes_by_counts_zscore']>diversity_cutoff
        
    fig = plt.figure(figsize=(4,4))
    ax0 = plt.subplot(111)
    
    hist_plot_div = ax0.hist(adata_in.obs['log1p_n_genes_by_counts_zscore'],bins=np.sqrt(adata_in.n_obs).astype(int))
    ax0.set_ylim([0,max(hist_plot_div[0])])
    
    ax0.axvspan(diversity_cutoff,max(hist_plot_div[1]), alpha=0.25, color='green')
    ax0.vlines(diversity_cutoff,0,max(hist_plot_div[0]),color='green')

    ax0.set_title("Relative Transcript Diversity")
    
    return(adata_in)

def relative_diversity_leiden(adata_in,cluster_hrd_cutoff):
    clust_percent_hrd = adata_in[adata_in.obs['high_relative_diversity']].obs['leiden'].value_counts()/adata_in.obs['leiden'].value_counts()
    clust_percent_hrd[~np.isnan(clust_percent_hrd)] = stats.zscore(clust_percent_hrd[~np.isnan(clust_percent_hrd)])
    adata_in.obs['hrd_leiden'] = np.isin(adata_in.obs['leiden'],clust_percent_hrd[clust_percent_hrd>cluster_hrd_cutoff].index)
    
    return(adata_in,clust_percent_hrd)
    
def filter_by_relative_diversity(adata_in,zscore_minimum = 0, manual_diversity = False, mito_cutoff = 1, mito_tag = 'MT-',threshold_method = 'minimum'):
    #threshold_minimum, threshold_mean, threshold_otsu
    adata_in = adata_in.copy()
    cell_qc_metrics = sc.pp.calculate_qc_metrics(adata_in)[0]
    adata_in.obs['total_counts'] = cell_qc_metrics['total_counts']
    
    adata_in.obs['mito_counts'] = adata_in.X[:,np.where(adata_in.var.index.str.startswith(mito_tag))[0]].sum(axis=1)
    adata_in.obs['mito_percent'] = adata_in.obs['mito_counts']/adata_in.obs['total_counts']

    adata_in.obs['log1p_n_genes_by_counts'] = cell_qc_metrics['log1p_n_genes_by_counts']
    adata_in.obs['log1p_n_genes_by_counts_zscore'] = stats.zscore(adata_in.obs['log1p_n_genes_by_counts'])
    
    if(threshold_method == 'minimum'):
        diversity_cutoff = threshold_minimum(adata_in.obs['log1p_n_genes_by_counts_zscore'][adata_in.obs['log1p_n_genes_by_counts_zscore'] >= zscore_minimum])
    if(threshold_method == 'mean'):
        diversity_cutoff = threshold_mean(adata_in.obs['log1p_n_genes_by_counts_zscore'][adata_in.obs['log1p_n_genes_by_counts_zscore'] >= zscore_minimum])
    if(threshold_method == 'otsu'):
        diversity_cutoff = threshold_otsu(adata_in.obs['log1p_n_genes_by_counts_zscore'][adata_in.obs['log1p_n_genes_by_counts_zscore'] >= zscore_minimum])
    #something here to utilize user input for diversity cutoff
    
    if(manual_diversity):
        diversity_cutoff = manual_diversity
        adata_in.obs['high_relative_diversity'] = adata_in.obs['log1p_n_genes_by_counts_zscore']>diversity_cutoff
    else:
        adata_in.obs['high_relative_diversity'] = adata_in.obs['log1p_n_genes_by_counts_zscore']>diversity_cutoff
    
    subset_adata_in = adata_in[adata_in.obs['high_relative_diversity']].copy()
    
    subset_adata_in.obs['mito_percent_zscore'] = stats.zscore(subset_adata_in.obs['mito_percent'])

    subset_adata_in.obs['low_relative_mito_percent'] = subset_adata_in.obs['mito_percent_zscore']<mito_cutoff
    
    fig = plt.figure(figsize=(8,4))
    ax0 = plt.subplot(121)
    ax1 = plt.subplot(122)
    
    hist_plot_div = ax0.hist(adata_in.obs['log1p_n_genes_by_counts_zscore'],bins=np.sqrt(adata_in.n_obs).astype(int))
    ax0.set_ylim([0,max(hist_plot_div[0])])
    
    ax0.axvspan(diversity_cutoff,max(hist_plot_div[1]), alpha=0.25, color='green')
    ax0.vlines(diversity_cutoff,0,max(hist_plot_div[0]),color='green')

    ax0.set_title("Relative Transcript Diversity")
    
    hist_plot_mito = ax1.hist(subset_adata_in.obs['mito_percent_zscore'],bins=np.sqrt(subset_adata_in.n_obs).astype(int),color='orange')
    ax1.set_ylim([0,max(hist_plot_mito[0])])
    
    ax1.axvspan(min(hist_plot_mito[1]),mito_cutoff,color='green',alpha=0.25)
    ax1.vlines(mito_cutoff,0,max(hist_plot_mito[0]),color='green')
    ax1.set_title("Relative Mito Percent")
    
    subset_adata_in = subset_adata_in[subset_adata_in.obs['low_relative_mito_percent']].copy()
    
    adata_in.obs['high_relative_diversity_and_low_mito'] = np.isin(adata_in.obs.index,subset_adata_in.obs.index)
    
    return(adata_in,subset_adata_in,diversity_cutoff)

def gf_icf(adata, layer=None, norm=None):
    """
    return GF-ICF scores for each element in anndata counts matrix
    Parameters:
        adata (AnnData.AnnData): AnnData object
        layer (str or None): name of layer to perform GF-ICF normalization on. if None, use AnnData.X
        norm (str or None): normalization strategy following GF-ICF transform.
            None: do not normalize GF-ICF scores
            'l1': divide each score by sum of scores for each cell
            'l2': divide each score by sqrt of sum of squares of scores for each cell
    Returns:
        AnnData.AnnData: adata is edited in place to add GF-ICF normalization to .layers
    """
    if layer is None:
        # gene frequency in each cell (l1 norm along rows)
        tf = adata.X / adata.X.sum(axis=1)[:,None]
        # number of cells containing each gene (sum nonzero along columns)
        nt = adata.X.astype(bool).sum(axis=0)
    else:
        # gene frequency in each cell (l1 norm along rows)
        tf = adata.layers[layer] / adata.layers[layer].sum(axis=1)[:,None]
        # number of cells containing each gene (sum nonzero along columns)
        nt = adata.layers[layer].astype(bool).sum(axis=0)
    # inverse cell frequency (total cells / number of cells containing each gene)
    # use "classic" approach here of adding 1 to prevent zero divisions
    idf = np.log((adata.n_obs + 1) / (nt + 1))
    # save normalized GF-ICF scores to .layers
    if norm is None:
        adata.layers["gf-icf"] = tf * idf
    else:
        adata.layers["gf-icf"] = normalize(tf * idf, norm=norm, axis=1)

def find_inflection(ann_data, inflection_percentiles = [0,15,30,100],output_prefix='Output'):
    ann_data_cumsum = np.cumsum(ann_data.obs['total_counts'])
    x_vals=np.arange(0,ann_data.obs.shape[0])
    secant_coef=ann_data_cumsum[ann_data.obs.shape[0]-1]/ann_data.obs.shape[0]
    secant_line=secant_coef*x_vals
    secant_dist=abs(ann_data_cumsum-secant_line)
    inflection_percentiles_inds = np.percentile(x_vals,inflection_percentiles).astype(int)
    inflection_points = secant_dist.argsort()[::-1]
    percentile_points = inflection_points[inflection_percentiles_inds]
    color=plt.cm.tab10(np.linspace(0,1,ann_data.obs.shape[0]))
    plt.figure(figsize=(20,10))
    plt.plot(np.array(ann_data_cumsum), label="Cumulative Sum")
    #plt.plot(np.array(secant_line), label="Secant Line")
    plt.plot(np.array(secant_dist), label="Secant Distance")
    for percentile in percentile_points:
        plt.axvline(x=percentile,ymin=0,c=color[percentile],linestyle='--',linewidth=2,label="Inflection point {}".format(percentile))
    plt.legend()
    #save to file
    if(output_prefix!=''):
        plt.savefig(output_prefix+'_inflectionCheck.png',bbox_inches='tight')
    else:
        plt.show()
    print("Inflection point at {} for {} percentiles of greatest secant distances".format(percentile_points,inflection_percentiles))
    #SJCG: added the ability to return a dictionary of points
    return(dict(zip(inflection_percentiles, percentile_points)))
    
def reorder_AnnData(AnnData, descending = True):
    AnnData.obs['total_counts'] = AnnData.X.sum(axis=1)
    if(descending==True):
        AnnData = AnnData[np.argsort(AnnData.obs['total_counts'])[::-1]].copy()
    elif(descending==False):
        AnnData = AnnData[np.argsort(AnnData.obs['total_counts'])[:]].copy()
    return(AnnData)
    
def arcsinh_transform(AnnData, cofactor = 1000):
    AnnData.X = np.arcsinh(AnnData.X*cofactor,dtype='float')

def cluster_summary_stats(AnnData,raw=False):
    cluster_means = np.zeros((len(np.unique(AnnData.obs['louvain'])),AnnData.n_vars))
    cluster_medians = np.zeros((len(np.unique(AnnData.obs['louvain'])),AnnData.n_vars))
    cluster_stdev = np.zeros((len(np.unique(AnnData.obs['louvain'])),AnnData.n_vars))
    if(raw == True):
        for cluster in range(len(np.unique(AnnData.obs['louvain']))):
            cluster_means[cluster]=np.array(np.mean(AnnData[AnnData.obs['louvain'].isin([str(cluster)])].X,axis = 0))
            cluster_medians[cluster]=np.array(np.median(AnnData[AnnData.obs['louvain'].isin([str(cluster)])].X,axis = 0))
            cluster_stdev[cluster]=np.array(np.std(AnnData[AnnData.obs['louvain'].isin([str(cluster)])].X,axis = 0))
    elif(raw == False):    
        for cluster in range(len(np.unique(AnnData.obs['louvain']))):
            cluster_means[cluster]=np.array(np.mean(AnnData[AnnData.obs['louvain'].isin([str(cluster)])].raw.X,axis = 0))
            cluster_medians[cluster]=np.array(np.median(AnnData[AnnData.obs['louvain'].isin([str(cluster)])].raw.X,axis = 0))
            cluster_stdev[cluster]=np.array(np.std(AnnData[AnnData.obs['louvain'].isin([str(cluster)])].raw.X,axis = 0))
    AnnData.layers['Cluster_Medians'] = np.array(cluster_medians[AnnData.obs['louvain'].astype(int)])
    AnnData.layers['Cluster_Means'] = cluster_means[AnnData.obs['louvain'].astype(int)]
    AnnData.layers['Cluster_Stdevs'] = cluster_stdev[AnnData.obs['louvain'].astype(int)]

def cluster_wilcoxon_rank_sum(AnnData,feature_list,alternative='greater'):
    cluster_list = AnnData.obs['louvain']
    p_values = np.zeros((len(np.unique(cluster_list)),len(feature_list)))
    for cluster in range(len(np.unique(cluster_list))):
        for feature in range(len(feature_list)):
            p_values[cluster,feature]=stats.mannwhitneyu(AnnData[cluster_list.isin([str(cluster)])].obs_vector(feature_list[feature]),AnnData.obs_vector(feature_list[feature]),alternative='greater',use_continuity=True)[1]
    AnnData.uns['Cluster_p_values'] = pd.DataFrame(p_values,np.arange(len(np.unique(cluster_list))),feature_list)

def cluster_p_threshold(AnnData,threshold = 0.05):
    for columns in AnnData.uns['Cluster_p_values']:
        AnnData.obs[columns+'_threshold'] = (AnnData.uns['Cluster_p_values']<threshold)[columns][AnnData.obs['louvain'].astype(int)].astype(int).values
        AnnData.obs[columns+'_enrichment'] = (AnnData.uns['Cluster_p_values'])[columns][AnnData.obs['louvain'].astype(int)].values