import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from math import log10
import seaborn as sns # type: ignore
import os
import argparse

def signature_proportion_heatmap(sig_matrix,label_df,order,legend_prefix,output,format = 'png') :
    
    metadata_df = label_df.sort_values(by=order)
    plot_df = sig_matrix.loc[:,metadata_df.index].T
    
    cmap = sns.light_palette("darksalmon", as_cmap=True)
    #cmap = sns.diverging_palette(240, 10, n=9,as_cmap=True)
    cluster1 = metadata_df[order[0]].values
    cluster1_lut = dict(zip( set(cluster1),  [cm.Set2(x) for x in range(10)] )) # type: ignore
    cluster1_colors = pd.Series(cluster1).map(cluster1_lut)
    cluster1_colors.index = metadata_df.index
    
    if len(order) == 1 :
    
        g = sns.clustermap(plot_df, cmap=cmap,vmin=0, vmax=1,  cbar_kws={"shrink": .3},
                    row_cluster=False, 
                    col_cluster=False,
                    row_colors=[cluster1_colors[plot_df.index]],
                    linewidths=0, figsize=(6,12))
        g.ax_cbar.set_position((0.05, .3, .03, .4))# type: ignore
        g.cax.set_title("Signature proportion",fontsize = 10)# type: ignore

        for label in sorted(set(cluster1)):
            g.ax_col_dendrogram.bar(0, 0, color=cluster1_lut[label], label=label, linewidth=0.5)
        l1 = g.ax_col_dendrogram.legend(title=legend_prefix, loc="center", ncol=5, bbox_to_anchor=(0.47, .9), bbox_transform=plt.gcf().transFigure) # type: ignore
    else :

        cluster2 = metadata_df.loc[:,order[1]].values
        cluster2_lut = dict(zip( set(cluster2),  [cm.Set3(x) for x in range(10)] ))# type: ignore
        cluster2_colors = pd.Series(cluster2).map(cluster2_lut)
        cluster2_colors.index = metadata_df.index

        cmap = sns.light_palette("darksalmon", as_cmap=True)
        #cmap = sns.diverging_palette(240, 10, n=9,as_cmap=True)

        g = sns.clustermap(plot_df, cmap=cmap,vmin=0, vmax=1,  cbar_kws={"shrink": .3},
                        row_cluster=False, 
                        col_cluster=False,
                        row_colors=[cluster1_colors[plot_df.index],cluster2_colors[plot_df.index]],
                        linewidths=0, figsize=(6,10))
        g.ax_cbar.set_position((0.05, .3, .03, .4))# type: ignore
        g.cax.set_title("Signature proportion",fontsize = 10)# type: ignore

        
        for label in sorted(set(cluster1)):
            g.ax_col_dendrogram.bar(0, 0, color=cluster1_lut[label], label=label, linewidth=0)
        l1 = g.ax_col_dendrogram.legend(title=order[0] + ' ' + legend_prefix, loc="center", ncol=5, bbox_to_anchor=(0.47, 1.0), bbox_transform=plt.gcf().transFigure) # type: ignore
    
        for label in sorted(set(cluster2)):
            g.ax_row_dendrogram.bar(0, 0, color=cluster2_lut[label], label=label, linewidth=0)
        l2 = g.ax_row_dendrogram.legend(title=order[1] + ' ' + legend_prefix, loc="center", ncol=5, bbox_to_anchor=(0.47, .9), bbox_transform=plt.gcf().transFigure) # type: ignore

    plt.savefig(output,dpi = 300,bbox_inches = 'tight',format = format)

def main() :
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",help="path of signature proportion matrix and coefficient matrix")
    parser.add_argument("-c", "--cluster",help="path of consensus cluster label")
    parser.add_argument("-l", "--lacto_prefix",help="prefix of lacto signature proportion matrix")
    parser.add_argument("-b", "--bifido_prefix",help="prefix of bifido signature proportion matrix")
    parser.add_argument("-o","--output",help="output path of probiotics signature")
    args = parser.parse_args()
    
    if os.path.isdir(args.output) == False :
        os.mkdir(args.output)
        
    ## merge lacto & bifido signature proportion matrix
    lacto_matrix = pd.read_csv(args.input + args.lacto_prefix + '_signature_proportion_matrix.txt',sep = '\t',index_col = 0)
    bifido_matrix = pd.read_csv(args.input + args.bifido_prefix + '_signature_proportion_matrix.txt',sep = '\t',index_col = 0)
    sig_matrix = pd.concat([lacto_matrix,bifido_matrix],axis=0)
    sig_matrix.index = [x.capitalize() for x in sig_matrix.index]  # type: ignore
    ## 
    lacto_cluster = pd.read_csv(args.cluster + 'consensus_clustering_label_4cluster.txt',sep = '\t',index_col = 0)
    lacto_cluster.columns = ['Lacto']
    bifido_cluster = pd.read_csv(args.cluster + 'consensus_clustering_label_4cluster.txt',sep = '\t',index_col = 0)
    #bifido_cluster = pd.read_csv(args.cluster + 'bifido_clustering.txt',sep = '\t',index_col = 0)

    metadata = pd.concat([lacto_cluster,bifido_cluster],axis=1)
    metadata.columns = ['Lactobacillus','Bifidobacterium']  # type: ignore
    
    signature_proportion_heatmap(sig_matrix.T,metadata,['Lactobacillus'],"consensus cluster",args.output + 'sig_proportion_heatmap_ordr_by_lacto.png')
        
if __name__ == '__main__' :
    main()