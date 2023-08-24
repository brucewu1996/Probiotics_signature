import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from math import log10
from matplotlib.pyplot import gcf
from statannot import add_stat_annotation
from itertools import combinations
import plotly.express as px
import plotly.io as pio
import seaborn as sns
import os
import argparse


def signature_proportion_heatmap(sig_matrix,label_df,order,legend_prefix,output,format = 'png') :
    
    metadata_df = label_df.sort_values(by=order)
    plot_df = sig_matrix.loc[:,metadata_df.index]
    
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
                    col_colors=[cluster1_colors[plot_df.columns]],
                    linewidths=0, figsize=(10, 6))
        g.ax_cbar.set_position((0.05, .3, .03, .4))# type: ignore
        g.cax.set_title("Signature proportion",fontsize = 10)# type: ignore

        for label in set(cluster1):
            g.ax_col_dendrogram.bar(0, 0, color=cluster1_lut[label], label=label, linewidth=0)
        l1 = g.ax_col_dendrogram.legend(title=legend_prefix, loc="center", ncol=5, bbox_to_anchor=(0.47, .9), bbox_transform=gcf().transFigure) # type: ignore
    else :

        cluster2 = metadata_df[order[1]].values
        cluster2_lut = dict(zip( set(cluster2),  [cm.Set3(x) for x in range(10)] ))# type: ignore
        cluster2_colors = pd.Series(cluster2).map(cluster2_lut)
        cluster2_colors.index = metadata_df.index

        cmap = sns.light_palette("darksalmon", as_cmap=True)
        #cmap = sns.diverging_palette(240, 10, n=9,as_cmap=True)

        g = sns.clustermap(plot_df, cmap=cmap,vmin=0, vmax=1,  cbar_kws={"shrink": .3},
                        row_cluster=False, 
                        col_cluster=False,
                        col_colors=[cluster1_colors[plot_df.columns],cluster2_colors[plot_df.columns]],
                        linewidths=0, figsize=(10, 6))
        g.ax_cbar.set_position((0.05, .3, .03, .4))# type: ignore
        g.cax.set_title("Signature proportion",fontsize = 10)# type: ignore

        for label in set(cluster1):
            g.ax_col_dendrogram.bar(0, 0, color=cluster1_lut[label], label=label, linewidth=0)
        l1 = g.ax_col_dendrogram.legend(title=order[0] + ' ' + legend_prefix, loc="center", ncol=5, bbox_to_anchor=(0.47, 1.0), bbox_transform=gcf().transFigure) # type: ignore

        for label in set(cluster2):
            g.ax_row_dendrogram.bar(0, 0, color=cluster2_lut[label], label=label, linewidth=0)
        l2 = g.ax_row_dendrogram.legend(title=order[1] + ' ' + legend_prefix, loc="center", ncol=5, bbox_to_anchor=(0.47, .9), bbox_transform=gcf().transFigure) # type: ignore

    plt.savefig(output,dpi = 300,bbox_inches = 'tight',format = format)
    
    
def signature_composition_piechart(input_df,output) :

    df = input_df.T.copy()
    df['Species'] = df.index
    df = df.melt(id_vars='Species')

    for prefix in np.unique(df['variable']) :
        plot_df = df.loc[np.where(df['variable'] == prefix,True,False),:]
        fig = px.pie(plot_df, values='value', names='Species', title='Signature composition : ' + prefix)
        fig.update_traces(textposition='inside',sort=False)
        fig.update_layout(uniformtext_minsize=20, uniformtext_mode=False)
        fig.write_image(output + prefix.replace(' ','_') + '.png')  # type: ignore
        
        
def signature_coefficient_relative_importance_stack_barplot(sig_coef_matrix,output_path,output_prefix,xlabel = 'Relative importantance'):
   sig_coef_matrix = sig_coef_matrix.T
   sig_coef_matrix.columns = [x.capitalize() for x in sig_coef_matrix.columns]
   plot_df = pd.DataFrame()
   plot_df['Subtype'] = [x.replace('_',' ') for x in sig_coef_matrix.index]
   for c in sig_coef_matrix.columns :
      relative_coef = sig_coef_matrix[c].transform(lambda x: 100 * x/x.sum()).values
      log_relative_coef = list(map(lambda x : log10(x+1),relative_coef))
      relative_importance = [100 * x / sum(log_relative_coef) for x in log_relative_coef]
      plot_df[c] = relative_importance
   plot_df = plot_df.melt(id_vars='Subtype')
   plot_df['Subtype'] = [x.capitalize() for x in plot_df['Subtype'].values]  # type: ignore
   
   legend_order = np.unique(plot_df['variable'])
   bar_order = np.unique(plot_df['Subtype'])
   fig = px.bar(plot_df, x="value", y="variable", color="Subtype",width=1200, height=500,template='simple_white',
               color_discrete_sequence=px.colors.qualitative.Pastel,category_orders={"variable": legend_order},labels={'value' : xlabel})
   fig.update_layout(barmode='stack', xaxis={'categoryorder':'array','categoryarray':bar_order}, title_x=0.5,yaxis_title=None)
   pio.write_image(fig,output_path + output_prefix + '_' +'signature_composition_barplot'  + '.png',format = 'png',scale = 2)
   pio.write_image(fig,output_path + output_prefix + '_' +'signature_composition_barplot'  + '.svg',format = 'svg',scale = 2)
   fig.show() 

def signature_coefficient_relative_importance_barplot(sig_matrix,output_path,xlabel = 'Species',xticks = None,format = 'png',cmap = 'Set3',fig_size = (7,3)) :
    '''
    sig_matrix : dataframe; x is signature , y is element in signature
    output_path : str ; output folder of fig output
    '''
    # convert sig matrix coef into relative importance
    plot_df = sig_matrix.T.copy()
    for c in plot_df.columns :
        relative_coef = plot_df[c].transform(lambda x: 100 * x/x.sum()).values
        log_relative_coef = list(map(lambda x : log10(x+1),relative_coef))
        plot_df[c] = log_relative_coef
    # plot setting
    n_element = plot_df.shape[0]
    plot_df['Species'] = list(plot_df.index)
    # plot loliplot
    for c in plot_df.columns[:-1] :
        plt.figure(figsize=fig_size)
        sns.barplot(data=plot_df,y=c,x='Species',palette=cmap)
        plt.title(c.capitalize())
        if c != plot_df.columns[-2] :
            plt.xticks([])
            plt.xlabel('')
            plt.ylabel('')
        else :
            if not xticks :
                xticks = list(plot_df.index)
                #xlabel = [x.split('_')[0][0] + '.' + x.split('_')[1] for x in plot_df.index]
            x = np.arange(n_element)
            plt.xticks(x,xticks,rotation=60)
            plt.xlabel(xlabel)
            plt.ylabel('')
            #plt.ylabel("Relative importance")
        plt.savefig("%s%s_sig_relative_importance_barplot.%s" % (output_path,c.replace(' ','_'),format),dpi=300,format = 'svg')
        plt.show()

def sig_proportion_compoarison(sig_matrix,metadata,output_path,format='png') :

    sig_matrix_proportion = pd.concat([sig_matrix.T,metadata],axis=1).melt(id_vars='cluster')
    sig_matrix_proportion.columns = ['Cluster','Signature','Proportion']  # type: ignore

    fig, axes = plt.subplots(2,4,figsize = (20,12))
    for idx,ax in enumerate(axes.ravel()) : # type: ignore
        
        sig_idx = np.where(sig_matrix_proportion['Signature'] == sig_matrix.index[idx],True,False)
        df = sig_matrix_proportion.loc[sig_idx,:] 
        cluster = df['Cluster'].unique()
        box_pair = list(combinations(cluster,2))  # type: ignore
        sns.boxplot(data=df, x='Cluster', y= 'Proportion',palette='rainbow_r',ax = ax)
        ax.set_title(sig_matrix.index[idx])
        ax.set_xlabel('Cluster')
        ax.set_ylabel('Signature proportion')
        add_stat_annotation(ax = ax,data=df, x='Cluster', y='Proportion',box_pairs=box_pair,test='Mann-Whitney', text_format='star', loc='inside', verbose=2)
    plt.tight_layout()
    plt.savefig(output_path,dpi = 300,bbox_inches = 'tight',format=format)
    

def main() :
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",help="path of signature proportion matrix and coefficient matrix")
    parser.add_argument("-l","--lacto-prefix",help="prefix of signature matrix")
    parser.add_argument("-l","--lacto-prefix",help="prefix of signature matrix")
    parser.add_argument("-c", "--cluster",help="label of sample cluster")
    parser.add_argument("-o","--output",help="output path of probiotics signature")
    args = parser.parse_args()
    ## merge lacto & bifido signature proportion matrix
    lacto_matrix = pd.read_csv(args.input + '%s_signature_proportion_matrix.txt',sep = '\t',index_col = 0)
    bifido_matrix = pd.read_csv(args.input + '%s_signature_proportion_matrix.txt',sep = '\t',index_col = 0)
    sig_matrix = pd.concat([lacto_matrix,bifido_matrix],axis=0)
    sig_matrix.index = [x.capitalize() for x in sig_matrix.index]  # type: ignore
    ## 
    lacto_cluster = pd.read_csv(args.cluster + 'lacto_sig_cluster.txt',sep = '\t',index_col = 0)
    bifido_cluster = pd.read_csv(args.cluster + 'bifido_sig_cluster.txt',sep = '\t',index_col = 0)
    metadata = pd.concat([lacto_cluster,bifido_cluster],axis=1)
    metadata.columns = ['Lactobacillus','Bifidobacterium']  # type: ignore
    
    if os.path.isdir(args.output) == False :
        os.mkdir(args.output)
        
    signature_proportion_heatmap(sig_matrix,metadata,['Lactobacillus'],"consensus cluster",args.output + 'sig_proportion_heatmap_ordr_by_lacto.png')
    #### sig coefficient pie chart
    lacto_coef = pd.read_csv(args.input + 'lactobacillus_signature_coefficient_matrix.txt',sep = '\t',index_col = 0)
    bifido_coef = pd.read_csv(args.input + 'bifidobacterium_signature_coefficient_matrix.txt',sep = '\t',index_col = 0)
    
    signature_composition_piechart(lacto_coef,args.output)
    signature_composition_piechart(bifido_coef,args.output)
    
        
if __name__ == '__main__' :
    main()