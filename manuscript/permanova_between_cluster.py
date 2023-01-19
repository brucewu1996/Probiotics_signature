import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm,colors
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
import seaborn as sns# type: ignore
import os
import argparse
from itertools import combinations
from skbio.stats.ordination import pcoa # type: ignore
from skbio.stats.distance import permanova,DistanceMatrix# type: ignore

def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    **kwargs
        Forwarded to `~matplotlib.patches.Ellipse`

    Returns
    -------
    matplotlib.patches.Ellipse
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1]) # type: ignore
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)# type: ignore
    ell_radius_y = np.sqrt(1 - pearson)# type: ignore
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, **kwargs)
    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std# type: ignore
    mean_x = np.mean(x)

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std# type: ignore
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)

def pcoa_scatterplot(dmatrix,metadata,hue,title,output_path,fig_size = (8,6),style=None) :
    '''
    dmatrix : dataframe, distance matrix with index & columns
    metadata : dataframe, row as sample ,columns as feature
    '''
    x = dmatrix.to_numpy()
    pcoa_r = pcoa(x,number_of_dimensions=2)
    pcoa_df = pcoa_r.samples
    pcoa_df.index = dmatrix.index
    pc1_exp,pc2_exp = round(100*pcoa_r.proportion_explained,2)
    pcoa_df = pd.concat([pcoa_df,metadata],axis=1)

    plt.figure(figsize=fig_size)
    if style != None :
        sns.scatterplot(data=pcoa_df,x= 'PC1',y='PC2',hue=hue,style='Diagnosis',palette="Set2",s=100)
    else :
        sns.scatterplot(data=pcoa_df,x= 'PC1',y='PC2',hue=hue,palette="Set2")

    xlabel = "PCoA1 (" + str(pc1_exp) + "%)"
    ylabel = "PCoA2 (" + str(pc2_exp) + "%)"
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(bbox_to_anchor=(1.1, 1))
    plt.title(title)
    plt.savefig(output_path,dpi = 300,bbox_inches = 'tight')

def mds_scatterplot(dmatrix,metadata,hue,title,output_path,style=None) :
    '''
    dmatrix : dataframe, distance matrix with index & columns
    metadata : dataframe, row as sample ,columns as feature
    '''
    x = dmatrix.to_numpy()
    mds = MDS(n_components=2,dissimilarity='precomputed')# type: ignore
    mds_r = mds.fit_transform(x)
    mds_df = pd.DataFrame(mds_r,index = dmatrix.index)
    mds_df.columns = ['MDS1','MDS2']# type: ignore
    mds_df = pd.concat([mds_df,metadata],axis=1)

    plt.figure(figsize=(8,6))
    if style != None :
        sns.scatterplot(data=mds_df,x= 'MDS1',y='MDS2',hue=hue,style='Diagnosis',palette="Set2",s=100)
    else :
        sns.scatterplot(data=mds_df,x= 'MDS1',y='MDS2',hue=hue,palette="Set2")

    plt.legend(bbox_to_anchor=(1.1, 1))
    plt.title(title)
    plt.savefig(output_path,dpi = 300,bbox_inches = 'tight')

def pcoa_with_permanova_scatterplot(dmatrix,metadata,hue,condition,title,output_path,color = ['#66c2a5','#fc8d62'],fig_size = (8,6),format = 'png') :
    '''
    dmatrix : dataframe, distance matrix with index & columns
    metadata : dataframe, row as sample ,columns as feature
    '''
    x = dmatrix.to_numpy()
    pcoa_r = pcoa(x,number_of_dimensions=2)
    pcoa_df = pcoa_r.samples
    pcoa_df.columns = ['PC1','PC2']
    pcoa_df.index = dmatrix.index
    pc1_exp,pc2_exp = round(100*pcoa_r.proportion_explained,2)
    ###centroid
    #centroid for first condition
    c1_idx = np.where(metadata[hue] == condition[0],True,False)
    x1= pcoa_df.loc[c1_idx,'PC1'].values
    y1 = pcoa_df.loc[c1_idx,'PC2'].values
    c2_idx = np.where(metadata[hue] == condition[1],True,False)
    x2= pcoa_df.loc[c2_idx,'PC1'].values
    y2 = pcoa_df.loc[c2_idx,'PC2'].values
    ###permanova
    m = DistanceMatrix(dmatrix,ids=dmatrix.index)
    perm_r = permanova(m,metadata,hue,permutations= 10000)
    pvalue = perm_r['p-value']
    text_pos = pcoa_df.mean().values
    text = 'Permanova p-value : ' + str(round(pvalue,3))
    #plot 
    pcoa_df = pd.concat([pcoa_df,metadata],axis=1)
    plt.figure(figsize=fig_size)
    ax = plt.gca()
    sns.scatterplot(data=pcoa_df,x= 'PC1',y='PC2',hue=hue,palette=color)
    text_x = min(pcoa_df['PC1'])
    text_y = max(pcoa_df['PC2'])
    plt.text(text_x,text_y,text)
    
    #plot centroid & ellipse for condition 1
    confidence_ellipse(x1, y1, ax, edgecolor=color[0],n_std=3)
    ax.scatter(np.mean(x1), np.mean(y1), c=color[0], s=100)
    #plot centroid & ellipse for condition 2
    confidence_ellipse(x2, y2, ax, edgecolor=color[1],n_std=3)
    ax.scatter(np.mean(x2), np.mean(y2), c=color[1], s=100)

    xlabel = "PCoA1 (" + str(pc1_exp) + "%)"
    ylabel = "PCoA2 (" + str(pc2_exp) + "%)"
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.savefig(output_path,dpi = 300,bbox_inches = 'tight',format = format)


def main() :
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",help="path of unifrac diatance matrix")
    parser.add_argument("-c", "--cluster",help="path of consensus clustering label")
    parser.add_argument("-o","--output",help="output path of probiotics signature")
    parser.add_argument("-e","--hue",help='hue by which label')
    args = parser.parse_args()
    #set color of each cluster
    cluster_color = dict()
    cmap = cm.get_cmap('Set2', 8) 
    n = 1
    for i in range(cmap.N):    # type: ignore
        if n < 6 :
            rgba = cmap(i)
            # rgb2hex accepts rgb or rgba
            color = colors.rgb2hex(rgba)
            cluster_color[n] = color
            n+=1
        else : 
            break

    if os.path.isdir(args.output) == False :
        os.mkdir(args.output)
    ### import metadata 
    lacto_cluster = pd.read_csv(args.cluster + 'lacto_clustering.txt',sep = '\t',index_col = 0)
    bifido_cluster = pd.read_csv(args.cluster + 'bifido_clustering.txt',sep = '\t',index_col = 0)
    plot_metadata = pd.concat([lacto_cluster,bifido_cluster],axis=1)
    plot_metadata.columns = ['Lactobacillus','Bifidobacterium']# type: ignore
    ### import unifrac distance
    matrix = pd.read_csv(args.input,sep = '\t',index_col=0)
    matrix.columns = [x.split('_')[0] for x in matrix.index]# type: ignore
    matrix.index = [x.split('_')[0] for x in matrix.index]# type: ignore
    
    matrix = matrix.loc[plot_metadata.index,plot_metadata.index]
    # Load the pandas matrix into skbio format
    dm = DistanceMatrix(matrix)
    # plot permanova
    hue = args.hue
    cb = list(combinations(np.unique(plot_metadata[args.hue].values),2))
    prefix = args.hue
    plot_matrix = matrix.copy()

    for i in range(len(cb)) :
        condition = cb[i]
        color = [cluster_color[x] for x in condition]
        idx = np.where((plot_metadata[hue] == condition[0]) | (plot_metadata[hue] == condition[1]),True,False)
        m = plot_matrix.loc[idx,idx]
        label = plot_metadata.loc[idx,:]
        #title = "Permanova with weighted unifrac distance between " + prefix +" cluster"+ str(condition[0]) + ' and '+prefix + ' cluster' + str(condition[1])
        title = ''
        file = args.output + 'permanova_' + prefix + '_' + 'c' +str(condition[0]) + '_c' + str(condition[1]) + '.svg' 
        try :
            pcoa_with_permanova_scatterplot(m,label,hue=hue,condition=condition,color= color,title=title,output_path=file,format='svg')
        except :
            pass

        
if __name__ == '__main__' :
    main()