import pandas as pd
import numpy as np
import re
import glob
import os
from itertools import compress
import matplotlib.pyplot as plt

def split_tax(mpa_report) :
    '''
    tax is a row name of metaphlan format taxonomy information
    ex : 'k__Bacteria|p__Actinobacteria|c__Actinomycetia|
    o__Streptomycetales|f__Streptomycetaceae|
    g__Streptomyces|s__Streptomyces_griseus_group'
    '''
    data = pd.read_csv(mpa_report,sep = "\t",index_col =0)
    tax = list(data.index)
    
    levels = ["Kingdom","Phylum","Class","Order","Family","Genus","Species"]
    levels_dict = {"k" : 0,"p" : 1,"c" : 2,"o" : 3,"f" : 4,"g" : 5,"s" : 6}
    tax_table = pd.DataFrame(columns = levels)
    tax_idx = list()
    for i in range(len(tax)) :
        tax_name = tax[i].split("|",7)
        tax_table.loc[i] = [np.nan] * 7
        tax_idx.append(tax_name[-1])
        for j in range(len(tax_name)) :
            name = tax_name[j]
            level = name.split("_")[0]
            n = name.count('_')
            if n > 2  and level != 's' :
                t = name.split("__",2)[1]
                tax_table.iloc[i,5] = t.split("_",2)[0]
                tax_table.iloc[i,6] = t
            else :
                idx = levels_dict[level]
                tax_table.iloc[i,idx] = name.split("__",2)[1]

    tax_table.index = tax_idx
    data.index = tax_idx
    return data,tax_table

def taxa_prevalence(otu_table,threshold = 0) :
    prevalence = list()
    for i in range(otu_table.shape[0]) :
        prevalence.append(np.mean(otu_table.iloc[i,:] > threshold)*100 )
        
    return prevalence

def filter_taxa_by_prevalence(ab_table,tax_table,prevalence_threshold) :
    '''
    if relative abundance is necessary, convert to ra must be complete before filter otu table.
    '''
    tax_idx = list(tax_table.index)
    #remove taxonomy which is unidentified or uncultured
    tmp = [bool(re.search("unidentified|uncultured",i)) for i in tax_idx]
    rm_idx = list(tax_table.iloc[tmp,:].index)
    tax_table = tax_table.drop(rm_idx)
    ab_table = ab_table.drop(rm_idx)
    #keep taxonomy only belong Bacteria
    tax_table[tax_table["Kingdom"] == "Bacteria"]
    ab_table[tax_table["Kingdom"] == "Bacteria"]
    #remove low prevalence(%) taxa 
    pre = taxa_prevalence(ab_table)
    rm_idx = [i > prevalence_threshold for i in pre]
    ab_table = ab_table.iloc[rm_idx,:]
    tax_table = tax_table.iloc[rm_idx,:]
    
    return ab_table,tax_table

def relative_abundance(otu_table) :
    
    df = pd.DataFrame()
    for i in otu_table.columns :
        df[i] = (otu_table[i] / otu_table[i]['k__Bacteria'])*100
        
    return df
    
def taxa_abundance(otu_table) :
    n_feature = otu_table.shape[0]
    abundance = np.zeros(n_feature)
    for i in range(n_feature) :
        abundance[i] = np.mean(otu_table.iloc[i,:])
    return abundance
    

def add_taxa_name(otu_table,tax_table):
    taxa_name = list()
    level_dict = {1:"k",2:"p",3:"c",4:"o",5:"f",6:"g"}
    for i in range(tax_table.shape[0]):
        tmp = list(tax_table.iloc[i,:])
        try:
            na_idx = tmp.index(np.NaN)
            name = str(level_dict[na_idx]) + "__" + str(tmp[na_idx-1])
        except ValueError:
            na_idx = -1
            name = "g__" + str(tmp[-1])
        taxa_name.append(name)
        
    return taxa_name

def binary_anomaly_detection_scatter_plot(threshold,score,label,class_label,path,title):
        """ 
        DESCRIPTION
        Plot the curve of distance
        --------------------------------------------------------------- 
        label parameter only for binary data
        class_label[0] equal to health data
        class_label[1] equal to non-health data
        """ 
        if len(class_label) == 1 :
            binary = False
        else :
            binary = True

        n = len(score)
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(1, 1, 1)
        radius = np.ones((n, 1)) * threshold

        score_dict = {"score" : score,"label" : label}
        df = pd.DataFrame(score_dict)

        if binary :
            df1 = df[df["label"] == class_label[0]]
            df2 = df[df["label"] == class_label[1]]
            x1 = np.arange(1,df1.shape[0]+1)
            x2 = np.arange(df1.shape[0],df.shape[0])
        else :
            df1 = df[df["label"] == class_label[0]]
            x1 = np.arange(1,df1.shape[0]+1)

        ax.plot(radius, 
                color='r',
                linestyle='-', 
                marker='None',
                linewidth=3, 
                markeredgecolor='k',
                markerfacecolor='w', 
                markersize=6)
        
        ax.plot(x1,
                df1["score"].values,
                color='k',
                linestyle=':',
                marker='o',
                linewidth=1,
                markeredgecolor='k',
                markerfacecolor='C3',
                markersize=6,
                label = class_label[0])
        if binary :
            ax.plot(x2,
                    df2["score"].values,
                    color='k',
                    linestyle=':',
                    marker='o',
                    linewidth=1,
                    markeredgecolor='k',
                    markerfacecolor='C4',
                    markersize=6,
                    label = class_label[1])

            ax.legend(["Threshold",class_label[0],class_label[1]], 
                  ncol=1, loc=0, 
                  edgecolor='black', 
                  markerscale=1, fancybox=True)
        else :
            ax.legend(["Threshold",class_label[0]], 
                  ncol=1, loc=0, 
                  edgecolor='black', 
                  markerscale=1, fancybox=True)
        
        ax.set_xlabel('Samples')
        ax.set_ylabel('Score')
        plt.title(title)
        
        ax.yaxis.grid()
        plt.savefig(path,dpi = 300)
        plt.show()

def merge_emu_output(path) :
    '''
    path : folder path of emu output
    sample format : sample + sample_barcode + .fastq
    emu output report format : sample name +'rel-abundance.tsv'
    '''
    #path = '/home/bruce1996/data/Yi-Fung-Chuang/within_individual/emu_output/'
    file_list = glob.glob(path +"*tsv")
    file_list.sort()
    for idx,file in enumerate(file_list) :
        _,filename = os.path.split(file)
        if re.search('threshold',filename) != None :
            continue
        sample_name = filename.split('_')[1] 
        df = pd.read_csv(file,sep = '\t')
        #only extract abundance & species columns
        df = df.iloc[:,0:2]
        df.columns = [sample_name,'species']
        if idx == 0 :
            merge_df = df
        else :
            merge_df = pd.merge(merge_df,df,how = 'outer',on = 'species')
    merge_df = merge_df.fillna(0)
    idx = [1,0] + [x for x in range(2,merge_df.shape[1])]
    merge_df = merge_df.iloc[:,idx]
    #merge_df.index = merge_df['species']

    return merge_df

