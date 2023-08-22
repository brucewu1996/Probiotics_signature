import pandas as pd
import numpy as np
import re,argparse

def metaphlan_subtype(meta_df,subtype_df,genus,species_col,group_col) :

    subtype_dict = dict()
    for i in range(subtype_df.shape[0]) :
        s = subtype_df[species_col][i]
        s = 's__' +s.replace(' ','_')
        type = subtype_df[group_col][i]
        if type == 'no phylogroup' :
            type = genus + '_others'
        else :
            type = type + '_' + 'subtype'
        
        subtype_dict[s] = type  

    idx = [bool(re.search(genus,x)) for x in meta_df.index]
    df = meta_df.loc[idx,:]

    idx = list()
    ## change species name to subtype name
    convert_dict = dict()
    for i in range(df.shape[0]) :
        if df.index[i] in subtype_dict.keys() :
            s = df.index[i]
            subtype = subtype_dict[s]
        elif df.index[i] == 's__Lactobacillus_casei_group' :
            s = df.index[i]
            subtype = 'Lacticaseibacillus_subtype'
        else :
            s = df.index[i]
            subtype = genus + '_others'

        if subtype in convert_dict.keys() :
            convert_dict[subtype] = convert_dict[subtype] + [s]
        else :
            convert_dict[subtype] = [s]
        idx.append(subtype)

    df.index = idx
    meta_df = df.copy()
    meta_df['subtype'] = df.index
    meta_df = meta_df.groupby('subtype').agg('sum')
    meta_df.index.name = None

    return meta_df,convert_dict

def main() :
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",help="path of subtype matrix")
    parser.add_argument("-r","--reference",help="path of subtype reference. ex: s__Lactobacillus_plantarum : Lactiplantibacillus_subtype")
    parser.add_argument("-f","--format",type=bool,default=True,help="species name is metaphlan format or not! metaphlan format : s__Lactobacillus_plantarum")
    parser.add_argument("-o","--output",help="output path of subtype matrix")
    args = parser.parse_args()
    
    matrix = pd.read_csv(args.input,sep='\t',index_col=0)
    subtype = pd.read_csv(args.reference)
    if args.format == False :
        matrix.index = ['s__' + x.replace(' ','_') for x in matrix.index]# type: ignore   
        
    subtype_matrix,_ = metaphlan_subtype(matrix,subtype,'Lactobacillus','species','phylogroup')
    subtype_matrix.to_csv(args.output,sep='\t')
 
if __name__ == '__main__' :
    main()
