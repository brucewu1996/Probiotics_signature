import pandas as pd
import re
import argparse

def format_metaphlan_output_2_species(meta_path) :
    df = pd.read_csv(meta_path,sep='\t',index_col=0)
    idx = [bool(re.search('s__(?!.*t__).*$',x)) for x in df.index]
    df = df.loc[idx,:]
    #idx = [bool(re.search('^(?!k__Archaea|k__Viruses).*$',x))  for x in tmp.index]
    idx = [bool(re.search('k__Bacteria',x)) for x in df.index]
    df = df.loc[idx,:]
    idx = [bool(re.search('unclassified',x)) == False for x in df.index]
    meta_df = df.loc[idx,:]
    meta_df.index = [x.split('|')[-1] for x in meta_df.index]
    return meta_df

def main() :

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="path of input metaphlan output",type = str)
    parser.add_argument("-o", "--output", help="output directory",type = str)

    args = parser.parse_args()
    input_path = args.input
    output_path = args.output

    format_meta_df = format_metaphlan_output_2_species(input_path)
    format_meta_df.to_csv(output_path,sep='\t')

if __name__ == '__main__':
    main()  