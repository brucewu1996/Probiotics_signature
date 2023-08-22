import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from math import log10
import seaborn as sns # type: ignore
import plotly.express as px # type: ignore
import plotly.io as pio # type: ignore
import os
import argparse


def signature_coefficient_relative_importance_barplot(sig_coef_matrix,output_path,output_prefix,xlabel = 'Relative importantance'):
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
               color_discrete_sequence=px.colors.qualitative.Set2,category_orders={"variable": legend_order},labels={'value' : xlabel})
   fig.update_layout(barmode='stack', xaxis={'categoryorder':'array','categoryarray':bar_order}, title_x=0.5,yaxis_title=None)
   pio.write_image(fig,output_path + output_prefix + '_' +'signature_composition_barplot'  + '.png',format = 'png',scale = 2)
   pio.write_image(fig,output_path + output_prefix + '_' +'signature_composition_barplot'  + '.svg',format = 'svg',scale = 2)
   #fig.show() 
   
def main() :
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",help="path of signature proportion matrix and coefficient matrix")
    parser.add_argument("-l", "--lacto_prefix",help="prefix of lacto signature proportion matrix")
    parser.add_argument("-b", "--bifido_prefix",help="prefix of bifido signature proportion matrix")
    parser.add_argument("-o","--output",help="output path of probiotics signature")
    args = parser.parse_args()
    
    if os.path.isdir(args.output) == False :
        os.mkdir(args.output)
    #### sig coefficient pie chart
    lacto_coef = pd.read_csv(args.input + args.lacto_prefix +'_signature_coefficient_matrix.txt',sep = '\t',index_col = 0)
    bifido_coef = pd.read_csv(args.input + args.bifido_prefix +'_signature_coefficient_matrix.txt',sep = '\t',index_col = 0)
    
    signature_coefficient_relative_importance_barplot(lacto_coef,args.output,args.lacto_prefix)
    signature_coefficient_relative_importance_barplot(bifido_coef,args.output,args.bifido_prefix)
        
if __name__ == '__main__' :
    main()