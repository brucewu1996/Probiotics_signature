import numpy as np
import pandas as pd
import argparse
import os
from sklearn.decomposition import NMF# type: ignore

def finger_print_proportion(x,w,h):
    n_finger_print = w.shape[1]
    n_sample = w.shape[0]
    proportion_matrix = np.zeros([n_sample,n_finger_print+1])

    for i in range(n_sample) :
        total = sum(x[i,:])
        accum = 0
        if total == 0 :
            pass
        else :
            for j in range(n_finger_print) :
                ab = sum(np.dot(w[i,j],h[j,:]))# type: ignore                
                proportion_matrix[i,j] = ab / total
                accum += (ab / total)

        proportion_matrix[i,-1] = 1 - accum

    return proportion_matrix

def fit_signature(denovo_df,fit_df,k,prefix) :
    X = denovo_df.T.to_numpy()
    x = fit_df.T.to_numpy()
    nmf_model = NMF(n_components=k, init='random', random_state=0)
    nmf_model.fit(X)
    w = nmf_model.transform(x)
    H = nmf_model.components_

    finger_print_matrix = finger_print_proportion(x,w,H)
    index = [prefix + ' signature' + str(x) for x in range(1,k+1)]
    sig_coefficient = pd.DataFrame(H,index = index , columns=fit_df.index)
    index.append(prefix + ' residual')

    finger_print_df = pd.DataFrame(finger_print_matrix.T,index=index,columns=denovo_df.columns)
    finger_print_df[finger_print_df > 1] = 1
    finger_print_df[finger_print_df < 0] = 0

    return finger_print_df,sig_coefficient

def main() :
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",help="path of input metaphlan table")
    parser.add_argument("-k", "--number",type=int,help="number of signature to be decomposed")
    parser.add_argument("-p","--prefix",type=str,help="prefix of output signature")
    parser.add_argument("-o","--output",help="output path of probiotics signature")
    args = parser.parse_args()
    ab_matrix = pd.read_csv(args.input,sep = '\t',index_col = 0)
    #ab_matrix.columns = [x.split('_')[0] for x in ab_matrix.columns]
    
    if os.path.isdir(args.output) == False :
        os.mkdir(args.output)
    
    sig_ab_matrix,sig_coef_matrix = fit_signature(ab_matrix,ab_matrix,args.number,args.prefix)
    sig_ab_matrix.to_csv(args.output + args.prefix + "_" + "signature_proportion_matrix.txt",sep = '\t')
    sig_coef_matrix.to_csv(args.output + args.prefix + "_" + "signature_coefficient_matrix.txt",sep = '\t')
    
if __name__ == '__main__' :
    main()