import pandas as pd
import numpy as np
import re
from numba import jit
from scipy.stats import spearmanr,pearsonr
import networkx as nx
import random
import time


@jit(nopython=True)
def enrich_score(taxa_array) : 
    '''
    pathway : name of pathway, ex : DTDPRHAMSYN-PWY: dTDP-&beta;-L-rhamnose biosynthesis
    taxonomy_list : taxonomy list with "order"
    '''
    #calculate enrichment score unit, enrichment score upper bond = 1, lower bond = -1
    es = 0
    max_es = 0
    es_pos = 1/sum(taxa_array)
    es_n = 1/(len(taxa_array) - sum(taxa_array))
    #calculate enrichment score by "ordered" taxonomy list
    for t in taxa_array :
        if t  == 1 :
            es += es_pos
        elif t == 0 :
            es -= es_n
        else :
            print('Value in taxa array must be 0 or 1')
            break
        if es > max_es :
            max_es = es
    return max_es

@jit(nopython=True)
def permutation(taxa_array,n_times = 1000) :
    
    origin_es = enrich_score(taxa_array)
    permutation_es = np.zeros(n_times)
    for i in range(n_times) :
        np.random.shuffle(taxa_array)
        permutation_es[i] = enrich_score(taxa_array)
    es_above_origin = sum(permutation_es > origin_es)
    if es_above_origin == 0 :
        pesudo_f = 1/ (n_times * 10)
    else :
        pesudo_f = es_above_origin / n_times
    return pesudo_f,origin_es

class taxonomy_enrich_analysis :

    def __init__(self,deseq_output,taxonomy,metaphlan_output,humann_index,graph) :
        '''
        self.deseq : dataframe, output of deseq2
        self.metaphlan_output : dataframe, output of metaphlan
        slef.taxonomy : list, list of overall taxonomy
        self.humann_index : index,
        '''
        self.deseq = deseq_output
        self.metaphlan_output = metaphlan_output
        self.taxonomy = taxonomy
        ### source from networkx or humann
        self.graph = graph
        self.humann_index = humann_index
        if self.humann_index != None and self.graph == None:
            p_idx = [bool(re.search("^(?!.*\|)",x)) for x in self.humann_index]
            self.pathway_list = self.humann_index[p_idx]
        elif self.graph != None and self.humann_index == None:
            nx.get_node_attributes(self.graph,'Cluster')
            cluster_list = list(set(nx.get_node_attributes(self.graph,'Cluster').values()))
            cluster_list.remove('No cluster')
            self.pathway_list = cluster_list
            
        self.pseudo_p_value = np.zeros(len(self.pathway_list))
        self.es_score = np.zeros(len(self.pathway_list))

    def sort_taxonomy(self,method='fold-change') :
        '''
        for TSEA node mode, input fold-change / p-value from deseq2 output
        '''
        if method == 'fold-change' :
            self.deseq = self.deseq.sort_values(by='log2FoldChange',ascending=False)
            self.taxonomy = list(self.deseq.index)
        elif method == 'p-value' :
            self.deseq = self.deseq.sort_values(by='padj',ascending=False)
            self.taxonomy = list(self.deseq.index)
        else :
            self.taxonomy = self.taxonomy

    def sort_taxonomy_edge(self,method = 'spearman') :
        '''
        For TSEA edge mode, this function will compute correlation between features.
        Specially, according to characteristic of metagenomic data, process of computing correlation for as least one feature is non-zero.
        '''
        corr_z_dict = dict()
        corr_df = self.metaphlan_output
        n = corr_df.shape[0]
        for col in range(n-1) :
            for row in range(col+1,n) :
                name = (corr_df.index[col],corr_df.index[row])
                t1 = corr_df.iloc[col,:]
                t2 = corr_df.iloc[row,:]
                idx = (t1 + t2) > 0
                if method == 'spearman' :
                    r,p = spearmanr(t1[idx],t2[idx])
                elif method == 'pearson' :
                    r,p = pearsonr(t1[idx],t2[idx])
                else :
                    print('Only prove pearson / spearman correlation method')
                z = np.arctanh(r)
                corr_z_dict[name] = z
        df = pd.DataFrame({'Name' : corr_z_dict.keys(),'z-score' : corr_z_dict.values()}).sort_values(by='z-score',ascending=False)
        self.taxonomy = list(df['Name'].values)
        del df
    
    def format_taxa_array(self,pathway,taxa_list,mode='Node') :
        if self.graph == None and self.humann_index == None :
            taxa_idx = [bool(re.search(pathway,x)) for x in self.humann_index]
            taxa_in_path = self.humann_index[taxa_idx]
            taxa_in_pathway = lambda x : x.split('|')[-1].split('.')[-1]
            taxa_in_path = list(map(taxa_in_pathway,taxa_in_path[1:]))
            taxa_array = np.zeros(len(taxa_list))

        elif self.graph != None and self.humann_index == None :
            taxa_in_path = [x for x,y in self.graph.nodes(data=True) if y['Cluster'] == pathway]
            taxa_array = np.zeros(len(taxa_list))

        for idx,t in enumerate(taxa_list) :
            if mode == 'Node' :
                if t in taxa_in_path :
                    taxa_array[idx] = 1
            else :
                if t[0] in taxa_in_path and t[1] in taxa_in_path :
                    taxa_array[idx] = 1

        return taxa_array

    def TSEA(self,mode='Node',source='humann',sorted=False) :
        if sorted == False :
            if mode == 'Node' :
                self.sort_taxonomy()
            elif mode == 'Edge' :
                self.sort_taxonomy_edge()
            else :
                print('Only have Node/Edge mode for taxonomy set enrichment analysis')
                return
        print('taxonomy list length : %d' % len(self.taxonomy))
        for idx,p in enumerate(self.pathway_list) :
            start_time = time.time()
            print('Calculate enrichment score of {idx} / {total}  pathway : {pathway}'.format(idx =idx+1,total = len(self.pathway_list),pathway =p))
     
            if mode == 'Node' :
                taxa_array = self.format_taxa_array(p,self.taxonomy)
            else :
                taxa_array = self.format_taxa_array(p,self.taxonomy,mode='Edge')

            if sum(taxa_array) == 0 :
                pass
            else :
                self.pseudo_p_value[idx],self.es_score[idx] = permutation(taxa_array)
            
            end_time = time.time()
            time_delta = round(end_time - start_time,2)
            #print('Pathway : {p} pesudo-F : {f}, enrich score : {es}'.format(p=p,f=self.pseudo_p_value[idx],es=self.es_score[idx] ))
            print('TSEA of pathway : %s cost %0.2f s' % (p,time_delta))

        tsea_result = pd.DataFrame({'Pathway' : self.pathway_list,'Enrich_score' : self.es_score,'p-value' : self.pseudo_p_value})
        return tsea_result

def main() :

    deseq_df = pd.read_csv('/home/bruce1996/ssd/mci_deseq_res.txt',sep='\t',index_col=0).sort_values(by='log2FoldChange',ascending=False)
    metaphlan_output = "/home/bruce1996/nvme2/mci_for_adlasso/data/MCI_species_relative_abundance.txt"
    meta_df = pd.read_csv(metaphlan_output,sep = '\t',index_col=0)
    humann_df = pd.read_csv('/home/bruce1996/data/MCI/merge_mci_huamnn_result.txt',sep='\t',index_col=0)
    ##remove un of humann result
    idx = [bool(re.search('UN',x)) == False for x in humann_df.index]
    humann_df = humann_df.loc[idx,:]
    ##remove special charcter from humann output index
    removeSpecialChars = lambda x : x.translate ({ord(c): "" for c in "[]"})
    humann_df.index = list(map(removeSpecialChars,list(humann_df.index) ))
    tsea_output =  '/home/bruce1996/data/MCI/tsea/'

    tsea = taxonomy_enrich_analysis(deseq_output=deseq_df,metaphlan_output=meta_df,humann_index=humann_df.index,graph=None,taxonomy=None)
    node_df = tsea.TSEA(mode='Node',sorted=False).sort_values(by='Pathway')
    edge_df = tsea.TSEA(mode='Edge',sorted=False).sort_values(by='Pathway')
    tsea_df = pd.DataFrame({'Name' : node_df['Pathway'],'Node-pseudo-F' : node_df['p-value'],'Edge-pseudo-F' : edge_df['p-value']})
    tsea_df.to_csv(tsea_output + 'Humann-tsea-result.txt',sep='\t')

if __name__ == '__main__' :
    main()
