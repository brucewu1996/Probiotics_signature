class probiotics_signature :

    def __init__(self,metadata,abundance_matrix,distance_matrix) :
        """
        Init a probiotics signature object.

        Args:
            metadata (pd.Dataframe): clinical information or demographic data of patients / participants. row is sample, column is condition.
            abundance_matrix (pd.Dataframe): patient abundance matrix. row is species, column is sample.
            distance_matrix (pd.Dataframe): n_sample * n_sample distance matrix. The value inside the matrix is the distance between samples. 
            (ex : Unifrac / Bray-curtis distance)
        """        
        self.sample_metadata = None
        self.abundance_matrix = None
        self.distance_matrix = None

    def species_abundance_prevelence_scatterplot(self,output_path,format='pdf') :
        """
        Scatterplot of target species abundance / prevalence
        """        
        if not self.abundance_matrix :
            print("Please input abundance matrix")
            return 

    def concat_species_into_subtype(self,subtype_dict) :
        """ 
        Merge abundance from rate species into subtype. (ex : Lactobacillus genus)

        Args:
            subtype_dict (dict): dict of subtype : species. key is species name, value is subtype name.
        """     
        if not self.abundance_matrix :
            print("Please input abundance matrix")
            return    

    def evaluate_nmf_component(self,k_min=1,k_max=10) :
        """
        Evaluate optimal k component for NMF decomposition processing.

        Args:
            k_min (int, optional): The minimum component number. Defaults to 1.
            k_max (int, optional): The maximum component number. Defaults to 10.
        """        
        pass

    def decompose_probiotics_signature(self,k) :
        """_summary_

        Args:
            k (int): Number of component expected to decompose.
        """        
        pass

    def plot_signature_relative_importance(self) :
        pass

    def plot_signature_relative_importance_heatmap(self) :
        pass

    def plot_pcoa_permanova_scatterplot(self) :
        pass

    def plot_mds_permanova_scatterplot(self) :
        pass

    def diversity_comparison_between_cluster(self,ordinate_method = 'PCoA') :
        pass