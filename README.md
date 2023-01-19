# Probiotics signature

![Untitled](Github%20repo%20readme%20f32432f57cf7401e822b5357555a072a/Untitled.png)

---

## Introduction

Probiotics signature provide *specific probiotics combination* via non-negative matrix factorization (NMF) based on relative abundance matrix from population. [In previous probiotics intervention research](https://www.aging-us.com/article/102810/text), probiotics for intervention were designed by bacterial function rather than abundance of bacteria in participants. We can categorize people into multiple sub-class according to their probiotics signature profile to achieve precise supplementation of probiotics.

---

## Requirement 
1. Taxonomy classification by [MetaPhlAn2](https://github.com/biobakery/MetaPhlAn).
2. Extract probiotics signature via [NMF](https://cran.r-project.org/web/packages/NMF/index.html) & *Scikit-learn* decomposition module.
3. Population grouping through **[ConsensusClusterPlus](https://bioconductor.org/packages/release/bioc/html/ConsensusClusterPlus.html)**.
4. Overall taxonomic similarity by *weighted-Unifrac distance* & *PERMANOVA.*
