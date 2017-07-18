# Intlim:  Integration through Linear Modeling

## Introduction

Metabolomics is playing an increasing role in clinical and translational research and has greatly assisted in identifying novel, putative biomarkers for cancers.  Studies of metabolomics involve the measurements of small molecules (<1500 Daltons) in biospecimens such as blood, tissue, and urine, among others.  Given that metabolites reflect disease phenotype and downstream effects of biochemical pathways/post-translational modifications, they are seen as ideal candidates for biomarker discovery.  Several studies have utilized metabolomics to identify biomarkers for various cancers, heart disease, and diabetes.  

However, despite this progress, there are many challenges to interpreting metabolomics data- especially for untargeted studies that involve many unidentified metabolites (as opposed to targeted studies involving known metabolites)[6].  One of the challenges is to interpret the metabolomics profiles and understand how they affect, and are affected by genes and the proteins they produce.  Approaches have been developed to integrate transcriptomic and metabolomics data.  Pathway-based approaches such as Metaboanalyst, INMEX, and IMPALA integrate transcriptomic and metabolomics data for pathway enrichment analysis.  One caveat is that these approaches rely on previously curated pathways and “expert” definitions of what constitutes a pathway[10].  Correlation based approaches have been developed to integrate gene expression and metabolomics data.  MixOmics and DiffCorr are examples of such tools available.  While identifying globally correlated gene-metabolite pairs will enhance our understanding of how genes affect metabolic phenotypes, these co-regulated gene-metabolite pairs may not represent biological interactions specific to a phenotype of interest (e.g. diagnosis, hormone receptor status, cancer type, etc.).  To identify these phenotype-specific gene-metabolite relationships, we propose a novel linear modeling approach. 
The novel linear model is:  E(m|g,t) = β1 + β2 g + β3 t + β4 (g:t) + ε where m and g are log-transformed metabolite abundances and gene levels respectively, t is phenotype (cancer type, patient diagnosis, treatment group, etc), (g:t) is the interaction between gene expression and phenotype, and ε is the error term.  A statistically significant p value of the (g:t) interaction term indicates that the slope relating gene expression and metabolite abundance is different from one phenotype compared to another.  Through this model, we can identify gene-metabolite correlations that are specific to a particular phenotype.  

This model has been applied in this study to the publically available NCI-60 cancer cell line data on gene expression and metabolomics (available in this package) as well as previously published gene expression and metabolomics data from a breast cancer study.  To increase usability of our approach, we have also implemented our approach as an R Bioconductor package, IntLim or Integration through Linear Modeling, that includes an RShiny web application (which does not require knowledge of R).  

## Contact

If you have any questions, comments, or concerns on how to use IntLim please contact Ewy Mathe at ewy.mathe@osumc.edu.  You can also contact Jalal Siddiqui at jalal.siddiqui@osumc.edu