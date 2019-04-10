# Generative modeling and latent space arithmetics predict single-cell perturbation response across cell types, studies and species.

<img align="center"  src="/sketch/sketch.png?raw=true">



This repository includes python scripts and notebooks in the [code](https://github.com/theislab/scGen/tree/master/code) folder to reproduce figures from the paper [(bioRxiv, 2018)](https://www.biorxiv.org/content/10.1101/478503v2) according to the table bellow.



figure       | notebook/script     
---------------| ---------------
| 2, s1, s2, s3, s4, s7  | scgen_kang, cvae, st_gan, vec_arith_pca, vec_arith, scgen_sal_ta | 
|        3          | scgen_hpoly, scgen_salmonella| 
|        4          | cross_study| 
|        5, s8      | cross_species|
|        6, s9      | pancrease, bbknn, mnn, cca, scanorama|
|        s6      |scgen_kang_multiple|
|        s10        |mouse_atlas| 
| [*Figure 2*](https://nbviewer.jupyter.org/github/M0hammadL/scGen_reproducibility/blob/master/Jupyter%20Notebooks/Fig2.ipynb)| Fig2.ipynb| 

To run the notebooks and scripts you need following packages :

tensorflow, scanpy, numpy, matplotlib, scipy, wget.



