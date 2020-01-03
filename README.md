# Deconvolution algorithms for diffusion MRI

**Implementation of the algorithms described here:**

> **Sparse Wars: A Survey and Comparative Study of Spherical Deconvolution Algorithms for Diffusion MRI.**
Erick J Canales-Rodríguez, Jon Haitz Legarreta, Marco Pizzolato, Gaëtan Rensonnet, Gabriel Girard, Jonathan Rafael Patiño, Muhamed Barakovic, David Romascano, Yasser Alemán-Gomez, Joaquim Radua, Edith Pomarol-Clotet, Raymond Salvador, Jean-Philippe Thiran, Alessandro Daducci. Neuroimage, 2019 (https://www.sciencedirect.com/science/article/abs/pii/S1053811918307699?via%3Dihub)

> **Spherical Deconvolution of Multichannel Diffusion MRI Data with Non-Gaussian Noise Models and Spatial Regularization.**
Erick J Canales-Rodríguez, Alessandro Daducci, Stamatios N Sotiropoulos, Emmanuel Caruyer, Santiago Aja-Fernández, Joaquim Radua, Yasser Iturria-Medina, Lester Melie-García, Yasser Alemán-Gómez, Jean-Philippe Thiran, Salvador Sarró, Edith Pomarol-Clotet, Raymond Salvador. PLoS ONE, 2015, 10(10): e0138910. (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0138910)


The current implementation is written in Matlab.

**List of included methods:**
1. `Best-subset selection based on the extended Bayesian information criterion (NNLS-BSS-EBIC)`
2. `LASSO based on the EBIC (LASSO-EBIC)`
3. `Non-negative iterative reweighted l1 minimization (IRL1)`
4. `Sparse Bayesian Learning (SBL)`
5. `Robust and unbiased model-based spherical deconvolution (RUMBA-TV)`


**Install dependencies:**
- SparseLab Toolbox (https://sparselab.stanford.edu/)
- SPArse Modeling Software (SPAMS: http://spams-devel.gforge.inria.fr/)
- SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
