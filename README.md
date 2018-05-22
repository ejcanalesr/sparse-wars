# sparse-wars

**Implementation of the spherical deconvolution (Sd) algorithms described here:**

> **Sparse Wars: A Survey and Comparative Study of Spherical Deconvolution Algorithms for Diffusion MRI.**
Erick Jorge Canales-Rodríguez, Jon Haitz Legarreta, Marco Pizzolato, Gaëtan Rensonnet, Gabriel Girard, Jonathan Rafael Patiño, Muhamed Barakovic, David Romascano, Yasser Alemán-Gomez, Joaquim Radua, Edith Pomarol-Clotet, Raymond Salvador, Jean-Philippe Thiran, Alessandro Daducci (*Submitted to Neuroimage*, 2018).

**Code implementation:**
The current implementation is written in Matlab.

**List of reconstruction methods:**
1. `Best-subset selection based on the extended Bayesian information criterion (EBIC)`
2. `LASSO based on the EBIC`
3. `Non-negative Cauchy deconvolution`
4. `Non-negative regularized FOCal Underdetermined System Solver (FOCUSS)`
5. `Non-negative iterative reweighted l1 minimization`
6. `Sparse Bayesian Learning (SBL)`
7. `Accelerated damped Richardson-Lucy deconvolution (AdRL-SD)`
8. `Accelerated robust and unbiased model-based spherical deconvolution (ARUMBA-SD)`


**Install dependencies:**
- SparseLab Toolbox (https://sparselab.stanford.edu/)
- SPArse Modeling Software (SPAMS: http://spams-devel.gforge.inria.fr/) 
