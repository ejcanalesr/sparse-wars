% Erick J. Canales-Rodríguez, Jon Haitz Legarreta-Gorroño, Yasser Alemán-Gomez, Marco Pizzolato, Jonathan Rafael Patiño, Gaëtan Olivier D Rensonnet, Muhamed Barakovic, Joaquim Radua, Edith Pomarol-Clotet, Raymond Salvador, Jean-Philippe Thiran, Alessandro Daducci. Sparse wars: A survey and comparative study of spherical deconvolution algorithms for diffusion MRI. Neuroimage. 2019 Jan 1;184:140-160.

function [beta] = deconv_nnIRL1_L0_SD(Signal, Kernel)
% ---------------------
ns     = size(Kernel,2);
vzeros = zeros(ns,1);
beta   = vzeros;
% ---------------------
tau    = 1e-3;
Niter  = 30;

inv_W  = diag(beta + tau);  % W = diag(1./(beta + tau));
X      = Kernel*inv_W;
% ---------------------

lamb0 = 0.0015;
param.lambda     = lamb0;  % regularization parameter for l0-problem
param.pos        = true;
param.intercept  = false;
param.ista       = false;
param.loss       = 'square';
param.regul      = 'l1';

for j=1:Niter
    q      = mexFistaFlat(Signal, X, vzeros, param);
    q      = abs(q);
    beta   = inv_W*q;
    inv_W  = diag(beta + tau);
    X      = Kernel*inv_W;
end
end
