% Erick J. Canales-Rodríguez, Jon Haitz Legarreta-Gorroño, Yasser Alemán-Gomez, Marco Pizzolato, Jonathan Rafael Patiño, Gaëtan Olivier D Rensonnet, Muhamed Barakovic, Joaquim Radua, Edith Pomarol-Clotet, Raymond Salvador, Jean-Philippe Thiran, Alessandro Daducci. Sparse wars: A survey and comparative study of spherical deconvolution algorithms for diffusion MRI. Neuroimage. 2019 Jan 1;184:140-160.

function [fODF] = deconv_best_subset(Kernel, Signal, fODF0)

[Ns, Nc] = size(Kernel);

ind_non_zero_initial = find(fODF0>0);
Nmax = length(ind_non_zero_initial);
if Nmax > 20
    % keep only the higher 20 values to avoid memory errors...
    Nmax = 20;
    [~, ind_non_zero_initial] = sort(fODF0, 'descend');
    ind_non_zero_initial      = ind_non_zero_initial(1:Nmax);
end

Kerneln = Kernel(:,ind_non_zero_initial);

% Log-Likelihood term
% ---------------- Gausssian-Rician Likelihood model ---------------- %
LogLike = @(Y, S, sig2) Ns*log((sum(( Y - sqrt(S.^2 + sig2) ).^2))/Ns);

% BIC equation
g    = 1;
BIC  = @(n, k, LogLike, Nvar) 2*k*g*log(Nvar) + k*log(n) + LogLike;

SNR    = 20; % expected SNR
sigma2 = (1/SNR)^2;

% --- All possible simple combinations of 1, 2 and 3 fibers (+ isotropic)
index = dec2bin(1:2^Nmax-1);
index = index == '1';
n_fibers   = sum(index,2);
max_fibers = 3;
if fODF0(end) > 0
    max_fibers = 4; % 3 fibers + isotropic compartment
end
index(n_fibers > max_fibers,:) = []; % deleting complex configurations!

% ---- testing the combinations via EBIC
Ncomb = size(index,1);
T_BIC = zeros(Ncomb,1);
for ni = 1:Ncomb
    Phi  = Kerneln(:,index(ni,:));
    fODF = pinv(Phi)*Signal;
    
    % -- BIC evaluation
    if min(fODF) >= 0 
        LogLike_i = LogLike(Signal,(Phi*fODF), sigma2);
        k_fib = sum(fODF>0);
        T_BIC(ni) = BIC (Ns, k_fib, LogLike_i, Ncomb);
    else
        T_BIC(ni)  = Inf;
    end
end

% estimation on the optimal subset
[~, ind_min] = min(T_BIC);
ind_min_relative = index(ind_min,:);
Phi = Kerneln(:,ind_min_relative);

fODF = zeros(Nc,1);
fODF(ind_non_zero_initial(ind_min_relative)) = pinv(Phi)*(Signal);
end
