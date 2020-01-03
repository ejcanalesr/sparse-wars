% Erick J. Canales-Rodríguez, Jon Haitz Legarreta-Gorroño, Yasser Alemán-Gomez, Marco Pizzolato, Jonathan Rafael Patiño, Gaëtan Olivier D Rensonnet, Muhamed Barakovic, Joaquim Radua, Edith Pomarol-Clotet, Raymond Salvador, Jean-Philippe Thiran, Alessandro Daducci. Sparse wars: A survey and comparative study of spherical deconvolution algorithms for diffusion MRI. Neuroimage. 2019 Jan 1;184:140-160.

function [fODF] = deconv_nnLASSO_BIC(Signal, Kernel, Vrec, neig, threshold)

[Ns, ~] = size(Kernel);
% --- BIC
g = 1;
Model_sel = @(RSS, n, k, nmodels) 2*k*g*log(nmodels) + k*log(n) + n*log(RSS/n);

lambda_max = norm( Kernel'*Signal, 'inf' );
factor = 1e-4;
% Niter  = 20;
mfODF  = SolveLasso(Kernel, Signal, size(Kernel,2), 'nnlasso', [], factor*lambda_max, [] ,1, []);
mfODF  = abs(mfODF);

k_fib = sum(mfODF > 0.05);
mfODF(:, k_fib > 6) = [];

Nlamb   = size(mfODF,2);
T_model = zeros(1,Nlamb);

SNR  = 20;
sig2 = (1/SNR)^2;

for i=1:Nlamb
    fodf_i =  mfODF(:,i);
    % --- extract peaks ---
    [~, ~, peaks] = find_discrete_peaks(Vrec, fodf_i(1:end-1), neig, threshold);
    fodf_i = fodf_i.*[peaks; 1];

    % --- debias step
    ind_p  = fodf_i > 0;
    Phi = Kernel(:, ind_p);
    fodf_i(ind_p) = pinv(Phi)*Signal;
    if min(fodf_i) >= 0
        % --- model selection
        k_compart  = sum(fodf_i > 0);
        RSS        = sum( ( Signal - sqrt( (Kernel*fodf_i).^2 + sig2 ) ).^2 );
        T_model(i) = Model_sel (RSS, Ns, k_compart, Nlamb);
    else
        T_model(i) = Inf;
    end
    mfODF(:,i) = fodf_i;
end

[~,   ind_opt]   = min(T_model);
fODF = mfODF(:,ind_opt(1));
end
