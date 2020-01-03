% Erick J. Canales-Rodríguez, Jon Haitz Legarreta-Gorroño, Yasser Alemán-Gomez, Marco Pizzolato, Jonathan Rafael Patiño, Gaëtan Olivier D Rensonnet, Muhamed Barakovic, Joaquim Radua, Edith Pomarol-Clotet, Raymond Salvador, Jean-Philippe Thiran, Alessandro Daducci. Sparse wars: A survey and comparative study of spherical deconvolution algorithms for diffusion MRI. Neuroimage. 2019 Jan 1;184:140-160.

function [fODF] = deconv_iter_nnSBL_fast_SD(Signal, Kernel, fODF0)
% ------------ Non-Negative SBL penalty using IRL2 approach ---------------
% -------------------------------------------------------------------------
% --- Regularization parameters
thresh  = sqrt(5e-4);  % -> equivalent to 5e-4 in fODF scale

% --- Precompute variables
KernelT = Kernel'; 

% --- Initialization
w = sqrt(fODF0);
N = size(Kernel,1);
I = eye(N);

lambda    = (1/20)^2; % initial value
Niter     = 30;
sig       = 1e-16;
W         = w.^4 + sig;
alpha     = 1e-6;
scalf     = 0.85;
% -------------------------------------

% --- Iterative estimation
for i = 1:Niter
    % --- find significant elements ( i.e., f > thresh^2 * max(f) )
    w(abs(w) <= thresh*max(abs(w))) = 0;
    ids = abs(w) > 0;
    natoms = sum(ids);
    
    if natoms < N
        % --- fast update on the selected elements (Non-Negative IRL2 scheme)
        % equivalent formulation which requires the inversion of a natomsxnatoms system

        Tinv = diag(1./W(ids));        
        Kc  = Kernel(:,ids);
        KTc = KernelT(ids,:);
        
        w(ids) = 1./w(ids) .* ( (KTc*Kc + lambda*Tinv)\KTc*Signal );
        
        % --- SBL-parameters update
        % Ref: "Solving Sparse Linear Inverse Problems: Analysis of Reweighted L1 and L2 Methods"
        %       David Wipf and Srikantan Nagarajan
        for kk = 1:Niter
            factor = W(ids) .* ( 1 - diag( ((KTc*Kc + alpha*diag(1./W(ids)))\KTc)*Kc ) );
            factor = abs(factor);
            W(ids) = w(ids).^4 + factor;
        end
    elseif natoms >= N
        % --- fast update on the selected elements (Non-Negative IRL2 scheme)
        % equivalent formulation which requires the inversion of a NxN system
        T   = diag(W(ids));
        Kc  = Kernel(:,ids);
        KTc = KernelT(ids,:);
        
        w(ids) = 1./w(ids) .* ( T*KTc * ( (Kc*T*KTc + lambda*I)\Signal) );
        
        % --- SBL-parameters update
        % Ref: "Solving Sparse Linear Inverse Problems: Analysis of Reweighted L1 and L2 Methods"
        %       David Wipf and Srikantan Nagarajan
        for kk = 1:Niter
            Bs     = ( alpha*I + Kc*diag(W(ids))*KTc );
            factor = W(ids) - W(ids).^2 .* diag( KTc * (Bs\Kc) );
            factor = abs(factor);
            W(ids) = w(ids).^4 + factor;
        end
    end
    % Ref: "Sparse Bayesian Learning for Basis Selection"
    %      IEEE TRANSACTIONS ON SIGNAL PROCESSING, VOL. 52, NO. 8, AUGUST 2004
    %      David P. Wipf and Bhaskar D. Rao

    s2 = sum((Signal - Kc*(w(ids).^2)).^2)/N;
    % correct the variance
    s2 = min((1/15)^2, max(s2,(1/30)^2)); % estimator in the interval sigma = [1/15, 1/30]
    % update parameter
    lambda = s2  + (lambda/N) * sum(1 - factor./W(ids));
    lambda = scalf * lambda; % better results using scalf < 1 (i.e., 0.8-0.9)
end
% calculate the final (non-negative) fODF
fODF = w.*w;
end
