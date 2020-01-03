function [fODF_4D, fgm_3D, fcsf_3D, SNR_est, std_SNR, var_3D] = sph_deconv_tv_motor_all(Signal, Kernel, fODF0, Niter, Mask, ncoils, MRI_recon_type, R, use_tv)

% Erick J Canales-Rodríguez, Alessandro Daducci, Stamatios N Sotiropoulos, Emmanuel Caruyer, Santiago Aja-Fernández, Joaquim Radua, Yasser Iturria-Medina, Lester Melie-García, Yasser Alemán-Gómez, Jean-Philippe Thiran, Salvador Sarró, Edith Pomarol-Clotet, Raymond Salvador. Spherical Deconvolution of Multichannel Diffusion MRI Data with Non-Gaussian Noise Models and Spatial Regularization. PLoS ONE, 2015, 10(10): e0138910.

% ----------------------------------------------------
% fODF_4D: purely anisotropic part
% fgm_3D : fraction of the isotropic component for GM
% fcsf_3D: fraction of the isotropic component for CSF
% SNR_est: Estimated mean SNR

if ( strcmp( MRI_recon_type, 'SMF-SENSE-based') ) || (ncoils == 1)
    n_order = 1;
elseif strcmp( MRI_recon_type, 'SoS-GRAPPA-based')
    n_order = ncoils;
end

% ------- Mask
Mask1D = Mask(:);
ind_mask  = find(Mask1D == 1);
% -------

Dim  = size(Signal);
Data_2d = zeros(prod(Dim(1:3)),Dim(4),'single');

for i=1:Dim(4)
    Data_3d = squeeze(Signal(:,:,:,i));
    Data_2d(:,i) = Data_3d(:);
end
clear Data_3d Signal;
Data_2d = Data_2d';

Data_2d = squeeze(Data_2d(:,ind_mask));

fODF = repmat(fODF0, [1, size(Data_2d, 2)]);

Reblurred = Kernel*fODF;

fzero = 0;
KernelT = Kernel'; % only one time

lambda_aux = zeros(prod(Dim(1:3)),1,'single');

%%________________________ Main Algorithm ______________________________ %% 
sigma0 = 1/15;
sigma2 = sigma0^2;
lambda = sigma2;
epsilon = 1e-7;
N = size(Data_2d, 1);
sigma2 = sigma2*ones(size(Data_2d));
% --------------
Reblurred_S = (Data_2d.*Reblurred)./sigma2;
for i = 1:Niter
    display(['Iter -> ' num2str(i) ' of ' num2str(Niter)]);
    % ------------------- R-L deconvolution part -------------------- %
    
    fODFi = fODF;
    
    Ratio = mBessel_ratio(n_order,Reblurred_S);
    
    RL_factor = KernelT*( Data_2d.*( Ratio ) )./( KernelT*(Reblurred) + eps);
        
    % -------------------- TV regularization part ------------------- %
    
    if strcmp( use_tv, 'yes')
        
        TV_factor = ones(size(fODFi),'single');
        
        for j = 1:size(fODFi,1)
            fODF_jv = zeros(prod(Dim(1:3)),1,'single');
            
            fODF_j = squeeze( fODFi(j,:) );
            fODF_jv(ind_mask) = fODF_j;
            
            fODF_3D_j = reshape(fODF_jv,[Dim(1) Dim(2) Dim(3)]);
            
            % --------------------------- %
            clear fODF_j fODF_jv
            
            Gr = grad(fODF_3D_j);
            d = sqrt(sum(Gr.^2,4));
            d = sqrt( epsilon^2 + d.^2 );
            div_f = div( Gr ./ repmat( d, [1 1 1 3]) );
            G0 = abs( 1 - lambda.*div_f ); % abs -> to avoid negative values
            TV_factor3D =  1./(G0 + eps); % eps -> to avoid zero values
            
            TV_factor1D = TV_factor3D(:)';
            
            TV_factor1D = squeeze(TV_factor1D(ind_mask));
            
            TV_factor(j,:) = TV_factor1D;
            % ---------------
            
            clear Gr d div_f G0 TV_factor3D TV_factor1D
        end
        % ------------------------------- %
        RL_factor = RL_factor.*TV_factor;
    end
    % --------------------------------------------------------------- %
    
    fODF = fODFi.*RL_factor;
    
    fODF = max(fzero, fODF); % positivity
    
    %if i <= 100
       %cte = sqrt(1 + 2*n_order*mean(sigma2(:)));
       %cte = sqrt(1 + n_order*mean(sigma2(:)));
 %      cte = 1;
 %      fODF =  cte*fODF./repmat( sum(fODF,1) + eps, [size(fODF,1), 1] ); 
       % Energy preservation at each step, which included the bias on the S0 image.
       % This step is used only during the first iterations to stabilize
       % the recovery at the begining...
    % end
        
    % --------------------- Update of variables --------------------- %
    Reblurred = Kernel*fODF;
    Reblurred_S = (Data_2d.*Reblurred)./sigma2;
    
    % ---------------- noise ---------------- %
    sigma2_i = (1/N)*sum( (Data_2d.^2 + Reblurred.^2)/2 - (sigma2.*Reblurred_S).*Ratio, 1)./n_order;
    sigma2_i = min((1/8)^2, max(sigma2_i,(1/80)^2)); % estimator in the interval sigma = [1/SNR_min, 1/SNR_max],
                                                      % where SNR_min = 8 and SNR_max = 80
    % ---------- automatic lambda ----------- %
    % (:,ind_mask)
    mean_SNR = mean( 1./sqrt( sigma2_i ) );
    std_SNR  =  std( 1./sqrt( sigma2_i ) );
    display(['Estimated Mean SNR (S0/sigma) -> ' num2str( mean_SNR ) ' (+-) ' num2str( std_SNR )]);
    display(['Reconstruction using -> ' num2str( n_order ) ' coils ']);
    
    meanfODF = mean(sum(fODF,1));
    display(['mean sum(fODF) -> ' num2str( meanfODF )]);

    % --------------------------------------- %
    sigma2 = repmat(sigma2_i,[size(Data_2d,1), 1]);
    
    if strcmp( use_tv, 'yes')
        % Discrepancy principle
        if R==1
            lambda = mean(sigma2_i); % using the same penalization for all the voxels
                                     % assuming equal variance across the image
            if mean(sigma2_i) < (1/30)^2
                lambda = (1/30)^2;   % force TV for low level of noise
            end
        elseif R > 1
            lambda_aux(ind_mask) = sigma2_i;
            lambda = reshape(lambda_aux, [Dim(1) Dim(2) Dim(3)]); % adaptive spatial regularization for fast pMRI
                                                                  % assuming spatial dependence for the variance
                                                                  % e.g., due to pMRI or different values for different tissue types
            % (in the near future, this parameter could be low-pass filtered for robust
            % estimation...)
        end
    end
end

fODF = fODF./repmat( sum(fODF,1) + eps, [size(fODF,1), 1] ); % energy preservation

SNR_est = mean_SNR;
var_3D = zeros(prod(Dim(1:3)),1,'single');
var_3D(ind_mask) = sigma2_i;
var_3D = reshape(var_3D, [Dim(1) Dim(2) Dim(3)]);
clear lambda_aux lambda
%%____________________ End of Main Algorithm ___________________________ %% 

% Reshaping the ODF from 2D to 4D
nc = length(fODF0);
fODF_4D = zeros(Dim(1), Dim(2), Dim(3), nc-2, 'single');
for i=1:nc-2
    fODF_1D_i = zeros(prod(Dim(1:3)),1,'single');
    
    fODF_1D = squeeze(fODF(i,:));
    
    fODF_1D_i(ind_mask) = fODF_1D;
    
    fODF_4D(:,:,:,i) = reshape(fODF_1D_i, [Dim(1) Dim(2) Dim(3)]); % purely anisotropic part
end
%
fODF_1D = zeros(prod(Dim(1:3)),1,'single');

fODF_jv = squeeze(fODF(nc,:));

fODF_1D(ind_mask) = fODF_jv;

fgm_3D  = reshape(fODF_1D,[Dim(1) Dim(2) Dim(3)]); % fraction of the isotropic component for GM
% -----
fODF_1D = zeros(prod(Dim(1:3)),1,'single');

fODF_jv = squeeze(fODF(nc-1,:));

fODF_1D(ind_mask) = fODF_jv;

fcsf_3D = reshape(fODF_1D,[Dim(1) Dim(2) Dim(3)]); % fraction of the isotropic component for CSF

end

%% ________________________ Private Functions ___________________________%%

function [fxyz] = grad(M)
% grad - gradient operator
fx = M([2:end end],:,:)-M;
fy = M(:,[2:end end],:)-M;
fz = M(:,:,[2:end end])-M;

fxyz = cat(4,fx,fy,fz);
end

function fd = div(P)
% div - divergence operator

Px = squeeze(P(:,:,:,1));
Py = squeeze(P(:,:,:,2));
Pz = squeeze(P(:,:,:,3));

fx = Px-Px([1 1:end-1],:,:);
fx(1,:,:)   = Px(1,:,:);    % boundary
fx(end,:,:) = -Px(end-1,:,:);
% ---
fy = Py-Py(:,[1 1:end-1],:);
fy(:,1,:)   = Py(:,1,:);    % boundary
fy(:,end,:) = -Py(:,end-1,:);
% ---
fz = Pz-Pz(:,:,[1 1:end-1]);
fz(:,:,1)   = Pz(:,:,1);    % boundary
fz(:,:,end) = -Pz(:,:,end-1);
% ---
fd = fx+fy+fz;
end

function y = mBessel_ratio(n,x)
% y = mBessel{n}(x)/mBessel{n-1}(x) = besseli(n,x)./besseli(n-1,x)
% Fast evaluation using the Perron's Continued Fraction equation.
% For more details see: http://www.ams.org/journals/mcom/1978-32-143/S0025-5718-1978-0470267-9/

y = x./( (2*n + x) - ( 2*x.*(n+1/2)./ ( 2*n + 1 + 2*x - ( 2*x.*(n+3/2)./ ( 2*n + 2 + 2*x - ( 2*x.*(n+5/2)./ ( 2*n + 3 + 2*x ) ) ) ) ) ) );
end
