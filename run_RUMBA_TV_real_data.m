% Erick J Canales-Rodríguez, Alessandro Daducci, Stamatios N Sotiropoulos, Emmanuel Caruyer, Santiago Aja-Fernández, Joaquim Radua, Yasser Iturria-Medina, Lester Melie-García, Yasser Alemán-Gómez, Jean-Philippe Thiran, Salvador Sarró, Edith Pomarol-Clotet, Raymond Salvador. Spherical Deconvolution of Multichannel Diffusion MRI Data with Non-Gaussian Noise Models and Spatial Regularization. PLoS ONE, 2015, 10(10): e0138910.

clear all; clc;

% ========================== Path to the data =========================== %
path_to_folder  = '/home/erick/Code/preprocesed_data_hcp';
path_to_data    = [path_to_folder '/dwihcpl_den_gibss_eddy.nii'];
b_vector        = [path_to_folder '/gradientshcpl.txt'];
b_value         = [path_to_folder '/bvalueshcpl.txt'];
path_to_mask    = [path_to_folder '/mask.nii'];

% ======= Kernel/Dictionary used for the model-based deconvolution ====== %
% Define the diffusivities in your data (i.e, DTI eigenvalues for regions of parallel fibers), for example:
lambda1 = 1.7*10^-3;% mm^2/s  (parallel to the fiber)  
lambda2 = 0.2*10^-3;% mm^2/s  (perpendicular to the fiber) 
% --- mean diffusivity assumed in CSF
lambda_csf = 3.0*10^-3;% mm^2/s
% --- mean diffusivity assumed in GM
lambda_gm  = 0.8*10^-4;% mm^2/s

% ===================== Parameters of the algorithm ===================== %
Niter = 600; % number of iterations for convergence 
              % (600 is my recommended value)

% --- Number of coils in your scanner
ncoils = 32;
%ncoils = 8;

% --- MRI reconstruction method 
% 'SMF-SENSE-based'  is valid for SENSE, mSENSE, ASSET and related approaches
% 'SoS-GRAPPA-based' is valid for GRAPPA, SMASH, AUTO-SMASH, PILS and related approaches
MRI_recon_type = 'SMF-SENSE-based';
%MRI_recon_type = 'SoS-GRAPPA-based';

% --- Acceleration factor (R) of the acquisition (i.e., level of the k-space subsampling)
% For SIEMENS, R = iPAT factor. For GE, R = ASSET factor. For PHILIPS, R = SENSE factor
% Typical values used in practice are R = 2 or R = 1
% When the image is acquired using full k-space coverage (without using fast pMRI options) R = 1
% R = 2;
R = 1;

% ----------------- Format to save the fODFs and maxima ----------------- %
format_to_save_fODF = 'ODF-FSL-dyads';

% ---------------- use Total Variation regularization ------------------- %
use_tv = 'yes';
%use_tv = 'no';

recon_folder = ['/home/erick/RUMBA_SD_recon_all_preproc'];
path_to_save = [recon_folder '/RUMBA_SD_TV' '_Recon_' num2str(Niter)];
mkdir(path_to_save);

tic
disp('=================================================');
if strcmp( use_tv, 'yes')
    disp ('The reconstruction method is RUMBA-SD+TV')
else
    disp ('The reconstruction method is RUMBA-SD');
end
RUMBA_sph_deconv(path_to_data, path_to_mask, b_vector, b_value, format_to_save_fODF, path_to_save, Niter, lambda1, lambda2, lambda_csf, lambda_gm, ncoils, MRI_recon_type, R, use_tv);

disp('=================================================');
total_time1 = toc;
disp(datestr(datenum(0,0,0,0,0,total_time1),'HH:MM:SS'));
