% Erick J. Canales-Rodríguez, Jon Haitz Legarreta-Gorroño, Yasser Alemán-Gomez, Marco Pizzolato, Jonathan Rafael Patiño, Gaëtan Olivier D Rensonnet, Muhamed Barakovic, Joaquim Radua, Edith Pomarol-Clotet, Raymond Salvador, Jean-Philippe Thiran, Alessandro Daducci. Sparse wars: A survey and comparative study of spherical deconvolution algorithms for diffusion MRI. Neuroimage. 2019 Jan 1;184:140-160.

clear; clc;

% ========================== Path to the data =========================== %

path_to_folder     = '/home/erick/Code/preprocesed_data_hcp';
path_to_data       = [path_to_folder '/dwihcpl_den_gibss_eddy.nii'];
b_vector           = [path_to_folder '/gradientshcpl.txt'];
b_value            = [path_to_folder '/bvalueshcpl.txt'];
path_to_mask       = [path_to_folder '/mask.nii'];

% ----------------- Format to save the fODFs and maxima ----------------- %
format_to_save_fODF = 'ODF-FSL-dyads';

path_to_save = '/home/erick/SBL_recon';
mkdir(path_to_save);

% ======= Kernel/Dictionary used for the model-based deconvolution ====== %
% Define the diffusivities in your data (i.e, DTI eigenvalues for regions of parallel fibers), for example:
lambda1    = 1.7e-3 % mm^2/s  (parallel to the fiber)
lambda2    = 0.23-3 % mm^2/s  (perpendicular to the fiber)
lambda_csf = 3.0e-3 % mm^2/s  (free diffusion coefficient)

tic
display('=================================================');
SBL_SD(path_to_data, path_to_mask, b_vector, b_value, format_to_save_fODF, path_to_save, lambda1, lambda2, lambda_csf);
display('=================================================');
total_time1 = toc;
disp(datestr(datenum(0,0,0,0,0,total_time1),'HH:MM:SS'));
