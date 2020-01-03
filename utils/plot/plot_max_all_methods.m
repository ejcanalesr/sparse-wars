% Erick J. Canales-Rodríguez, Jon Haitz Legarreta-Gorroño, Yasser Alemán-Gomez, Marco Pizzolato, Jonathan Rafael Patiño, Gaëtan Olivier D Rensonnet, Muhamed Barakovic, Joaquim Radua, Edith Pomarol-Clotet, Raymond Salvador, Jean-Philippe Thiran, Alessandro Daducci. Sparse wars: A survey and comparative study of spherical deconvolution algorithms for diffusion MRI. Neuroimage. 2019 Jan 1;184:140-160.

% Script to plot all the ODFs

clear; clc; close all

% For example ...

% path_to_mask = ['brain_mask.nii'];
path_to_mask = ['/home/erick/Documents/Sparse_wars/Real_data/dti_fa.nii'];

invert_z = 0;

type = 'y';
if strcmp( type, 'y')
    %slice = 73;
    %slice = 74;
    %slices = [69, 70, 74];
    %x1 = 44; x2 = 83;
    %y1 = 27; y2 = 45;
    slice = 69;
    
    x1 = 46; x2 = 62;
    y1 = 27; y2 = 44;
    plot_axis_limit = [x1, x2, slice-2, slice+2, y1, y2];
end

plot_SBL      = 0;
plot_RFOCUSS  = 0;
plot_LASSO    = 0;
plot_ARUMBA   = 0;
plot_IRL1     = 0;
plot_NNLS_EBI = 0;
plot_dRL_SD   = 0;
plot_Cauchy   = 1;


if plot_SBL == 1
    path_to_peaks = ['/home/erick/Documents/Sparse_wars/Real_data_Sant_Pau/SBL_recon/ODF-FSL-dyads_format/Matrix_Dir.mat'];
    
    plot_max(path_to_mask, path_to_peaks, slice, type, invert_z, plot_axis_limit);
    box on; grid off;
        
    figure_name = ['/home/erick/Documents/Sparse_wars/Real_data_Sant_Pau/SBL_recon/ODF-FSL-dyads_format/Peaks_' type '_slice' num2str(slice)];
    export_fig (figure_name, '-jpg', '-nocrop', '-transparent', '-a1', '-r300', '-nocrop');
    close all
end

if plot_RFOCUSS == 1
    path_to_peaks = ['/home/erick/Documents/Sparse_wars/Real_data_Sant_Pau/RFOCUSS_recon/ODF-FSL-dyads_format/Matrix_Dir.mat'];
    
    plot_max(path_to_mask, path_to_peaks, slice, type, invert_z, plot_axis_limit);
    box on; grid off;
        
    figure_name = ['/home/erick/Documents/Sparse_wars/Real_data_Sant_Pau/RFOCUSS_recon/ODF-FSL-dyads_format/Peaks_' type '_slice' num2str(slice)];
    export_fig (figure_name, '-jpg', '-nocrop', '-transparent', '-a1', '-r300', '-nocrop');
    close all
end

if plot_LASSO == 1
    path_to_peaks = ['/home/erick/Documents/Sparse_wars/Real_data_Sant_Pau/LASSO_EBIC_recon/ODF-FSL-dyads_format/Matrix_Dir.mat'];
    
    plot_max(path_to_mask, path_to_peaks, slice, type, invert_z, plot_axis_limit);
    box on; grid off;
        
    figure_name = ['/home/erick/Documents/Sparse_wars/Real_data_Sant_Pau/LASSO_EBIC_recon/ODF-FSL-dyads_format/Peaks_' type '_slice' num2str(slice)];
    export_fig (figure_name, '-jpg', '-nocrop', '-transparent', '-a1', '-r300', '-nocrop');
    close all
end

if plot_ARUMBA == 1
    path_to_peaks = ['/home/erick/Documents/Sparse_wars/Real_data_Sant_Pau/ARUMBA_recon/ODF-FSL-dyads_format/Matrix_Dir.mat'];
    
    plot_max(path_to_mask, path_to_peaks, slice, type, invert_z, plot_axis_limit);
    box on; grid off;
        
    figure_name = ['/home/erick/Documents/Sparse_wars/Real_data_Sant_Pau/ARUMBA_recon/ODF-FSL-dyads_format/Peaks_' type '_slice' num2str(slice)];
    export_fig (figure_name, '-jpg', '-nocrop', '-transparent', '-a1', '-r300', '-nocrop');
    close all
end


if plot_IRL1 == 1
    path_to_peaks = ['/home/erick/Documents/Sparse_wars/Real_data_Sant_Pau/IRL1_recon/ODF-FSL-dyads_format/Matrix_Dir.mat'];
    
    plot_max(path_to_mask, path_to_peaks, slice, type, invert_z, plot_axis_limit);
    box on; grid off;
        
    figure_name = ['/home/erick/Documents/Sparse_wars/Real_data_Sant_Pau/IRL1_recon/ODF-FSL-dyads_format/Peaks_' type '_slice' num2str(slice)];
    export_fig (figure_name, '-jpg', '-nocrop', '-transparent', '-a1', '-r300', '-nocrop');
    close all
end

if plot_NNLS_EBI == 1
    path_to_peaks = ['/home/erick/Documents/Sparse_wars/Real_data_Sant_Pau/NNLS_BSS_EBIC_recon/ODF-FSL-dyads_format/Matrix_Dir.mat'];
    
    plot_max(path_to_mask, path_to_peaks, slice, type, invert_z, plot_axis_limit);
    box on; grid off;
        
    figure_name = ['/home/erick/Documents/Sparse_wars/Real_data_Sant_Pau/NNLS_BSS_EBIC_recon/ODF-FSL-dyads_format/Peaks_' type '_slice' num2str(slice)];
    export_fig (figure_name, '-jpg', '-nocrop', '-transparent', '-a1', '-r300', '-nocrop');
    close all
end

if plot_dRL_SD == 1
    path_to_peaks = ['/home/erick/Documents/Sparse_wars/Real_data_Sant_Pau/dRLSD_recon/ODF-FSL-dyads_format/Matrix_Dir.mat'];
    
    plot_max(path_to_mask, path_to_peaks, slice, type, invert_z, plot_axis_limit);
    box on; grid off;
        
    figure_name = ['/home/erick/Documents/Sparse_wars/Real_data_Sant_Pau/dRLSD_recon/ODF-FSL-dyads_format/Peaks_' type '_slice' num2str(slice)];
    export_fig (figure_name, '-jpg', '-nocrop', '-transparent', '-a1', '-r300', '-nocrop');
    close all
end

if plot_Cauchy == 1
    path_to_peaks = ['/home/erick/Documents/Sparse_wars/Real_data_Sant_Pau/Cauchy_recon/ODF-FSL-dyads_format/Matrix_Dir.mat'];
    
    plot_max(path_to_mask, path_to_peaks, slice, type, invert_z, plot_axis_limit);
    box on; grid off;
        
    figure_name = ['/home/erick/Documents/Sparse_wars/Real_data_Sant_Pau/Cauchy_recon/ODF-FSL-dyads_format/Peaks_' type '_slice' num2str(slice)];
    export_fig (figure_name, '-jpg', '-nocrop', '-transparent', '-a1', '-r300', '-nocrop');
    close all
end
