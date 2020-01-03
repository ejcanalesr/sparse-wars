% Erick J. Canales-Rodríguez, Jon Haitz Legarreta-Gorroño, Yasser Alemán-Gomez, Marco Pizzolato, Jonathan Rafael Patiño, Gaëtan Olivier D Rensonnet, Muhamed Barakovic, Joaquim Radua, Edith Pomarol-Clotet, Raymond Salvador, Jean-Philippe Thiran, Alessandro Daducci. Sparse wars: A survey and comparative study of spherical deconvolution algorithms for diffusion MRI. Neuroimage. 2019 Jan 1;184:140-160.

% Script to plot the fODF in a brain slice
clear; clc
path_to_mask = ['brain_mask.nii'];

% For example ...
folder  = {'/home/erick/Documents/Sparse_wars/Real_data_Sant_Pau/IRL1_recon/ODF-FSL-dyads_format/'};


%_________________________________________________________________________%

plot_fODF = 'yes'; plot_fODFmax = 'no';
for i=1:length(folder)
    folder_i = char(folder(i));
    display (['Creating figure for path ' num2str(i)]);
    display (['Folder -> ' folder_i]);
    display (' ');
    % ---------------------------------------------------------------------
    path_to_ODF =    [folder_i '/ODF_matrix_4D.nii'];
    path_to_Matdir = [folder_i '/Matrix_Dir.mat'];
    %path_to_Matdir = [folder_i '/peaks.nii'];

    
    %path_to_FWM =    [folder_i '/data-vf_wm.nii'];
    %path_to_GFA =    [folder_i '/data-vf_gfa.nii'];
    
    path_to_FWM =    [folder_i '/F_WM.nii'];
    path_to_GFA =    ['dti_fa.nii'];

    % ---------------------------------------------------------------------
    type = 'y'; slice = 74;
    
    x1 = 38; x2 = 89;
    y1 = 15; y2 = 48;

    plot_ODFs_SF(path_to_ODF, path_to_mask, path_to_GFA, path_to_FWM, path_to_Matdir, slice, type, plot_fODF, plot_fODFmax);

    axis([x1 x2 slice-2 slice+2 y1 y2]); box on; axis off
    line([x1 x2],[slice slice],[y1 y1],'LineWidth',12,'Color',[.2 .2 .8]')
    line([x1 x2],[slice slice],[y2 y2],'LineWidth',12,'Color',[.2 .2 .8]')
    line([x1 x1],[slice slice],[y1 y2],'LineWidth',12,'Color',[.2 .2 .8]')
    line([x2 x2],[slice slice],[y1 y2],'LineWidth',12,'Color',[.2 .2 .8]')
    %figure_name = [folder_i 'plot_ODF_roi_slice_' type '_' type num2str(slice)];
    %export_fig (figure_name, '-jpg', '-nocrop', '-transparent', '-a1', '-r300', '-nocrop');
    %close all;
end
