function [] = RUMBA_sph_deconv(path_to_data, path_to_mask, b_vector, b_value, format_to_save_fODF, path_to_save_data, Niter, lambda1, lambda2, lambda_csf, lambda_gm, ncoils, MRI_recon_type, R, use_tv)
 
% Erick J Canales-Rodríguez, Alessandro Daducci, Stamatios N Sotiropoulos, Emmanuel Caruyer, Santiago Aja-Fernández, Joaquim Radua, Yasser Iturria-Medina, Lester Melie-García, Yasser Alemán-Gómez, Jean-Philippe Thiran, Salvador Sarró, Edith Pomarol-Clotet, Raymond Salvador. Spherical Deconvolution of Multichannel Diffusion MRI Data with Non-Gaussian Noise Models and Spatial Regularization. PLoS ONE, 2015, 10(10): e0138910. 

%% ========================================================================
%        Loading the measurement scheme
%  ========================================================================

grad = load(b_vector);
if size(grad,2) > size(grad,1)
    disp(' - the gradient table is compatible with the FSL format')
    disp(' - it will be transformed to our format')
    grad = grad';
else
    disp(' - the gradient table is compatible with our format')
end

if isfloat(b_value)
    b = b_value*ones(length(grad),1);
else
    b = load(b_value);
    b = b(:);
end

disp(['The measurement scheme contains ' num2str(size(grad,1)) ' volumes']);
display(['The max b-value is ' num2str(max(b))]);

% --- vector normalization: obtaining unit vectors (0-vectors are preserved)
% norm_factor = sqrt(sum(grad.^2,2));
% grad = grad./repmat(norm_factor + eps,[1 3]);

% --- save gradients and b-values in the FSL format
% grad_fsl = grad';
% save('bvecs_for_fsl.txt','grad_fsl','-ASCII');
% b_for_fsl = b';
% save('bvals_for_fsl.txt','b_for_fsl','-ASCII');

%% ========================================================================
%        Loading the data and mask + data correction
%  ========================================================================

Vdata = spm_vol(path_to_data);
E = spm_read_vols(Vdata);
% E = single(E); % convert to single to save memory in this dataset !!!
E(E<0) = eps;

Vmask = spm_vol(path_to_mask);
Mask = spm_read_vols(Vmask);
Mask(Mask>0) = 1.0;

ind_S0 = find(b.*sum(grad.^2,2) == 0);
disp(['From the ' num2str(size(grad,1)) ' volumes, ' num2str(length(ind_S0)) ' are b0 images']);

S0_est = squeeze( E(:,:,:,ind_S0) );
% case of several b0 images in the data
if length(ind_S0) > 1
    S0_est = mean(S0_est,4).*Mask;
end

% - mean b0 image first... (reordering the data)
E(:,:,:,ind_S0) = [];
E = cat(4,S0_est,E);

grad(ind_S0,:) = [];
grad = [0 0 0; grad];

b(ind_S0) = [];
b = [0; b];

Dim = size(E);
n1 = Dim(1);
n2 = Dim(2);
n3 = Dim(3);
n4 = Dim(4);

display(['The number of volumes in the 4D diffusion image is: ' num2str(n4)]);

for i = 1:n4
    E(:,:,:,i) = squeeze(E(:,:,:,i).*Mask)./(S0_est + eps);
end

% Repair signal
E(E>1)=1;

% --- crop the images in (x,y,z) to save time and memory --- %
% Don't worry, be happy, at the end the original image size will be preserved...
[xmin, xmax] = crop_ind(Mask,3,2);
[ymin, ymax] = crop_ind(Mask,3,1);
[zmin, zmax] = crop_ind(Mask,1,2);

E_red    = squeeze(E(xmin:xmax,ymin:ymax,zmin:zmax,:));
Mask_red = squeeze(Mask(xmin:xmax,ymin:ymax,zmin:zmax,:));
% ---------------------------------------------------------- %

if strcmp( format_to_save_fODF, 'DTK-Trackvis')
    % save new signal and grad in the DTK format
    % mkdir('DTK-Trackvis_format');
    mkdir([path_to_save_data filesep 'DTK-Trackvis_format'])
    % ----------------------------------------------------
    Vaux = Vdata(1:n4,:);
    Vaux = rmfield(Vaux,'private');
    % ---
    for i=1:n4
        %Vaux(i).fname = 'DTK-Trackvis_format/Signal_for_DTK.nii';
        Vaux(i).fname = [path_to_save_data filesep 'DTK-Trackvis_format/Signal_for_DTK.nii'];
        Vaux(i).dt = [16 0];
        Vaux(i).n = [i 1];
        Vaux(i).dim = [n1 n2 n3];
        spm_write_vol(Vaux(i),squeeze(E(:,:,:,i)));
    end
    clear Vaux VE
    % ----------------------------------------------------
    %save('DTK-Trackvis_format/grad_for_DTK.txt','grad','-ASCII');
    save([path_to_save_data filesep 'DTK-Trackvis_format/grad_for_DTK.txt'],'grad','-ASCII');
    clear grad_for_DTK;
end
clear E Mask

%% ========================================================================
%        Loading the reconstruction scheme
%  ========================================================================

if strcmp( format_to_save_fODF, 'DTK-Trackvis')
    V = load('sampling_and_reconstruction_schemes/On_the_sphere/362_shell_trackvis.txt');
    Vrec = V(1:181,:);
    % Obtaining the faces on one hemisphere, but considering the symmetry
    % about the plane: by construction in this scheme V(k,:) = -V(181+k,:) for all k=1:181;
    % The faces are used to extract the ODF maxima.
    Vfaces = convhulln(V);
elseif strcmp( format_to_save_fODF, 'ODF-FSL-dyads')
    %V = load('sampling_and_reconstruction_schemes/On_the_sphere/724_shell.txt');
    %V = load('sampling_and_reconstruction_schemes/On_the_sphere/vertices724.txt'); % *** Gab sampling scheme ***
    V = load('sampling_and_reconstruction_schemes/On_the_sphere/repulsion724.txt'); % *** Gab sampling scheme *** 

    Vrec = V(1:362,:);
    % Obtaining the faces on one hemisphere, but considering the symmetry
    % about the plane: by construction in this scheme V(k,:) = -V(362+k,:) for all k=1:362;
    % The faces are used to extract the ODF maxima.
    Vfaces = convhulln(V);
elseif strcmp( format_to_save_fODF, 'DSI-Studio')
    Vstruc = load('sampling_and_reconstruction_schemes/On_the_sphere/odf8.mat');
    V = Vstruc.odf_vertices';
    Vrec = V(1:321,:);
    % Obtaining the faces on one hemisphere, but considering the symmetry
    % about the plane: by construction in this scheme V(k,:) = -V(321+k,:) for all k=1:321;
    % The faces are used to extract the ODF maxima.
    Vfaces = Vstruc.odf_faces;
end
disp(['The reconstruction scheme contains ' num2str(size(V,1)) ' directions']);

% F = convhulln(V); command to obtain the facets from any of these spherical grids

%% ========================================================================
%             Precomputation for fODF maxima extraction
%  ========================================================================

if strcmp( format_to_save_fODF, 'DTK-Trackvis')
    ang_neig = 16;
elseif strcmp( format_to_save_fODF, 'ODF-FSL-dyads')
    ang_neig = 12.5;
elseif strcmp( format_to_save_fODF, 'DSI-Studio')
    ang_neig = 12.5;
end

% ang_neig = 12.5;

% --- parameters to find_discrete_peaks
cos_ang = Vrec*Vrec';
% --- correction due to numerical errors
cos_ang(cos_ang >  1)  =  1;
cos_ang(cos_ang < -1)  = -1;
% --------------------------------------
ang = acosd(cos_ang);
% --- considering the antipodal symmetry
ang = min(ang,180-ang);
% --------------------------------------
% --- select all dirs within "ang_neig" degrees
for i = 1:length(Vrec)
    ang_vect = ang(i,:);
    neig_i = find(ang_vect < ang_neig);
    neig(i,1:length(neig_i)) = neig_i;
end

%% ========================================================================
% - Creating the Kernel for the reconstruction. 
%  ========================================================================

% Kernel is a Dictionary of diffusion tensor signals oriented along several
% directions on the sphere (in our example V, which is the set of orientations where the fODF
% will be computed).
% It also includes two compartments with isotropic diffusion to fit real Brain data: CSF and GM

[phi, theta] = cart2sph(V(:,1),V(:,2),V(:,3)); % set of directions
theta = -theta;

Kernel = zeros(length(grad),length(V) );
S0 = 1; % The Dictionary is created under the assumption S0 = 1;
fi = 1; % volume fraction
for i=1:length(phi)
    anglesFi = [phi(i), theta(i)]*(180/pi); % in degrees
    Kernel(:,i) = create_signal_multi_tensor(anglesFi, fi, [lambda1, lambda2, lambda2], b, grad, S0);
end
% adding the isotropic compartment for CSF
Kernel(:,length(phi) + 1) = create_signal_multi_tensor([0 0], fi, [lambda_csf, lambda_csf, lambda_csf], b, grad, S0);
% adding the isotropic compartment for GM
Kernel(:,length(phi) + 2) = create_signal_multi_tensor([0 0], fi, [lambda_gm, lambda_gm, lambda_gm], b, grad, S0);
                                         

%% ========================================================================
%                     ODF-Reconstruction
%  ========================================================================
disp (' ');
disp ('fODF-Reconstruction');

% -- Initial solution: Constant value on the sphere, non-negative iso-probable spherical function
fODF0 = ones(size(Kernel,2),1);
fODF0 = fODF0/sum(fODF0);

Kernel0 = Kernel;
fODF00 = fODF0;

% Reconstruction based on only one hemisphere to optimize resources...
if strcmp( format_to_save_fODF, 'DTK-Trackvis')
    Kernel0(:,1:181) = [];
    fODF00(1:181) = [];
elseif strcmp( format_to_save_fODF, 'ODF-FSL-dyads')
    Kernel0(:,1:362) = [];
    fODF00(1:362) = [];
elseif strcmp( format_to_save_fODF, 'DSI-Studio')
    Kernel0(:,1:321) = [];
    fODF00(1:321) = [];
end

fODF00 = fODF00/sum(fODF00);

tic
[fODF_4D, fgm_3D, fcsf_3D, SNR_est, std_SNR, var_3D] = sph_deconv_tv_motor_all(E_red, Kernel0, fODF00, Niter, Mask_red, ncoils, MRI_recon_type, R, use_tv);
toc

SNR_est = double(SNR_est);
std_SNR = double(std_SNR);

% --- restoring the images to their original size --- %
fODF_4D_aux = zeros(n1,n2,n3,size(fODF_4D,4));
fODF_4D_aux(xmin:xmax,ymin:ymax,zmin:zmax,:) = fODF_4D;
fODF_4D = fODF_4D_aux;
clear fODF_4D_aux;

fgm_3D_aux = zeros(n1,n2,n3);
fgm_3D_aux(xmin:xmax,ymin:ymax,zmin:zmax) = fgm_3D;
fgm_3D = fgm_3D_aux;
clear fgm_3D_aux;

fcsf_3D_aux = zeros(n1,n2,n3);
fcsf_3D_aux(xmin:xmax,ymin:ymax,zmin:zmax) = fcsf_3D;
fcsf_3D = fcsf_3D_aux;
clear fcsf_3D_aux;

var_3D_aux = zeros(n1,n2,n3);
var_3D_aux(xmin:xmax,ymin:ymax,zmin:zmax) = var_3D;
var_3D = var_3D_aux;
clear var_3D_aux
% --------------------------------------------------- %

Fiso = (fgm_3D + fcsf_3D);
fwm_3D = sum(fODF_4D,4);

%% ========================================================================
%                     ODF Maxima Extraction
%  ========================================================================

disp ('ODF Maxima Extraction');

% max number of peaks at each voxel
Num_fmax = 5;
fract_max = 0.1; % threshold to consider a fiber component as real
[Matrix_Dir]= obtain_peaks_info(fODF_4D, V, neig, fract_max, format_to_save_fODF, Fiso, Num_fmax);

%% ========================================================================
%                   Adding the isotropic components to the ODF
%  ========================================================================

% for x = 1:n1
%     for y = 1:n2
%         for z = 1:n3
%             fODF_b = squeeze(fODF_4D(x,y,z,:));            
%             % ----------------------------------------------------------- %
%             fiso = Fiso(x,y,z);
%             fODF_b = fODF_b + fiso/length(fODF_b); % iso+ani
%             fODF_b = fODF_b/( sum(fODF_b) + eps ); % (eps -> avoid NaN if fODF_b = 0)
%             % ----------------------------------------------------------- %
%             fODF_4D(x,y,z,:) = fODF_b;
%             % ----------------------------------------------------------- %
%         end
%     end
% end

fiso = repmat(Fiso/size(fODF_4D,4),[1 1 1 size(fODF_4D,4)]);
fODF_4D = fODF_4D + fiso;
clear fiso;
normf = repmat(sum(fODF_4D,4),[1 1 1 size(fODF_4D,4)]);
fODF_4D = fODF_4D./(normf + eps);
clear normf;

%% ========================================================================
%                     Save Results
%  ========================================================================

disp ('Save all Results');

if strcmp( format_to_save_fODF, 'DTK-Trackvis')
    path_to_save_data = [path_to_save_data filesep 'DTK-Trackvis_format'];
    mkdir(path_to_save_data);

    Viso = Vmask;
    Viso.fname = [path_to_save_data '/F_ISO.nii'];
    Viso.dt = [16 0];
    spm_write_vol(Viso,Fiso);
    
    Vcsf = Vmask;
    Vcsf.fname = [path_to_save_data '/F_CSF.nii'];
    Vcsf.dt = [16 0];
    spm_write_vol(Vcsf,fcsf_3D);
    
    Vgm = Vmask;
    Vgm.fname = [path_to_save_data '/F_GM.nii'];
    Vgm.dt = [16 0];
    spm_write_vol(Vgm,fgm_3D);
    
    Vwm = Vmask;
    Vwm.fname = [path_to_save_data '/F_WM.nii'];
    Vwm.dt = [16 0];
    spm_write_vol(Vwm,fwm_3D);
    
    Vvar = Vmask;
    Vvar.fname = [path_to_save_data '/variance_3d.nii'];
    Vvar.dt = [16 0];
    spm_write_vol(Vvar,var_3D);
    
    save ([path_to_save_data '/SNR_est'], 'SNR_est','std_SNR', '-ASCII');
    
    save ([path_to_save_data '/Matrix_Dir'], 'Matrix_Dir');

    % ---- Scalar metrics ---- %
    % --- gfa
    rms = sqrt(mean(fODF_4D.^2,4));
    stdev = std(fODF_4D,0,4);
    gfa = stdev./(rms + eps);
    
    
    Vgfa = Vmask;
    Vgfa.fname = [path_to_save_data '/gfa.nii'];
    Vgfa.dt = [16 0];
    spm_write_vol(Vgfa,gfa);

    display ('Voxelwise Standard HARDI-Reconstruction via "DTK"');
    try
        % --- Preparing the data to compute the ODF with DTK & Trackvis
        % --- The resulting hardi_odf and hardi_max images are used to save
        %     our fiber ODF estimates using the same .nii format. At the end, the original files
        %     produced by DTK can be deleted...
        % -------------------------------------------------------------
        system(['hardi_mat ' path_to_save_data '/grad_for_DTK.txt ' path_to_save_data '/temp_mat.dat -ref ' path_to_data ' -iop 1 0 0 0 1 0']);
        system(['odf_recon ' path_to_data ' ' num2str(size(grad,1)) ' 181 ' path_to_save_data '/hardi -b0 1 -iop 1 0 0 0 1 0 -mat ' path_to_save_data '/temp_mat.dat -p 3 -sn 1 -ot nii']);
        system(['odf_tracker ' path_to_save_data '/hardi ' path_to_save_data '/hardi_raw.trk -at 35 -m ' path_to_mask ' -it nii']);
        system(['spline_filter ' path_to_save_data '/hardi_raw.trk 1 ' path_to_save_data '/hardi.trk']);
        
        % save fODF to .nii
        % -------------------
        ODF_matrix = fODF_4D(:,:,:,1:181);
        
        ODFname = [path_to_save_data '/sphdeconv_lr_odf.nii'];
        
        data = [path_to_save_data '/hardi_odf.nii'];
        Vxx = spm_vol(data);
        Vtkv = Vxx;
        
        for i=1:n3
            Vtkv(i).fname = ODFname;
            Vtkv(i).private.dat.fname = ODFname;
        end
        
        ODF_rotated = permute(ODF_matrix,[4 1 2 3]);
        
        for i = 1:n3
            disp(i);
            Vtkv_i = Vtkv(i);
            Vtkv_i = spm_create_vol(Vtkv_i);
            spm_write_plane(Vtkv_i,squeeze( ODF_rotated(:,:,:,i) ),':');
        end
        clear Vxx Vtkv Vtkv_i;
        
        % save maxima to .nii
        % -------------------
        
        ODFname = [path_to_save_data '/sphdeconv_lr_max.nii'];
        
        data = [path_to_save_data '/hardi_max.nii'];
        Vxx = spm_vol(data);
        Vtkv = Vxx;
        
        for i=1:n3
            Vtkv(i).fname = ODFname;
            Vtkv(i).private.dat.fname = ODFname;
        end
        
        ODF_matrix_peaks_rotated = permute(ODF_matrix_peaks,[4 1 2 3]);
        
        for i = 1:n3
            disp(i);
            Vtkv_i = Vtkv(i);
            Vtkv_i = spm_create_vol(Vtkv_i);
            spm_write_plane(Vtkv_i,squeeze( ODF_matrix_peaks_rotated(:,:,:,i) ),':');
        end
        clear Vxx Vtkv Vtkv_i;
        
        % fiber tracking
        % -------------------
        system(['odf_tracker ' path_to_save_data '/sphdeconv_lr ' path_to_save_data '/sphdeconv_lr_raw.trk -at 45 -iz -m ' path_to_mask ' -it nii']);
        system(['spline_filter ' path_to_save_data '/sphdeconv_lr_raw.trk 1 ' path_to_save_data '/sphdeconv_lr.trk']);
        % --- delete some files ----
        delete([path_to_save_data '/hardi.trk']);
        delete([path_to_save_data '/hardi_raw.trk']);
        delete([path_to_save_data '/hardi_b0.nii']);
        delete([path_to_save_data '/hardi_dwi.nii']);
        delete([path_to_save_data '/hardi_odf.nii']);
        delete([path_to_save_data '/hardi_max.nii']);
        delete([path_to_save_data '/Signal_for_DTK.nii']);
    end
end

if strcmp( format_to_save_fODF, 'ODF-FSL-dyads')
    path_to_save_data = [path_to_save_data filesep 'ODF-FSL-dyads_format'];
    mkdir(path_to_save_data);
        
    Viso = Vmask;
    Viso.fname = [path_to_save_data '/F_ISO.nii'];
    Viso.dt = [16 0];
    spm_write_vol(Viso,Fiso);
    
    Vcsf = Vmask;
    Vcsf.fname = [path_to_save_data '/F_CSF.nii'];
    Vcsf.dt = [16 0];
    spm_write_vol(Vcsf,fcsf_3D);
    
    Vgm = Vmask;
    Vgm.fname = [path_to_save_data '/F_GM.nii'];
    Vgm.dt = [16 0];
    spm_write_vol(Vgm,fgm_3D);
    
    Vwm = Vmask;
    Vwm.fname = [path_to_save_data '/F_WM.nii'];
    Vwm.dt = [16 0];
    spm_write_vol(Vwm,fwm_3D);
    
    Vvar = Vmask;
    Vvar.fname = [path_to_save_data '/variance_3d.nii'];
    Vvar.dt = [16 0];
    spm_write_vol(Vvar,var_3D);
    
    save ([path_to_save_data '/SNR_est'], 'SNR_est','std_SNR', '-ASCII');
    
    % ---- Scalar metrics ---- %
    % --- gfa
    rms = sqrt(mean(fODF_4D.^2,4));
    stdev = std(fODF_4D,0,4);
    gfa = stdev./(rms + eps);
    
    Vgfa = Vmask;
    Vgfa.fname = [path_to_save_data '/gfa.nii'];
    Vgfa.dt = [16 0];
    spm_write_vol(Vgfa,gfa);

    
    Vdata_red = Vdata(1:3,:);
    Vfn = Vdata_red;
    
    % ---- Separate peaks ---- %
    for j=1:Num_fmax
        ODFname = [path_to_save_data '/Peaks_f' num2str(j-1) '.nii'];
        for i=1:3
            Vfn(i).fname = ODFname;
            Vfn(i).private.dat.fname = ODFname;
            Vfn(i).dt = [16 0];
        end
        
        nf = 3*j;
        nf_ind = nf-2:nf;
        for i = 1:3
            disp(i);
            Vtkv_i = Vfn(i);
            Vtkv_i = spm_create_vol(Vtkv_i);
            spm_write_plane(Vtkv_i,squeeze( Matrix_Dir(:,:,:,nf_ind(i)) ),':');
        end
    end
    clear Vfn Vtkv_i;
    
    % ---- 4D all peaks ---- %
    Vaux = Vdata(1:Num_fmax*3,:);
    Vaux = rmfield(Vaux,'private');
    S = size(fODF_4D);
    for i=1:Num_fmax*3
        VE(i) = Vaux(1);
        VE(i).fname = [path_to_save_data '/all_peaks.nii'];
        VE(i).dt = [16 0];
        VE(i).n = [i 1];
        VE(i).dim = [S(1) S(2) S(3)];
        spm_write_vol(VE(i),squeeze(Matrix_Dir(:,:,:,i)));
    end
    clear VE Vaux;
    
    gzip([path_to_save_data '/all_peaks.nii']);
    
    save ([path_to_save_data '/Matrix_Dir'], 'Matrix_Dir');
    % ----------------------------------
    
    Vaux = Vdata;
    Vaux = rmfield(Vaux,'private');
    S = size(fODF_4D);
    % ---
    for i=1:S(4)
        VE(i) = Vaux(1);
        VE(i).fname = [path_to_save_data '/ODF_matrix_4D.nii'];
        VE(i).dt = [16 0];
        VE(i).n = [i 1];
        VE(i).dim = [S(1) S(2) S(3)];
        spm_write_vol(VE(i),squeeze(fODF_4D(:,:,:,i)));
    end
    
end

if strcmp( format_to_save_fODF, 'DSI-Studio')
    
    path_to_save_data = [path_to_save_data filesep 'DSI-Studio_format'];
    mkdir(path_to_save_data);
    
    % ------------------------------------------------ %
    % ---- save in my format to visualize the ODFs --- %
    
    Viso = Vmask;
    Viso.fname = [path_to_save_data '/F_ISO.nii'];
    Viso.dt = [16 0];
    spm_write_vol(Viso,Fiso);
    
    Vcsf = Vmask;
    Vcsf.fname = [path_to_save_data '/F_CSF.nii'];
    Vcsf.dt = [16 0];
    spm_write_vol(Vcsf,fcsf_3D);
    
    Vgm = Vmask;
    Vgm.fname = [path_to_save_data '/F_GM.nii'];
    Vgm.dt = [16 0];
    spm_write_vol(Vgm,fgm_3D);
    
    Vwm = Vmask;
    Vwm.fname = [path_to_save_data '/F_WM.nii'];
    Vwm.dt = [16 0];
    spm_write_vol(Vwm,fwm_3D);
    
    Vvar = Vmask;
    Vvar.fname = [path_to_save_data '/variance_3d.nii'];
    Vvar.dt = [16 0];
    spm_write_vol(Vvar,var_3D);
    
    save ([path_to_save_data '/SNR_est'], 'SNR_est','std_SNR', '-ASCII');
    
    Vaux = Vdata;
    Vaux = rmfield(Vaux,'private');
    S = size(fODF_4D);
    % ---
    for i=1:S(4)
        VE(i) = Vaux(1);
        VE(i).fname = [path_to_save_data '/ODF_matrix_4D.nii'];
        VE(i).dt = [16 0];
        VE(i).n = [i 1];
        VE(i).dim = [S(1) S(2) S(3)];
        spm_write_vol(VE(i),squeeze(fODF_4D(:,:,:,i)));
    end
    
    % ---- Scalar metrics ---- %
    % --- gfa
    rms = sqrt(mean(fODF_4D.^2,4));
    stdev = std(fODF_4D,0,4);
    gfa = stdev./(rms + eps);
    
    Vgfa = Vmask;
    Vgfa.fname = [path_to_save_data '/gfa.nii'];
    Vgfa.dt = [16 0];
    spm_write_vol(Vgfa,gfa);
    
    save ([path_to_save_data '/Matrix_Dir'], 'Matrix_Dir');

    % ------------------ Save in DSI-Studio ------------------ %    
    voxel_size_vector = abs([Vmask.mat(1,1), Vmask.mat(2,2), Vmask.mat(3,3)]);
    rumba_sd2fib(Matrix_Dir, voxel_size_vector, gfa);
    % -------------------------------------------------------- %
end

end

% ==================================================================
%                     *** Private functions ***
% ==================================================================

%function [ODF_matrix_peaks, Matrix_Dir] = obtain_peaks_info(fODF_4D, V, neig, thr, format_to_save_fODF, Fiso, Num_fmax)
function [Matrix_Dir] = obtain_peaks_info(fODF_4D, V, neig, thr, format_to_save_fODF, Fiso, Num_fmax)


%% --------------------------------------
% fODF_4D = cat(4,fODF_4D,fODF_4D);
Dim = size(fODF_4D);

Matrix_Dir = zeros(Dim(1),Dim(2),Dim(3), 3*Num_fmax);

for x = 1:Dim(1)
    for y = 1:Dim(2)
        for z = 1:Dim(3)
            % ----------------------------------------------------------------------------------- %
            fiso = Fiso(x,y,z);
            %thr_xyz = thr + fiso;
            thr_xyz = thr/(1-fiso); %higher threshold for voxels with high fiso
            
            if strcmp( format_to_save_fODF, 'DTK-Trackvis')
                ODF = squeeze(fODF_4D(x,y,z,1:181));
                Vred = V(1:181,:);
            elseif strcmp( format_to_save_fODF, 'ODF-FSL-dyads')
                ODF = squeeze(fODF_4D(x,y,z,1:362));
                Vred = V(1:362,:);
            elseif strcmp( format_to_save_fODF, 'DSI-Studio')
                ODF = squeeze(fODF_4D(x,y,z,1:321));
                Vred = V(1:321,:);
            end
            
            if sum(ODF) > 0
                [peaks_dir, peaks_f, peaks] = find_discrete_peaks(Vred, ODF, neig, thr_xyz);
                
                peaks_ind = find(peaks);
                
                % ---
                M = length(peaks_f);
                if M > Num_fmax
                    M = Num_fmax;
                    [peaks_f, ind_ord] = sort(peaks_f,'descend');
                    peaks_dir = peaks_dir(ind_ord,:);
                    peaks_ind = peaks_ind(ind_ord,:);
                    % ---
                    peaks_f = peaks_f(1:Num_fmax);
                    peaks_dir = peaks_dir(1:Num_fmax,:);
                    peaks_ind = peaks_ind(1:Num_fmax);
                end
                
                peaks_f = (1-fiso)*( peaks_f/sum(peaks_f) ); % real volume fraction of each fiber.
                % (This takes into account the volume fraction of the isotropic terms)
                
                if strcmp( format_to_save_fODF, 'DTK-Trackvis')
                    ODF_matrix_peaks(x,y,z,peaks_ind) = 1;
                end
                                
                peaks_dirT = ([peaks_f, peaks_f, peaks_f].*peaks_dir)';
                Matrix_Dir(x,y,z,1:3*M) = peaks_dirT(:);
            end
            % ----------------------------------------------------------------------------------- %
        end
    end
end
end


function [peaks_dir, peaks_f, peaks] = find_discrete_peaks(V, ODF, neig, thr)

%% --- discrete maxima ---
peaks = zeros(size(ODF));
thr_value = thr*max(ODF);
for i = 1:length(ODF)
    neig_i = neig(i,:);
    ind_valid_i =  (neig_i ~= i) & (neig_i > 0) ;
    neig_final = neig_i(ind_valid_i);
    if ( ODF(i) > max(ODF(neig_final)) ) && ( ODF(i) >= thr_value )
        peaks(i) = 1;
    end
end
peaks_dir = V(find(peaks),:);
peaks_f = ODF(find(peaks));
end


function [amin, amax] = crop_ind(Data, dim1, dim2)
I = sum(sum(Data,dim1),dim2);
ind = find( I > 0);

amin = min(ind);
if amin > 1
    amin = amin - 1;
end

amax = max(ind);
if amax < size(I)
    amax = amax + 1;
end
end


function [S, D] = create_signal_multi_tensor (ang, f, Eigenvalues, b, grad, S0)

A = diag(Eigenvalues);

S = 0;
Nfibers = length(f);
f = f/sum(f);
for i = 1:Nfibers
    phi(i) = ang(i, 1);
    theta(i) = ang(i, 2);
    R = RotMatrix(phi(i),theta(i));
    D = R*A*R';
    S = S + f(i)*exp(-b.*diag(grad*D*grad'));
end
S = S0*S;
end

function R = RotMatrix(phi,theta)

c = pi/180;
phi = phi*c;
theta = theta*c;

Rz = [ cos(phi)  -sin(phi)  0
       sin(phi)   cos(phi)  0
           0         0      1];


Ry = [cos(theta)   0   sin(theta)
          0        1         0
     -sin(theta)   0   cos(theta)];

R =  Rz*Ry;
end


function rumba_sd2fib(Matrix_Dir, voxel_size_vector, gfa)
% Function provided by Ping-Hong Yeh

% Assuming 5 fibers per voxel...

Matrix_Dir(isnan(Matrix_Dir)) = 0;
fib.dimension = size(Matrix_Dir(:,:,:,1));
fib.voxel_size = voxel_size_vector;
fib.dir0 = Matrix_Dir(:,:,:,1:3);
fib.dir1 = Matrix_Dir(:,:,:,4:6);
fib.dir2 = Matrix_Dir(:,:,:,7:9);
fib.dir3 = Matrix_Dir(:,:,:,10:12);
fib.dir4 = Matrix_Dir(:,:,:,13:15);

fib.fa0 = sqrt(fib.dir0(:,:,:,1).*fib.dir0(:,:,:,1) + fib.dir0(:,:,:,2).*fib.dir0(:,:,:,2) + fib.dir0(:,:,:,3).*fib.dir0(:,:,:,3));
fib.fa1 = sqrt(fib.dir1(:,:,:,1).*fib.dir1(:,:,:,1) + fib.dir1(:,:,:,2).*fib.dir1(:,:,:,2) + fib.dir1(:,:,:,3).*fib.dir1(:,:,:,3));
fib.fa2 = sqrt(fib.dir2(:,:,:,1).*fib.dir2(:,:,:,1) + fib.dir2(:,:,:,2).*fib.dir2(:,:,:,2) + fib.dir2(:,:,:,3).*fib.dir2(:,:,:,3));
fib.fa3 = sqrt(fib.dir3(:,:,:,1).*fib.dir3(:,:,:,1) + fib.dir3(:,:,:,2).*fib.dir3(:,:,:,2) + fib.dir3(:,:,:,3).*fib.dir3(:,:,:,3));
fib.fa4 = sqrt(fib.dir4(:,:,:,1).*fib.dir4(:,:,:,1) + fib.dir4(:,:,:,2).*fib.dir4(:,:,:,2) + fib.dir4(:,:,:,3).*fib.dir4(:,:,:,3));

max_fib_fa0 = max(reshape(fib.fa0,1,[]));
mask0 = fib.fa0 < 0.01*max_fib_fa0; 
mask1 = fib.fa1 < 0.01*max_fib_fa0;  
mask2 = fib.fa2 < 0.01*max_fib_fa0;  
mask3 = fib.fa3 < 0.01*max_fib_fa0;  
mask4 = fib.fa4 < 0.01*max_fib_fa0;  

fib.fa0(mask0) = 1.0;
fib.fa1(mask0) = 1.0;
fib.fa2(mask0) = 1.0;
fib.fa3(mask0) = 1.0;
fib.fa4(mask0) = 1.0;

fib.dir0(:,:,:,1) = fib.dir0(:,:,:,1)./fib.fa0;
fib.dir0(:,:,:,2) = fib.dir0(:,:,:,2)./fib.fa0;
fib.dir0(:,:,:,3) = fib.dir0(:,:,:,3)./fib.fa0;
fib.dir1(:,:,:,1) = fib.dir1(:,:,:,1)./fib.fa1;
fib.dir1(:,:,:,2) = fib.dir1(:,:,:,2)./fib.fa1;
fib.dir1(:,:,:,3) = fib.dir1(:,:,:,3)./fib.fa1;
fib.dir2(:,:,:,1) = fib.dir2(:,:,:,1)./fib.fa2;
fib.dir2(:,:,:,2) = fib.dir2(:,:,:,2)./fib.fa2;
fib.dir2(:,:,:,3) = fib.dir2(:,:,:,3)./fib.fa2;
fib.dir3(:,:,:,1) = fib.dir3(:,:,:,1)./fib.fa3;
fib.dir3(:,:,:,2) = fib.dir3(:,:,:,2)./fib.fa3;
fib.dir3(:,:,:,3) = fib.dir3(:,:,:,3)./fib.fa3;
fib.dir4(:,:,:,1) = fib.dir4(:,:,:,1)./fib.fa4;
fib.dir4(:,:,:,2) = fib.dir4(:,:,:,2)./fib.fa4;
fib.dir4(:,:,:,3) = fib.dir4(:,:,:,3)./fib.fa4;


fib.fa0(mask0) = 0.0;
fib.fa1(mask0) = 0.0;
fib.fa2(mask0) = 0.0;
fib.fa3(mask0) = 0.0;
fib.fa4(mask0) = 0.0;

fib.fa0 = reshape(fib.fa0,1,[]);
fib.fa1 = reshape(fib.fa1,1,[]);
fib.fa2 = reshape(fib.fa2,1,[]);
fib.fa3 = reshape(fib.fa3,1,[]);
fib.fa4 = reshape(fib.fa4,1,[]);

% ----------------------------
gfa(isnan(gfa)) = 0;
fib.gfa = reshape(gfa,1,[]);
% ----------------------------

fib.dir0 = permute(fib.dir0,[4 1 2 3]);
fib.dir1 = permute(fib.dir1,[4 1 2 3]);
fib.dir2 = permute(fib.dir2,[4 1 2 3]);
fib.dir3 = permute(fib.dir3,[4 1 2 3]);
fib.dir4 = permute(fib.dir4,[4 1 2 3]);

fib.dir0 = reshape(fib.dir0,3,[]);
fib.dir1 = reshape(fib.dir1,3,[]);
fib.dir2 = reshape(fib.dir2,3,[]);
fib.dir3 = reshape(fib.dir3,3,[]);
fib.dir4 = reshape(fib.dir4,3,[]);

fib.dir0(isnan(fib.dir0)) = 0;
fib.dir1(isnan(fib.dir1)) = 0;
fib.dir2(isnan(fib.dir2)) = 0;
fib.dir3(isnan(fib.dir3)) = 0;
fib.dir4(isnan(fib.dir4)) = 0;

clear mask0;
clear mask1;
clear mask2;
clear mask3;
clear mask4;

save('out.fib','-struct', 'fib', '-v4');
gzip('out.fib');
delete('out.fib');
end
