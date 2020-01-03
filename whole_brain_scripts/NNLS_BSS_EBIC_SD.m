function [] = NNLS_BSS_EBIC_SD(path_to_data, path_to_mask, b_vector, b_value, format_to_save_fODF, path_to_save_data, lambda1, lambda2, lambda_csf)

% Erick J. Canales-Rodríguez, Jon Haitz Legarreta-Gorroño, Yasser Alemán-Gomez, Marco Pizzolato, Jonathan Rafael Patiño, Gaëtan Olivier D Rensonnet, Muhamed Barakovic, Joaquim Radua, Edith Pomarol-Clotet, Raymond Salvador, Jean-Philippe Thiran, Alessandro Daducci. Sparse wars: A survey and comparative study of spherical deconvolution algorithms for diffusion MRI. Neuroimage. 2019 Jan 1;184:140-160.

%% ========================================================================
%        Loading the measurement scheme
%  ========================================================================

grad = load(b_vector);
if size(grad,2) > size(grad,1)
    disp(' - the gradient table is compatible with the FSL format')
    disp(' - it will be transformed (transposed) to our format')
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

display(['The measurement scheme contains ' num2str(size(grad,1)) ' volumes']);
display(['The max b-value is ' num2str(max(b))]);

% --- vector normalization: obtaining unit vectors (0-vectors are preserved)
normalize = 0;
if normalize == 1
   norm_factor = sqrt(sum(grad.^2,2));
   grad = grad./repmat(norm_factor + eps,[1 3]);
end

invertz = 0;
if invertz == 1
    grad(:,3) = -grad(:,3);
end

%% ========================================================================
%        Loading the data and mask + data correction
%  ========================================================================

Vdata = spm_vol(path_to_data);
E     = spm_read_vols(Vdata);

Vmask = spm_vol(path_to_mask);
Mask  = spm_read_vols(Vmask);
Mask(Mask<0.1) = 0;

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
E(E>1) = 1;
E(E<0) = 0;

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
        Vaux(i).dt    = [16 0];
        Vaux(i).n     = [i 1];
        Vaux(i).dim   = [n1 n2 n3];
        spm_write_vol(Vaux(i),squeeze(E(:,:,:,i)));
    end
    clear Vaux VE
    % ----------------------------------------------------
    %save('DTK-Trackvis_format/grad_for_DTK.txt','grad','-ASCII');
    save([path_to_save_data filesep 'DTK-Trackvis_format/grad_for_DTK.txt'],'grad','-ASCII');
    clear grad_for_DTK;
end

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
display(['The reconstruction scheme contains ' num2str(size(V,1)) ' directions']);

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

% -------------------------------------------------------------------------
% ----------------- Parameters for computing the peaks --------------------
% -------------------------------------------------------------------------
Vrec = V(length(V)/2 + 1:end,:);
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
    neig_i(neig_i == i) = 0;
    neig(i,1:length(neig_i)) = neig_i;
end
% -------------------------------------------------------------------------


%% - Creating the Kernel for the reconstruction. 
%   ---------------------------------
% Kernel is a Dictionary of diffusion tensor signals oriented along several
% directions fODF_ron the sphere (in our example V, which is the set of orientations where the fODF
% will be computed).
% It also includes two compartments with isotropic diffusion to fit real Brain data: CSF and GM

% diffusivities assumed for the Dictionary

%lambda1 = 1.3*10^-3;% mm^2/s
%lambda2 = 0.3*10^-3;% mm^2/s (Like the ball and stick model of T. Behrens)

noise_param.add_rician_noise = 0;
noise_param.SNR = [];

[phi, theta] = cart2sph(V(:,1),V(:,2),V(:,3)); % set of directions
theta = -theta;

Kernel = zeros(length(grad),length(V) );
S0 = 1; % The Dictionary is created under the assumption S0 = 1;
fi = 1; % volume fraction
for i=1:length(phi)
    anglesFi = [phi(i), theta(i)]*(180/pi); % in degrees
    Kernel(:,i) = create_signal_multi_tensor(anglesFi, fi, [lambda1, lambda2, lambda2], ...
                                             b, grad, S0, noise_param);
end
% adding the isotropic compartment for CSF
%lambda_csf = 2.5*10^-3;% mm^2/s
%lambda_csf = 1.0*10^-4;% mm^2/s
Kernel(:,length(phi) + 1) = create_signal_multi_tensor([0 0], fi, lambda_csf*[1, 1, 1], ...
                                             b, grad, S0, noise_param);
%   ---------------------------------


%% ========================================================================
%                     ODF-Reconstruction
%  ========================================================================

%% - Initial estimate of the S0 component from the signal
disp ('Starting the reconstruction ...')
disp(' ')

% Reconstruction based on only one hemisphere to optimize resources...
Kernel0 = Kernel;
if strcmp( format_to_save_fODF, 'DTK-Trackvis')
    Kernel0(:,1:181) = [];
elseif strcmp( format_to_save_fODF, 'ODF-FSL-dyads')
    Kernel0(:,1:362) = [];
elseif strcmp( format_to_save_fODF, 'DSI-Studio')
    Kernel0(:,1:321) = [];
end
fODF0 = ones(size(Kernel0,2),1)/size(Kernel0,2);

% threshold to identify peaks
threshold  = 0.1;
% max-number of peaks at each voxel
Num_fmax   = 5;
Matrix_Dir = zeros(n1, n2, n3, 3*Num_fmax);
Fiso       = zeros(n1, n2, n3);
Fwm        = zeros(n1, n2, n3);

poolobj = parpool(4);
% ----
%NumWorkers     = 6;
%obj            = parcluster;
%obj.NumWorkers = NumWorkers;
%poolobj = parpool(NumWorkers);

tic
parfor ii=1:n1
%for ii=1:n1
    %disp (['Slice ' num2str(ii) ' off ' num2str(n1)]);
    for jj= 1:n2
        for kk = 1:n3
            if Mask(ii,jj,kk) > 0
                Signal_ijk = E(ii,jj,kk,:);
                [Matrix_Dir_ijk, Fiso_ijk, Fw_ijk] = par_NNLS_BSS_EBIC(Signal_ijk(:), Kernel0, fODF0, Vrec, neig, threshold, Num_fmax);

                Matrix_Dir(ii, jj, kk,:) = Matrix_Dir_ijk;
                Fiso(ii,jj,kk,:) = Fiso_ijk;
                Fwm (ii,jj,kk,:) = Fw_ijk;
            end
        end
    end
end
toc

try
delete(poolobj);
delete(gcp('nocreate'))
end

%% ========================================================================
%                     Save Results
%  ========================================================================

display ('Save all Results');

if strcmp( format_to_save_fODF, 'DTK-Trackvis')
    path_to_save_data = [path_to_save_data filesep 'DTK-Trackvis_format'];
    mkdir(path_to_save_data);

    Viso = Vmask;
    Viso.fname = [path_to_save_data '/F_ISO.nii'];
    Viso.dt = [16 0];
    spm_write_vol(Viso,Fiso);
    
    Vwm = Vmask;
    Vwm.fname = [path_to_save_data '/F_WM.nii'];
    Vwm.dt = [16 0];
    spm_write_vol(Vwm,Fwm);
        
    save ([path_to_save_data '/Matrix_Dir'], 'Matrix_Dir');

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
    
    Vwm = Vmask;
    Vwm.fname = [path_to_save_data '/F_WM.nii'];
    Vwm.dt = [16 0];
    spm_write_vol(Vwm,Fwm);
    
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
    S = size(Matrix_Dir);
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
    
    save ([path_to_save_data '/Matrix_Dir'], 'Matrix_Dir');
    
    % ------------------ Save in DSI-Studio ------------------ %    
    voxel_size_vector = abs([Vmask.mat(1,1), Vmask.mat(2,2), Vmask.mat(3,3)]);
    rumba_sd2fib(Matrix_Dir, voxel_size_vector, Mask);
    % -------------------------------------------------------- %
end

end

% ==================================================================
%                     *** Private functions ***
% ==================================================================

function [Matrix_Dir_aux, fcsf, fw] = par_NNLS_BSS_EBIC(Signal, Kernel0, fODF0, Vrec, neig, threshold, Num_fmax)

% --- NNLS estimation
fODF_initial = lsqnonneg(Kernel0, Signal);
fODF = deconv_best_subset(Kernel0, Signal(:), fODF_initial);

% --- debias by NNLS
% ind_pos = fODF>0;
% ind_pos(end) = 1;
% fODF(ind_pos) = lsqnonneg(Kernel0(:,ind_pos), Signal(:));

% ----------------------------------
fcsf = fODF(end);    % fraction of the isotropic component for CSF
fODF(end) = [];      % purely anisotropic part

fw =  sum(fODF);

% =========================================================
%                  ODF Maxima Extraction
% =========================================================
Matrix_Dir_aux = zeros(3*Num_fmax, 1);

if sum(fODF(:)) > 0
    [peaks_dir, peaks_f, peaks] = find_discrete_peaks(Vrec, fODF, neig, threshold);
    peaks_ind = find(peaks);
    
    M = length(peaks_f);
    if M > Num_fmax
        M = Num_fmax;
        [peaks_f, ind_ord] = sort(peaks_f,'descend');
        peaks_dir = peaks_dir(ind_ord,:);
        peaks_ind = peaks_ind(ind_ord,:);
        % ---
        peaks_f   = peaks_f(1:Num_fmax);
        peaks_dir = peaks_dir(1:Num_fmax,:);
        peaks_ind = peaks_ind(1:Num_fmax);
    end
    
    peaks_f = peaks_f/(sum(peaks_f) + fcsf); % normalized volume fraction of each fiber.
    peaks_dirT = ([peaks_f, peaks_f, peaks_f].*peaks_dir)';
    Matrix_Dir_aux(1:3*M, 1) = peaks_dirT(:);
end
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
