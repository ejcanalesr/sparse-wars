% Erick J. Canales-Rodríguez, Jon Haitz Legarreta-Gorroño, Yasser Alemán-Gomez, Marco Pizzolato, Jonathan Rafael Patiño, Gaëtan Olivier D Rensonnet, Muhamed Barakovic, Joaquim Radua, Edith Pomarol-Clotet, Raymond Salvador, Jean-Philippe Thiran, Alessandro Daducci. Sparse wars: A survey and comparative study of spherical deconvolution algorithms for diffusion MRI. Neuroimage. 2019 Jan 1;184:140-160.

function [] = plot_ODFs_SF(path_to_ODF, path_to_mask, path_to_GFA, path_to_FWM, path_to_Matdir, slices, type, plot_fODF, plot_fODFmax)

Vmask = spm_vol(path_to_mask);
Mask = spm_read_vols(Vmask);
Mask = single(Mask);

try
Vgfa = spm_vol(path_to_GFA);
f_gfa = spm_read_vols(Vgfa);
f_gfa = single(f_gfa);
[nx, ny, nz] = size(f_gfa);
end

if strcmp(plot_fODF, 'yes')
Vodf = spm_vol(path_to_ODF);
ODFmatrix = spm_read_vols(Vodf);
ODFmatrix = single(ODFmatrix);
[nx, ny, nz, ndir] = size(ODFmatrix);
end

try
VWM = spm_vol(path_to_FWM);
f_wm = spm_read_vols(VWM);
f_wm = single(f_wm);
[nx, ny, nz] = size(f_wm);
end

if strcmp( plot_fODFmax, 'yes')
    try
    load (path_to_Matdir);
    catch err
        [pathstr, name, ext] = fileparts(path_to_Matdir);
        Matrix_Dir = spm_read_vols(spm_vol([pathstr filesep 'peaks' '.nii']));
    end
    % Vdir = spm_vol(path_to_Matdir);
    % Matrix_Dir = spm_read_vols(Vdir);
end

invert_z = 1;
if invert_z == 1
    Matrix_Dir(:,3) = -Matrix_Dir(:,3);
    Matrix_Dir(:,6) = -Matrix_Dir(:,6);
    Matrix_Dir(:,9) = -Matrix_Dir(:,9); 
end

if strcmp(plot_fODF, 'yes')
if nx == 181
    V = load('sampling_and_reconstruction_schemes/On_the_sphere/362_shell_trackvis.txt');
    % Obtaining the faces on one hemisphere, but considering the symmetry
    % about the plane: by construction in this scheme V(k,:) = -V(181+k,:) for all k=1:181;
    % The faces are used to extract the ODF maxima.
    % -------------------------------------
    ODFmatrix = permute(ODFmatrix,[2 3 4 1]);
    ODFmatrix = cat(4,ODFmatrix,ODFmatrix);
    [nx, ny, nz, ndir] = size(ODFmatrix);
elseif ndir == 362
    V = load('sampling_and_reconstruction_schemes/On_the_sphere/724_shell.txt');
    % Obtaining the faces on one hemisphere, but considering the symmetry
    % about the plane: by construction in this scheme V(k,:) = -V(362+k,:) for all k=1:362;
    % The faces are used to extract the ODF maxima.
    ODFmatrix = cat(4,ODFmatrix,ODFmatrix);
    ndir = 2*ndir;
elseif ndir == 724
    V = load('sampling_and_reconstruction_schemes/On_the_sphere/724_shell.txt');
    % Obtaining the faces on one hemisphere, but considering the symmetry
    % about the plane: by construction in this scheme V(k,:) = -V(362+k,:) for all k=1:362;
    % The faces are used to extract the ODF maxima.
    % ODFmatrix = cat(4,ODFmatrix,ODFmatrix);
    % ndir = 2*ndir;
elseif ndir == 321
    Vstruc = load('sampling_and_reconstruction_schemes/On_the_sphere/odf8.mat');
    V = Vstruc.odf_vertices';
    % Obtaining the faces on one hemisphere, but considering the symmetry
    % about the plane: by construction in this scheme V(k,:) = -V(321+k,:) for all k=1:321;
    % The faces are used to extract the ODF maxima.
    ODFmatrix = cat(4,ODFmatrix,ODFmatrix);
    ndir = 2*ndir;
end

F = convhulln(V);

TF = isnan(ODFmatrix);
ODFmatrix(TF) = 0;
end


h = figure;
hold on
set(h,'Color',[1 1 1], 'position', [1, 25, 1920, 949],'DefaultLineLineSmoothing','on');

%figure;
if strcmp( type, 'z')
    f_gfa = permute(f_gfa,[2 1 3]);
    zslice = slices;
    slice(double(f_gfa), [],[], zslice);
    shading interp;
    hold on;
    for k = 1:length(slices)
        subplot(1,length(slices),k);
        z = slices(k);
        if strcmp( plot_fODF, 'yes')
            for x = 1:nx
                for y = 1:ny
                    if Mask(x,y,z) == 1
                        ODF = ODFmatrix(x,y,z,:);
                        hold on
                        % fwm = f_wm(x,y,z);
                        plot_ODF(0.5*(ODF(:)./max(ODF)),V, F, [(x) (y) (z + 1)]);
                    end
                end
            end
        end
        set(gcf, 'color', 'white');
        axis equal;
        axis ([0, nx+1, 0, ny+1, z-2, z+2]);
        grid on;
        xlabel(['Slice: ' num2str(z)]);
        if strcmp( plot_fODFmax, 'yes')
            hold on
            plot_peaks(z, Matrix_Dir, Mask, type);
        else
            shading interp;
        end
        view(0,90);
    end
end

if strcmp( type, 'y')
    f_gfa = permute(f_gfa,[2 1 3]);
    yslice = slices;
    slice(double(f_gfa), [], yslice, []);
    %slice(double(permute(f_wm,[2 1 3])), [], yslice, []);
    % colormap gray;
    shading interp;
    hold on;
    for k = 1:length(slices)
        subplot(1,length(slices),k);
        y = slices(k);
        if strcmp( plot_fODF, 'yes')
            for x = 1:nx
                for z = 1:nz
                    if Mask(x,y,z) == 1
                        ODF = ODFmatrix(x,y,z,:);
                        hold on
                        % fwm = f_wm(x,y,z);
                        %plot_ODF(fwm*0.5*(ODF(:)./max(ODF)),V, F, [(x) (y - 1) (z)]);
                        plot_ODF(0.5*(ODF(:)./max(ODF)),V, F, [(x) (y - 1) (z)]);
                    end
                end
            end
        end
        %shading interp;
        set(gcf, 'color', 'white');
        axis equal;
        axis ([0, nx+1,  y-2, y+2, 0, nz+1]);
        grid on;
        xlabel(['Slice: ' num2str(y)]);
        % shading interp;
        if strcmp( plot_fODFmax, 'yes')
            hold on
            plot_peaks(y, Matrix_Dir, Mask, type);
        else
            shading interp;
        end
        view(0,0);
    end
end

if strcmp( type, 'x')
    f_gfa = permute(f_gfa,[2 1 3]);
    xslice = slices;
    slice(double(f_gfa), xslice, [], []);
    hold on;
    shading interp;
    for k = 1:length(slices)
        subplot(1,length(slices),k);
        x = slices(k);
        if strcmp( plot_fODF, 'yes')
            for y = 1:ny
                for z = 1:nz
                    if Mask(x,y,z) == 1
                        ODF = ODFmatrix(xslice,y,z,:);
                        hold on
                        fwm = f_wm(x,y,z);
                        plot_ODF(fwm*0.5*(ODF(:)./max(ODF)),V, F, [(x + 1) (y) (z)]);
                    end
                end
            end
        end
        set(gcf, 'color', 'white');
        axis equal;
        axis ([xslice-2, xslice+2, 0, ny+1, 0, nz+1]);
        grid on;
        xlabel(['Slice: ' num2str(xslice)]);
        if strcmp( plot_fODFmax, 'yes')
            hold on
            plot_peaks(xslice, Matrix_Dir, Mask, type);
        else
            shading interp;
        end
        view (90,0)
    end
end

%camlight left
%camlight rigth
% camlight head
% material metal
% lighting  gouraud

%shading interp;
colormap gray;
camlight

% set(gcf,'PaperPositionMode', 'auto');
% set(gcf,'units','pixels');
% set(gca,'units','pixels');
% w_pos = get(gcf, 'position');
% set(gca, 'position', [10 10 w_pos(3)-10 w_pos(4)-10]);

end

function [] = plot_ODF(odf,gradn,F, Origen)
if nargin ==3
    Origen = [0 0 0];
end

%- plot_figures -%
% colordef black;
colordef white;
plotVSD(odf,gradn,F,1,Origen);
end

function plotVSD(S,V,F,scale,Origin)
D = S;
DV = repmat(D,[1 3]) .* V;
patchsignal(DV, F, Origin, scale);
end

function patchsignal(V, F, Origin, scale)

OV = (scale * V) + repmat(Origin,[size(V,1) 1]);
a = V ./ ...
    (repmat(sqrt(dot(V,V,2)),1,3) + eps);

patch('Vertices', OV, 'Faces', F, ...
    'FaceVertexCData', abs(a),...
    'FaceColor', 'interp',...
    'FaceLighting','phong',...
    'LineStyle','-',...
    'LineWidth',.05,...
    'EdgeColor',[.3 .3 .3],...
    'AmbientStrength',.4,...
    'FaceLighting','phong',...
    'SpecularColorReflectance',.2,...
    'DiffuseStrength',.5,...
    'BackFaceLighting', 'reverselit');

% patch('vertices', V, 'faces', F, ...
%       'facecolor', [0.7 0.8 1], 'edgecolor', [.2 .2 .6]);
end

function [] = plot_peaks(z, Matrix_Dir, Mask, type)
% --- Plot fiber directions ------
LineWidth = 1;
Rad = 0.04; Nf = 16;
[n1,n2,n3, kkk] = size( Matrix_Dir );
if strcmp( type, 'z')
    for i=1:n1
        for j=1:n2
            for h=z
                num_fib = length(Matrix_Dir(i,j,h,:))/3;
                for k=1:num_fib
                    if Mask(i,j,h) == 1
                        nf = 3*k;
                        nf_ind = nf-2:nf;
                        v1 = Matrix_Dir(i,j,h,nf_ind);
                        f1 = sqrt(sum(v1.^2));
                        if sum(f1 > 0)
                            if f1 > 0.90
                                f1 = 0.6;
                            else
                                f1 = 1.2;
                            end
                            iv = i; jv = j; hv = h;
                            hz = hv + 1;
                            line([iv-v1(1)*f1, iv+v1(1)*f1],[jv-v1(2)*f1, jv+v1(2)*f1],[hz-v1(3)*f1, hz+v1(3)*f1],'Color', abs([v1(1) v1(2) v1(3)]),'LineWidth',LineWidth);
                        end
                    end
                end
            end
        end
    end
end
if strcmp( type, 'y')
    for i=1:n1
        for j=z
            for h=1:n3
                num_fib = length(Matrix_Dir(i,j,h,:))/3;
                if Mask(i,j,h) == 1
                    for k=1:num_fib
                        nf = 3*k;
                        nf_ind = nf-2:nf;
                        v1 = Matrix_Dir(i,j,h,nf_ind);
                        f1 = sqrt(sum(v1.^2));
                        v1 = v1/f1;
                        if f1 > 0
                            %f1 = 1.2*f1;
                            if f1 >= 0.80
                                f1 = 0.5*f1;
                            else
                                f1 = 0.1*f1;
                            end
                            
                            %elseif ( (f1 >= 1/3) && (f1 < 0.70) )
                             %   f1 = 2*0.5*f1;
                            %else
                              %  f1 = 3*0.5*f1;
                            %end
                            
                            iv = i; jv = j; hv = h;
                            jz = jv-1;
                            
%                             %% Cylinders                            
%                             r1 = [iv-v1(1)*f1, jz-v1(2)*f1, hv-v1(3)*f1];
%                             r2 = [iv+v1(1)*f1, jz+v1(2)*f1, hv+v1(3)*f1];
%                             RGBcolor = max(0, min(1, abs((r2-r1)/norm(r2-r1))));
%                             
%                             [XX YY ZZ] = plot_cylinder2P(Rad, Nf, r1, r2);
%                             [F, V] = surf2patch(XX,YY,ZZ);
%                             patch('vertices', V, 'faces', F, 'facecolor', RGBcolor, 'edgecolor', 'none');
%                             
%                             %  fill the cylinders
%                             patch(XX(1,:)  ,YY(1,:)  ,ZZ(1,:)  , ZZ(1,:),   'facecolor', RGBcolor);
%                             patch(XX(end,:),YY(end,:),ZZ(end,:), ZZ(end,:), 'facecolor', RGBcolor);
%                             
%                             %Color_i = [v1(1) v1(2) v1(3)]/sqrt( v1(1)^2 + v1(2)^2  + v1(3)^2 );
%                             %Color_i = Color_i.*sign(Color_i);
%                             %hhh = streamtube({[iv-v1(1)*f1 jz-v1(2)*f1 hv-v1(3)*f1; iv+v1(1)*f1 jz+v1(2)*f1 hv+v1(3)*f1]}, 0.1);
%                             %set(hhh,'EdgeColor','none','FaceColor',Color_i, 'CDataMapping','direct','AmbientStrength', 0.3) %0.3
                            
                            % Line
                            
                            line([iv-v1(1)*f1, iv+v1(1)*f1],[jz-v1(2)*f1, jz+v1(2)*f1],[hv-v1(3)*f1, hv+v1(3)*f1],'Color', abs([v1(1) v1(2) v1(3)]),'LineWidth',LineWidth);
                        end
                    end
                end
            end
        end
    end
end
if strcmp( type, 'x')
    for i=z
        for j=1:n2
            for h=1:n3
                num_fib = length(Matrix_Dir(i,j,h,:))/3;
                for k=1:num_fib
                    if Mask(i,j,h) == 1
                        nf = 3*k;
                        nf_ind = nf-2:nf;
                        v1 = Matrix_Dir(i,j,h,nf_ind);
                        f1 = sqrt(sum(v1.^2));
                        if sum(f1 > 0)
                            if f1 > 0.90
                                f1 = 0.6;
                            else
                                f1 = 1.2;
                            end
                            iv = i; jv = j; hv = h;
                            iz = iv + 1;
                            line([iz-v1(1)*f1, iz+v1(1)*f1],[jv-v1(2)*f1, jv+v1(2)*f1],[hv-v1(3)*f1, hv+v1(3)*f1],'Color', abs([v1(1) v1(2) v1(3)]),'LineWidth',LineWidth);
                        end
                    end
                end
            end
        end
    end
end
end

function [X Y Z] = plot_cylinder2P(R, N, r1, r2)

    % The parametric surface will consist of a series of N-sided
    % polygons with successive radii given by the array R.
    % Z increases in equal sized steps from 0 to 1.
    
    % r1=[x1 y1 z1];r2=[x2 y2 z2]; and R, the radius could be a parameterized using a vector R(i) 
    % instead of a scalar (as for cylinder). 

    % Set up an array of angles for the polygon.
    theta = linspace(0,2*pi,N);

    m = length(R);                 % Number of radius values
                                   % supplied.

    if m == 1                      % Only one radius value supplied.
        R = [R; R];                % Add a duplicate radius to make
        m = 2;                     % a cylinder.
    end


    X = zeros(m, N);             % Preallocate memory.
    Y = zeros(m, N);
    Z = zeros(m, N);
    
    v=(r2-r1)/sqrt((r2-r1)*(r2-r1)');    %Normalized vector;
    %cylinder axis described by: r(t)=r1+v*t for 0<t<1
    R2=rand(1,3);              %linear independent vector (of v)
    x2=v-R2/(R2*v');    %orthogonal vector to v
    x2=x2/sqrt(x2*x2');     %orthonormal vector to v
    x3=cross(v,x2);     %vector orthonormal to v and x2
    x3=x3/sqrt(x3*x3');
    
    r1x=r1(1);r1y=r1(2);r1z=r1(3);
    r2x=r2(1);r2y=r2(2);r2z=r2(3);
    vx=v(1);vy=v(2);vz=v(3);
    x2x=x2(1);x2y=x2(2);x2z=x2(3);
    x3x=x3(1);x3y=x3(2);x3z=x3(3);
    
    time=linspace(0,1,m);
    for j = 1 : m
      t=time(j);
      X(j, :) = r1x+(r2x-r1x)*t+R(j)*cos(theta)*x2x+R(j)*sin(theta)*x3x; 
      Y(j, :) = r1y+(r2y-r1y)*t+R(j)*cos(theta)*x2y+R(j)*sin(theta)*x3y; 
      Z(j, :) = r1z+(r2z-r1z)*t+R(j)*cos(theta)*x2z+R(j)*sin(theta)*x3z;
    end

    % surf(X, Y, Z);
    % The parametric surface will consist of a series of N-sided
    % polygons with successive radii given by the array R.
    % Z increases in equal sized steps from 0 to 1.

    % Set up an array of angles for the polygon.
    theta = linspace(0,2*pi,N);

    m = length(R);                 % Number of radius values
                                   % supplied.

    if m == 1                      % Only one radius value supplied.
        R = [R; R];                % Add a duplicate radius to make
        m = 2;                     % a cylinder.
    end


    X = zeros(m, N);             % Preallocate memory.
    Y = zeros(m, N);
    Z = zeros(m, N);
    
    v=(r2-r1)/sqrt((r2-r1)*(r2-r1)');    %Normalized vector;
    %cylinder axis described by: r(t)=r1+v*t for 0<t<1
    R2=rand(1,3);              %linear independent vector (of v)
    x2=v-R2/(R2*v');    %orthogonal vector to v
    x2=x2/sqrt(x2*x2');     %orthonormal vector to v
    x3=cross(v,x2);     %vector orthonormal to v and x2
    x3=x3/sqrt(x3*x3');
    
    r1x=r1(1);r1y=r1(2);r1z=r1(3);
    r2x=r2(1);r2y=r2(2);r2z=r2(3);
    vx=v(1);vy=v(2);vz=v(3);
    x2x=x2(1);x2y=x2(2);x2z=x2(3);
    x3x=x3(1);x3y=x3(2);x3z=x3(3);
    
    time=linspace(0,1,m);
    for j = 1 : m
      t=time(j);
      X(j, :) = r1x+(r2x-r1x)*t+R(j)*cos(theta)*x2x+R(j)*sin(theta)*x3x; 
      Y(j, :) = r1y+(r2y-r1y)*t+R(j)*cos(theta)*x2y+R(j)*sin(theta)*x3y; 
      Z(j, :) = r1z+(r2z-r1z)*t+R(j)*cos(theta)*x2z+R(j)*sin(theta)*x3z;
    end

    % surf(X, Y, Z,'FaceColor','red','EdgeColor','none');
end
