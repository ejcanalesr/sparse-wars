% Erick J. Canales-Rodríguez, Jon Haitz Legarreta-Gorroño, Yasser Alemán-Gomez, Marco Pizzolato, Jonathan Rafael Patiño, Gaëtan Olivier D Rensonnet, Muhamed Barakovic, Joaquim Radua, Edith Pomarol-Clotet, Raymond Salvador, Jean-Philippe Thiran, Alessandro Daducci. Sparse wars: A survey and comparative study of spherical deconvolution algorithms for diffusion MRI. Neuroimage. 2019 Jan 1;184:140-160.

function [] = plot_max(path_to_mask, path_to_peaks, slices, type, invert_z, plot_axis_limit)

try
    Mask = load(path_to_mask);
    Mask = Mask.flag3;
    Mask(Mask>0) = 1;
end

try
    Vol = spm_vol(path_to_mask);
    Mask = spm_read_vols(Vol);
end

try
    Matrix_Dir = load (path_to_peaks);
    Matrix_Dir = Matrix_Dir.Matrix_Dir;
catch err
    Matrix_Dir = spm_read_vols(spm_vol(path_to_peaks));
end

if invert_z == 1
    Matrix_Dir(:,:,:,3) = -Matrix_Dir(:,:,:,3);
    Matrix_Dir(:,:,:,6) = -Matrix_Dir(:,:,:,6);
    Matrix_Dir(:,:,:,9) = -Matrix_Dir(:,:,:,9);
end


[nx, ny, nz, ~] = size(Matrix_Dir);

iso_f1 = sqrt(sum(Matrix_Dir(:,:,:,1:3).^2,4));
iso_f2 = sqrt(sum(Matrix_Dir(:,:,:,4:6).^2,4));
iso_f3 = sqrt(sum(Matrix_Dir(:,:,:,7:9).^2,4));

%iso_f = iso_f1 + iso_f2 + iso_f3 + iso_f4;
iso_f = Mask;

h = figure;
hold on
set(h,'Color',[1 1 1], 'position', [1, 25, 1920, 949], 'DefaultLineLineSmoothing','on');

if strcmp( type, 'z')
    % iso_f = permute(iso_f,[1 2 3]);
    zslice = slices;
    slice(double(iso_f), [],[], zslice);
    shading interp;
    hold on;
    for k = 1:length(slices)
        subplot(1,length(slices),k);
        z = slices(k);
        set(gcf, 'color', 'white');
        axis equal;
        axis ([0, nx+1, 0, ny+1, z-4, z+4]);
        grid on;
        xlabel(['Slice: ' num2str(z)]);
        hold on
        plot_peaks(z, Matrix_Dir, Mask, type);
        view(0,90);
    end
end

if strcmp( type, 'y')
    iso_f = permute(iso_f,[2 1 3]);
    yslice = slices;
    slice(double(iso_f), [], yslice, []);
    shading interp;
    hold on;
    x_min = round(plot_axis_limit(1));
    x_max = round(plot_axis_limit(2));
    y_min = round(plot_axis_limit(5));
    y_max = round(plot_axis_limit(6));

    for k = 1:length(slices)
        subplot(1,length(slices),k);
        y = slices(k);
        set(gcf, 'color', 'white');
        axis equal;
        axis ([0, nx+1,  y-2, y+2, 0, nz+1]);
        grid on;
        xlabel(['Slice: ' num2str(y)]);
        hold on
        plot_peaks(y, Matrix_Dir, Mask, type, x_min, x_max, y_min, y_max);
        view(0,0);
    end
    axis(plot_axis_limit);
end

if strcmp( type, 'x')
    iso_f = permute(iso_f,[2 1 3]);
    xslice = slices;
    slice(double(iso_f), xslice, [], []);
    shading interp;
    hold on;
    for k = 1:length(slices)
        subplot(1,length(slices),k);
        x = slices(k);
        set(gcf, 'color', 'white');
        axis equal;
        axis ([x-2, x+2, 0, ny+1, 0, nz+1]);
        grid on;
        xlabel(['Slice: ' num2str(x)]);
        hold on
        plot_peaks(x, Matrix_Dir, Mask, type);
        view (90,0)
    end
end

colormap gray;
camlight
set(gca,'DefaultLineLineSmoothing','on');
caxis([0,0.5]);

box on; grid off;
set(gcf,'PaperPositionMode', 'auto');
set(gcf,'units','pixels');
set(gca,'units','pixels');
w_pos = get(gcf, 'position');
set(gca, 'position', [10 10 w_pos(3)-10 w_pos(4)-10]);

end

function [] = plot_peaks(z, Matrix_Dir, Mask, type, x_min, x_max, y_min, y_max)
% --- Plot fiber directions ------
[n1,n2,n3, ~] = size( Matrix_Dir );
Rad = 0.08; Nf = 32;
if strcmp( type, 'z')
    for i = x_min:x_max
        for j = y_min:y_max
            for h=z
                num_fib = length(Matrix_Dir(i,j,h,:))/3;
                for k=1:num_fib
                    if Mask(i,j,h) > 0.1
                        nf = 3*k;
                        nf_ind = nf-2:nf;
                        v1 = Matrix_Dir(i,j,h,nf_ind);
                        if sum(abs(v1) > 0.1)
                            f1 = 0.5;
                            iv = i; jv = j; hv = h;
                            hz = hv + 2;
                            %% Line
                            % line([iv-v1(1)*f1, iv+v1(1)*f1],[jv-v1(2)*f1, jv+v1(2)*f1],[hz-v1(3)*f1, hz+v1(3)*f1],'Color', abs([v1(1) v1(2) v1(3)]),'LineWidth',3);                        
                            
                            %% Cylinders                            
                            % r1 = [iv-v1(1)*f1, jv-v1(2)*f1, hz-v1(3)*f1];
                            % r2 = [iv+v1(1)*f1, jv+v1(2)*f1, hz+v1(3)*f1];
                            
                            r1 = [jv iv hz] - [v1(1) v1(2) v1(3)];
                            r2 = [jv iv hz] + [v1(1) v1(2) v1(3)];

                            RGBcolor = abs((r2-r1)/norm(r2-r1));
                            
                            [XX, YY, ZZ] = plot_cylinder2P(Rad, Nf, r1, r2);
                            [F, V] = surf2patch(XX,YY,ZZ);
                            patch('vertices', V, 'faces', F, 'facecolor', RGBcolor, 'edgecolor', 'none');
                            
                            %  fill the cylinders
                            patch(XX(1,:)  ,YY(1,:)  ,ZZ(1,:)  , ZZ(1,:),   'facecolor', RGBcolor);
                            patch(XX(end,:),YY(end,:),ZZ(end,:), ZZ(end,:), 'facecolor', RGBcolor);
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
                for k=1:num_fib
                    if Mask(i,j,h) > 0.1
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
                            jz = jv-1;
                            %% Line
                            %line([iv-v1(1)*f1, iv+v1(1)*f1],[jz-v1(2)*f1, jz+v1(2)*f1],[hv-v1(3)*f1, hv+v1(3)*f1],'Color', abs([v1(1) v1(2) v1(3)]),'LineWidth',3);
                            
                            %% Cylinders                            
                            r1 = [iv-v1(1)*f1, jz-v1(2)*f1, hv-v1(3)*f1];
                            r2 = [iv+v1(1)*f1, jz+v1(2)*f1, hv+v1(3)*f1];
                            RGBcolor = abs((r2-r1)/norm(r2-r1));
                            
                            [XX, YY, ZZ] = plot_cylinder2P(Rad, Nf, r1, r2);
                            [F, V] = surf2patch(XX,YY,ZZ);
                            patch('vertices', V, 'faces', F, 'facecolor', RGBcolor, 'edgecolor', 'none');
                            
                            %  fill the cylinders
                            patch(XX(1,:)  ,YY(1,:)  ,ZZ(1,:)  , ZZ(1,:),   'facecolor', RGBcolor);
                            patch(XX(end,:),YY(end,:),ZZ(end,:), ZZ(end,:), 'facecolor', RGBcolor);
                        end
                    end
                end
            end
        end
    end
end
if strcmp( type, 'x')
    for i=x
        for j=1:n2
            for h=1:n3
                num_fib = length(Matrix_Dir(i,j,h,:))/3;
                for k=1:num_fib
                    if Mask(i,j,h) > 0.1
                        nf = 3*k;
                        nf_ind = nf-2:nf;
                        v1 = Matrix_Dir(i,j,h,nf_ind);
                        if sum(abs(v1) > 0)
                            f1 = 1;
                            iv = i; jv = j; hv = h;
                            iz = iv + 1;
                            %% Line
                            line([iz-v1(1)*f1, iz+v1(1)*f1],[jv-v1(2)*f1, jv+v1(2)*f1],[hv-v1(3)*f1, hv+v1(3)*f1],'Color', abs([v1(1) v1(2) v1(3)]),'LineWidth',3);
                            
                            %% Cylinders                            
                            r1 = [iz-v1(1)*f1, jv-v1(2)*f1, hv-v1(3)*f1];
                            r2 = [iz+v1(1)*f1, jv+v1(2)*f1, hv+v1(3)*f1];
                            RGBcolor = abs((r2-r1)/norm(r2-r1));
                            
                            [XX, YY, ZZ] = plot_cylinder2P(Rad, Nf, r1, r2);
                            [F, V] = surf2patch(XX,YY,ZZ);
                            patch('vertices', V, 'faces', F, 'facecolor', RGBcolor, 'edgecolor', 'none');
                            
                            %  fill the cylinders
                            patch(XX(1,:)  ,YY(1,:)  ,ZZ(1,:)  , ZZ(1,:),   'facecolor', RGBcolor);
                            patch(XX(end,:),YY(end,:),ZZ(end,:), ZZ(end,:), 'facecolor', RGBcolor); 
                        end
                    end
                end
            end
        end
    end
end
end


function [X, Y, Z] = plot_cylinder2P(R, N, r1, r2)

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
