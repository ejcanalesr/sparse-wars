% Erick J. Canales-Rodríguez, Jon Haitz Legarreta-Gorroño, Yasser Alemán-Gomez, Marco Pizzolato, Jonathan Rafael Patiño, Gaëtan Olivier D Rensonnet, Muhamed Barakovic, Joaquim Radua, Edith Pomarol-Clotet, Raymond Salvador, Jean-Philippe Thiran, Alessandro Daducci. Sparse wars: A survey and comparative study of spherical deconvolution algorithms for diffusion MRI. Neuroimage. 2019 Jan 1;184:140-160.                  

function [S, D] = create_signal_multi_tensor (ang, f, Eigenvalues, b, grad, S0, noise_param)
% -Normalizing the gradient vector and then transforming the b-value.
% -This part is only for non-unitary gradient vectors
% Transform_b_value_and_grad = repmat(sqrt(diag(grad*grad')+eps), [1 3]);
% grad = grad./Transform_b_value_and_grad;
% b = b.*(Transform_b_value_and_grad(:,1)).^2;

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

if noise_param.add_rician_noise == 1
    SNR = noise_param.SNR;
    sigma = S0/SNR;
    standar_deviation = sigma.*(ones(length(grad),1));
    med = zeros(length(grad),1);
    
    er1 = normrnd(med, standar_deviation);
    er2 = normrnd(med, standar_deviation);
    
    S = sqrt((S + er1).^2 + er2.^2);
end

return

% --- private funtions -----------
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
return
