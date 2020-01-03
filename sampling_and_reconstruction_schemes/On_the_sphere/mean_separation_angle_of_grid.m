% CosA = V*V';
% CosA(CosA>1)=1;
% CosA(CosA<-1)=-1;
% angles = 180*acos(CosA);
% size(angles)
load ODF_XYZ;
V = ODF_XYZ;
Vfaces = convhulln(V);

% --- parameters to find_discrete_peaks
cos_ang = V*V';
% --- correction due to numerical errors
cos_ang(cos_ang >  1)  =  1;
cos_ang(cos_ang < -1)  = -1;
% --------------------------------------
ang = acosd(cos_ang);
% --- considering the antipodal symmetry
ang = min(ang,180-ang);
% --------------------------------------
angles_all_pairs = [];
for i=1:length(Vfaces);
    angles_all_pairs = [angles_all_pairs; ang(Vfaces(i,1),Vfaces(i,2)); ang(Vfaces(i,1),Vfaces(i,3)); ang(Vfaces(i,2),Vfaces(i,3))];
end
mean(angles_all_pairs)
std(angles_all_pairs)