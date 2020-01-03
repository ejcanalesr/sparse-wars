% Erick J. Canales-Rodríguez, Jon Haitz Legarreta-Gorroño, Yasser Alemán-Gomez, Marco Pizzolato, Jonathan Rafael Patiño, Gaëtan Olivier D Rensonnet, Muhamed Barakovic, Joaquim Radua, Edith Pomarol-Clotet, Raymond Salvador, Jean-Philippe Thiran, Alessandro Daducci. Sparse wars: A survey and comparative study of spherical deconvolution algorithms for diffusion MRI. Neuroimage. 2019 Jan 1;184:140-160.

function [peaks_dir, peaks_f, peaks] = find_discrete_peaks(V, ODF, neig, thr)

%% --- discrete maxima ---
peaks = zeros(size(ODF));
thr_value = thr*max(ODF);
for i = 1:length(ODF)
    neig_i = neig(i,:);
    %ind_valid_i =  (neig_i ~= i) & (neig_i > 0) ;
    ind_valid_i = neig_i > 0 ;
    neig_final = neig_i(ind_valid_i);
    if ( ODF(i) > max(ODF(neig_final)) ) && ( ODF(i) >= thr_value )
        peaks(i) = 1;
    end
end
peaks_dir = V(find(peaks),:);
peaks_f = ODF(find(peaks));
end
