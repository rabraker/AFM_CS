clc
clear
pix = 256;
width = 5
pix_per_volt = (pix/width)*5

dat = csvread('C:\Users\arnold\Documents\labview\afm_imaging/data/cs-data-out02.csv');

% Drop all data corresponding to movement between points.
ind_meas = find(dat(:, 6) ~= 0);  % Index of the movements are all 0.
dat_meas = dat(ind_meas, :);

% Fit a line to the control data and subtract it. Helps remove some of
% the piezo drift. 
dat_meas(:,5) = detrend(dat_meas(:,5)); % detrend u_z as a long vector

% Split each measurement entities into a cell array.
%%
% bin all the data into pixels. 
for k = 1:max(dat_meas(:,end))
   
   % Get indexes corresponding to measurement entity k
   inds = find(dat_meas(:,end) == k); 
   
   % The x-direction has the most movement. Find the total spread and
   % approximate calculate the number of pixels that corresponds to. 
   X_ks = dat_meas(inds, 1);
   x_spread = max(X_ks) - min(X_ks);
   path_pix = ceil(x_spread*pix_per_volt);
   
   Y_ks = dat_meas(inds, 2);
   E_ks = dat_meas(inds, 3);
   U_ks = dat_meas(inds, 5);
   
   % Convert the x and y voltage data into bined and averaged pixel coordinates.
   x_pix = floor(pixelize(X_ks, npix)*pix);
   y_pix = floor(pixelize(Y_ks, npix)*pix);
   
   % Bin the height/deflection data into pixels and average. 
   e_pix = pixelize(E_ks, path_pix);
   u_pix = pixelize(U_ks, path_pix);
   
   % Finally, collect in the kth cell array entry. 
  pix_dat{k} = [x_pix, y_pix, e_pix, u_pix];
end
%%
close all
figure; hold on;
hole_depth = (20/7)*(1/1000)*(20)
I = zeros(pix, pix);
pixelifsampled = I;
for k = 1:length(pix_dat)
    u_pix = pix_dat{k}(:,4);
    x_pix = pix_dat{k}(:,1);
    y_pix = pix_dat{k}(:,2);

    if abs(max(u_pix) - min(u_pix))> hole_depth*.5
       % Then we have an edge. 
        have_edge = 1;
        cs = 'r';
    else
        have_edge = 0;
        cs = 'b';
    end
    
    % The mean doesn't really give us what we want. When a path contains
    % a hold and the flat area, the mean splits the difference and the
    % flat area data is "too high".
    %         u_pix = u_pix - mean(u_pix);

    % We assume that each path contains data from the flat surface.
    % Subtracting off the max is a crude way of registering all the path
    % data to the same height. 
    u_pix = u_pix - max(u_pix);
    
    % Just to get a feel for things:
    plot(u_pix, 'color', cs)
    ylim([-hole_depth, hole_depth])
    
    for jj = 1:length(u_pix)
        % Remember how we converted the x and y data into pixel
        % coordinates? Those now act as the indices to fill in the image
        % data.
            n_row = y_pix(jj)+1;
            m_col = x_pix(jj)+1;
            I(n_row, m_col) = u_pix(jj);
            pixelifsampled(n_row, m_col) = 1;
    end
    
    
end
figure

lo = min(min(I));
hi = max(max(I));
imshow(flipud(I), [lo, hi])


%%

addpath('C:\Users\arnold\Documents\labview\reconstruction\BP\')
addpath('C:\Users\arnold\Documents\labview\reconstruction\BP\l1magic');
[n m] = size(I);

I_vector = PixelMatrixToVector(I);
pixelifsampled = PixelMatrixToVector(pixelifsampled);
I_vector = I_vector(find(pixelifsampled>0.5));

A = @(x) IDCTfun(x,pixelifsampled);
At = @(x) DCTfun(x,pixelifsampled);


close all;


Ir = idct(l1qc_logbarrier(At(I_vector), A, At, I_vector, 0.1));


close all;


Ir = real(Ir);
subplot(1,2,1)
imshow(I.*PixelVectorToMatrix(pixelifsampled,[n m]),[lo, hi]);
title('sample');
subplot(1,2,2)
imshow(PixelVectorToMatrix(Ir,[n m]),[min(Ir), max(Ir)]);
title('BP reconstruction');







