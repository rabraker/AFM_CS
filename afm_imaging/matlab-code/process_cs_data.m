clc
clear
pix = 256;
width = 5;
microns_per_volt = 50/10;
pix_per_volt = (pix/width)*microns_per_volt;

Ts = 40e-6;

% This path should be set to the root of your github working directory. 
% root = 'C:\Users\arnold\Documents\labview';  % For windows
root = '/home/arnold/gradschool/afm-cs';
addpath(fullfile(root, 'afm_imaging/matlab-code/functions'));

dat = csvread(fullfile(root, 'afm_imaging/data/cs-data-out02.csv'));

% Drop all data corresponding to movement between points.
ind_meas = find(dat(:, 6) ~= 0);  % Index of the movements are all 0.
dat_meas = dat(ind_meas, :);

% Fit a line to the control data and subtract it. Helps remove some of
% the piezo drift. 
dat_meas(:,5) = detrend(dat_meas(:,5)); % detrend u_z as a long vector






% Split each measurement entities into a cell array.

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
   x_pix = floor(pixelize(X_ks, path_pix)*pix);
   y_pix = floor(pixelize(Y_ks, path_pix)*pix);
   
   
%    keyboard
   % Bin the height/deflection data into pixels and average. 
   e_pix = pixelize(E_ks, path_pix);
   u_pix = pixelize(U_ks, path_pix);
   
   % Finally, collect in the kth cell array entry. 
  pix_dat{k} = [x_pix, y_pix, e_pix, u_pix];
end
%
% close all
f1 = figure(1);
ylabel('u_z')
xlabel('pixel')
ax1 = gca();
hold on;

hole_depth = (20/7)*(1/1000)*(20)
I = zeros(pix, pix);
pixelifsampled = I;
for k = 1:length(pix_dat)
    u_pix = pix_dat{k}(:,4);
    x_pix = pix_dat{k}(:,1);
    y_pix = pix_dat{k}(:,2);

    % Hueristic to determine if a path contained a hole. At the moment,
    % just use this information to plot the uz data for each path. 
    if abs(max(u_pix) - min(u_pix))> hole_depth*.5
       % Then we have an edge. 
       % have_edge = 1;
        cs = 'r';
    else
        % have_edge = 0;
        cs = 'b';
    end
    
    % The mean doesn't really give us what we want. When a path contains
    % a hole and the flat area, the mean splits the difference and the
    % flat area data is "too high".
    %         u_pix = u_pix - mean(u_pix);

    % We assume that each path contains data from the flat surface.
    % Subtracting off the max is a crude way of registering all the path
    % data to the same height. 
    u_pix = u_pix - max(u_pix);
    
    % Just to get a feel for things:
    plot(ax1, u_pix, 'color', cs)
    ylim([-hole_depth, hole_depth])
    
    for jj = 1:length(u_pix)
        % Remember how we converted the x and y data into pixel
        % coordinates? Those now act as the indices to fill in the image
        % data.
            n_row = y_pix(jj)+1;  % pixels start at 0. matlab is 1-based indexing.
            m_col = x_pix(jj)+1;
            I(n_row, m_col) = u_pix(jj);
            pixelifsampled(n_row, m_col) = 1;
    end
    

end
%%
I = detrend_plane(I, pixelifsampled);
% f2 = figure(2); clf
% ax2 = gca();
% 
% imshow_sane(I, ax2, width, width);
% 

%%
% Make the image square, to use smp.
I = I(1:end-1, 19:end);
size(I)
pixelifsampled = pixelifsampled(1:end-1, 19:end);
size(pixelifsampled)

% ********* SMP *************
reconstruct_root = fullfile(root, 'reconstruction/SMP');
reconstruct_path = genpath(reconstruct_root)
addpath(reconstruct_path)

[n,m] = size(I);
Ir_smp = SMP(I.*pixelifsampled,pixelifsampled,round(0.02*n*m));


reconstruct_root = fullfile(root, 'reconstruction/BP');
reconstruct_path = genpath(reconstruct_root)
addpath(reconstruct_path)

% ********* BP *************
[n m] = size(I);

I_vector = PixelMatrixToVector(I);

pixelifsampled = PixelMatrixToVector(pixelifsampled);
I_vector = I_vector(find(pixelifsampled>0.5));

A = @(x) IDCTfun(x,pixelifsampled);
At = @(x) DCTfun(x,pixelifsampled);

Ir_bp = idct(l1qc_logbarrier(At(I_vector), A, At, I_vector, 0.1));
Ir_bp = real(Ir_bp);

% close all;
f5 = figure(5)
subplot(1,3,1)
ax3 = gca();
imshow_sane(I.*PixelVectorToMatrix(pixelifsampled,[n m]), ax3, width, width);
title('sample');

subplot(1,3,2)
ax4 = gca();
imshow_sane(PixelVectorToMatrix(Ir_bp,[n m]), ax4, width, width);
title('BP reconstruction');

subplot(1,3,3)
ax5 = gca();
imshow_sane(Ir_smp, ax5, width, width)
title('SMP reconstruction');





