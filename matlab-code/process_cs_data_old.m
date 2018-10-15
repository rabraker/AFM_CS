clc
clear


% This path should be set to the root of your github working directory. 
root = 'C:\Users\arnold\Documents\labview';  % For windows
% root = '/home/arnold/gradschool/afm-cs';
addpath(fullfile(root, 'afm_imaging/matlab-code/functions'));

Ts = 40e-6;




cs_exp_data_name ='cs-traj-512pix-8perc-500nm-5mic-01Hz_out_9-18-2017-02.csv';
cs_exp_data_name_s{1} = 'cs-traj-512pix-10perc-500nm-5mic-01Hz_out_10-14-2018-02.csv';
cs_exp_data_name_s{2} = 'cs-traj-512pix-10perc-500nm-5mic-01Hz_out_10-14-2018-03.csv';
sub_dir = '5microns/10-14-2018'



for i=1:length(cs_exp_data_name_s)
    cs_exp_data_name = cs_exp_data_name_s{i};

data_root = fullfile(PATHS.exp(), 'imaging', 'cs-imaging', sub_dir);

cs_exp_meta_name = strrep(cs_exp_data_name, '.csv', '-meta.mat');

cs_data_path = fullfile(data_root, cs_exp_data_name);
cs_meta_path = fullfile(data_root, cs_exp_meta_name);

dat = csvread(cs_data_path); 
load(cs_meta_path);  % Provides ExpMetaData

state_ticks = ExpMetaData.state_counts;
state_times = state_ticks*Ts
time_total = sum(state_times)

k = regexp(cs_data_path, 'Hz_out');
meta_in_path = sprintf('%s.mat', cs_data_path(1:k+1));
if exist(meta_in_path, 'file') ==2
    load(meta_in_path)
    width = CsExpMetaIn.width;
    npix =  CsExpMetaIn.npix;
else
    npix = 256;
    width = 10;
end
microns_per_volt = 50/10;
pix_per_volt = (npix/width)*microns_per_volt;
hole_depth = (20/7)*(1/1000)*(2*120);

%
% Drop all data corresponding to movement between points.
ind_meas = find(dat(:, 6) > 0);  % Index of the movements are all 0.
meta_ind = dat(:,6);

pix = 256;
width = 5;
microns_per_volt = 50/10;
pix_per_volt = (pix/width)*microns_per_volt;

Ts = 40e-6;

% % % % This path should be set to the root of your github working directory. 
% % % % root = 'C:\Users\arnold\Documents\labview';  % For windows
% % % root = '/home/arnold/gradschool/afm-cs';
% % % addpath(fullfile(root, 'afm_imaging/matlab-code/functions'));
% % % 
% % % dat = csvread(fullfile(root, 'afm_imaging/data/cs-data-out02.csv'));
% % % 
% % % % Drop all data corresponding to movement between points.
% % % ind_meas = find(dat(:, 6) ~= 0);  % Index of the movements are all 0.


% Fit a line to the control data and subtract it. Helps remove some of
% the piezo drift. 

uz = dat_meas(:,5);
figure(1+10*i); clf;
plot(uz); hold on
% uz2 = detrend2(uz);


dat_meas(:,5) = detrend(dat_meas(:,5)); % detrend u_z as a long vector
plot(dat_meas(:,5))
%

f2 = figure(2+10*i); clf
ylabel('u_z')
xlabel('pixel')
ax2 = gca();
hold on;

% Since we want to convert xy-coords to pixels, move all data to the
% positive orthant:
dat_meas(:,1) = dat_meas(:,1) - min(dat_meas(:,1));  %x dir
dat_meas(:,2) = dat_meas(:,2) - min(dat_meas(:,2));  %y dir

I = zeros(npix,npix);

pixelifsampled=zeros(npix, npix);
% bin all the data into pixels. 
y_pixs = [];
for k = 1:max(dat_meas(:,end))

   % Get indexes corresponding to measurement entity k
   inds = find(dat_meas(:,end) == k); 
   
   % Slice out the data for the current mu-path.
   U_ks = dat_meas(inds, 4);
   Y_ks = dat_meas(inds, 2)*pix_per_volt;
   X_ks = dat_meas(inds, 1)*pix_per_volt;

   % Register the control data to zero. We can do this because we are
   % scanning long enough that we are always guaranteed to exit a hole.
   % This lets 
%    figure(f3);clf
%    plot(U_ks)
%    plot(detrend(U_ks))
%    keyboard
%    U_ks = U_ks - mean(U_ks);
if abs(min(U_ks)-max(U_ks)) < .5*hole_depth
%    U_ks = detrend(U_ks);
end
   U_ks = U_ks - max(U_ks);
   
%     if min(U_ks) < -0.12
% %         continue
%     end
    if min(U_ks) < -0.6
        continue
    end   
% U_ks = max(U_ks, -.11);
% U_ks = max(U_ks, -0.5);    

    
   % Make the assumption that the y-data for each path is constant enough.
   % Since we start at the (0,0) corner the xplane, we'll take the floor,
   % and add 1 for 1-based indexing.
   y_pix = min(npix-1, floor(mean(Y_ks)));
   if y_pix >npix-1
       keyboard
   end
   y_pixs = [y_pixs; y_pix];
   x_spread = max(X_ks) - min(X_ks);
   xpix_start = floor(min(X_ks));
   npix_path_k = ceil(x_spread); % number of bins for this path.
   
   % Now, define a set of bins for the x-direction. Each bin will have a
   % different number of data points. 
   % If we include 0, then there are three points for two pixel bins.
   xbins = linspace(min(X_ks), max(X_ks), npix_path_k+1); 
   U_k = [];
   for jj = 1:npix_path_k
       if xpix_start+jj > npix
           continue
       end
       % Get the indeces correspondiing to the current x-data bin.
       ind_x = find(X_ks >= xbins(jj) & X_ks < xbins(jj+1));
       if ~isempty(ind_x) % Avoid errors if it is empty.
           % Slice out the corresponding height data, and average it
           % together since this is a single pixel.
           u_pix_jj = mean(U_ks(ind_x)); 
           
           % The image index for the y-direction (rows) is y_pix.
           % For the x-direction, the start of the path is at xpix_start.
           % We are inching along the mu-path, so at each iteration of this
           % loop, increment the pixel index by 1, so just use loop iterator,
           % jj.
           I(y_pix+1, xpix_start+jj) = u_pix_jj;
           pixelifsampled(y_pix+1, xpix_start+jj) = 1;
           
           % Collect the uz data just for plotting/visualization
           U_k = [U_k; u_pix_jj];
           if size(I,2) ==514
               keyboard
           end
       end
       
   end
   U_k = U_k - max(U_k);
   % -------------------
   % ---- visualize ------
   if abs(max(U_k) - min(U_k))> hole_depth*.5
       % Then we have an edge. 

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
%     figure(f2)
   plot(ax2, U_k, 'color', cs)
%    keyboard
end

figure(100 + 10*i);
ax = gca();
imshow(I, [min(min(I)), max(max(I))])




% I = detrend_sampled_plane(I, pixelifsampled);
% I = (I - max(max(I))).*pixelifsampled; 
% imshow(I, [min(min(I)), max(max(I))]);

bp = 1;
% imshow_sane(I, ax2, width, width);
% Make the image square, to use smp.

size(I)
size(pixelifsampled)

% ********* SMP *************
reconstruct_root = fullfile(root, 'reconstruction', 'SMP');
reconstruct_path = genpath(reconstruct_root)
addpath(reconstruct_path)



[n,m] = size(I);

% Ir_smp = SMP(I.*pixelifsampled,pixelifsampled,round(0.02*n*m));

tic
samp_frac = sum(sum(pixelifsampled))/(n*m)
reduce_frac = 0.1;  %It's my understanding Yufan says 1/10 of sub-sample fraction.
maxiter = round(samp_frac*reduce_frac*n*m)
Ir_smp = SMP_1D(I, pixelifsampled, maxiter);
time_smp = toc;

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

f5 = figure(7 + 10*i); clf
subplot(2,3,1)
ax3 = gca();
imshow_sane(I, ax3, width, width);
title('sample');

subplot(2,3,3)
ax5 = gca();
imshow_sane(Ir_smp, ax5, width, width)
title('SMP reconstruction');
drawnow

if bp
    % ********* BP *************
    [n m] = size(I);
    tic
    I_vector = PixelMatrixToVector(I);

    pixelifsampled_vec = PixelMatrixToVector(pixelifsampled);
    I_vector = I_vector(find(pixelifsampled_vec>0.5));

    A = @(x) IDCTfun(x,pixelifsampled_vec);
    At = @(x) DCTfun(x,pixelifsampled_vec);

    Ir_bp = idct(l1qc_logbarrier(At(I_vector), A, At, I_vector, 0.1));
    Ir_bp = real(Ir_bp);
    bp_im = PixelVectorToMatrix(Ir_bp,[n m]);
    time_bp = toc;
    fprintf('SMP Time: %f \nBP Time: %f\n', time_smp, time_bp);

    bp_im = detrend_plane(bp_im);
    
    subplot(2,3,2)
    ax4 = gca();
    % imshow_sane(PixelVectorToMatrix(Ir_bp,[n m]), ax4, width, width);
    
    imshow_sane(bp_im, ax4, width, width);
% %%
% bpt = bp_im - mean(bp_im(:));
% bpt = min(bpt, hole_depth*.5);
% bpt = max(bpt, -hole_depth*.5);
% bpt = bp_im - mean(bpt(:));
% imshow(bpt, [-hole_depth*.2, .2*hole_depth]);
% figure, plot(bpt(23, :))
%%
% imshow_sane(bpt, ax4, width, width);
%%
    title('BP reconstruction');
end
%
clc
s = metadata2text(ExpMetaData, Ts)
s = sprintf('%s\nperc=%.3f', s, sum(sum(pixelifsampled))/npix/npix);
s2 = sprintf('%s', cs_exp_data_name)
disp(s)
subplot(2,3,[4,5,6])
ax4 = gca();
ax4.Visible = 'off';

t1 = text(0,.5, s, 'Units', 'normalized');
t2 = text(0, t1.Extent(2)-.1, s2, 'Units', 'normalized', 'interpreter', 'none');

% text(0,-1.2, s, 'Units', 'normalized')
%%
savedata = 1;
if savedata
   img_data.cs_im = I;
   img_data.bp_im = bp_im;

   img_data.smp_im = Ir_smp;
   img_data.pixelifsampled = pixelifsampled;
   
   img_data.width = width;
   img_data.meta = ExpMetaData;
   img_data.Ts = Ts;
   img_data_file_name = strrep(cs_exp_data_name, '.csv', '_img-data.mat');
   img_data_path = fullfile(data_root, img_data_file_name);
   
   save(img_data_path, 'img_data')
    
end
%
% fig_root = 'C:\Users\arnold\Documents\labview\afm_imaging\matlab-code\figures';
cs_exp_fig_name = strrep(cs_data_path, '.csv', '-fig.fig')

fig_path = fullfile(cs_exp_fig_name);
saveas(f5, fig_path)

end % MAIN LOOP

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






