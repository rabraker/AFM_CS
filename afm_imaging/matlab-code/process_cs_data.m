clc
clear
pix = 256;
width = 5;
microns_per_volt = 50/10;
pix_per_volt = (pix/width)*microns_per_volt;
Ts = 40e-6;
hole_depth = (20/7)*(1/1000)*(20);

% dat = csvread('C:\Users\arnold\Documents\labview\afm_imaging/data/cs-data-out02.csv');
% dat = csvread('C:\Users\arnold\Documents\labview\afm_imaging/data/data-out_ontable_csimage.csv');
% dat = csvread('C:\Users\arnold\Documents\labview\afm_imaging\data\data-out_ontable_csimage15-250.csv');
% dat = csvread('C:\Users\arnold\Documents\labview\afm_imaging\data\data-out_ontable_csimage10-500.csv');
% dat = csvread('C:\Users\arnold\Documents\labview\afm_imaging\data\cs-data-10-250.csv');
% dat = csvread('C:\Users\arnold\Documents\labview\afm_imaging\data\data-out_cs-traj10-500-02.csv');  %% WORKS


data_root = 'C:\Users\arnold\Documents\labview\afm_imaging\data\';
cs_exp_data_name = 'cs-traj10-500_8-22-2017_08.csv'; % BEst
cs_exp_data_name = 'cs-traj10-500_out_8-25-2017-14.csv';
% cs_exp_data_name = 'data-out.csv';
cs_exp_meta_name = strrep(cs_exp_data_name, '.csv', '-meta.mat');

cs_data_path = fullfile(data_root, cs_exp_data_name);
cs_meta_path = fullfile(data_root, cs_exp_meta_name);

dat = csvread(cs_data_path); 
load(cs_meta_path);  % Provides ExpMetaData

state_ticks = ExpMetaData.state_counts;
state_times = state_ticks*Ts
time_total = sum(state_times)
%
% This path should be set to the root of your github working directory. 

root = 'C:\Users\arnold\Documents\labview';  % For windows
% root = '/home/arnold/gradschool/afm-cs';
addpath(fullfile(root, 'afm_imaging/matlab-code/functions'));
% dat = csvread(fullfile(root, 'afm_imaging/data/cs-data-out02.csv'));

% Drop all data corresponding to movement between points.
ind_meas = find(dat(:, 6) > 0);  % Index of the movements are all 0.
meta_ind = dat(:,6);
dat_meas = dat(ind_meas, :);

% Fit a line to the control data and subtract it. Helps remove some of
% the piezo drift. 
uz = dat_meas(:,5);
figure(1); clf;
plot(uz); hold on
% uz2 = detrend2(uz);


dat_meas(:,5) = detrend(dat_meas(:,5)); % detrend u_z as a long vector
plot(dat_meas(:,5))
%
f2 = figure(2); clf
ylabel('u_z')
xlabel('pixel')
ax2 = gca();
hold on;

% Since we want to convert xy-coords to pixels, move all data to the
% positive orthant:
dat_meas(:,1) = dat_meas(:,1) - min(dat_meas(:,1));  %x dir
dat_meas(:,2) = dat_meas(:,2) - min(dat_meas(:,2));  %y dir

I = zeros(pix,pix);
pixelifsampled=zeros(pix, pix);
% bin all the data into pixels. 
for k = 1:max(dat_meas(:,end))

   % Get indexes corresponding to measurement entity k
   inds = find(dat_meas(:,end) == k); 
   
   % Slice out the data for the current mu-path.
   U_ks = dat_meas(inds, 5);
   Y_ks = dat_meas(inds, 2)*pix_per_volt;
   X_ks = dat_meas(inds, 1)*pix_per_volt;

   % Register the control data to zero. We can do this because we are
   % scanning long enough that we are always guaranteed to exit a hole.
   % This lets 
   U_ks = U_ks - max(U_ks);
   
   % Make the assumption that the y-data for each path is constant enough.
   % Since we start at the (0,0) corner the xplane, we'll take the floor,
   % and add 1 for 1-based indexing.
   y_pix = floor(mean(Y_ks));
   
   x_spread = max(X_ks) - min(X_ks);
   xpix_start = floor(min(X_ks));
   npix_path_k = ceil(x_spread); % number of bins for this path.
   
   % Now, define a set of bins for the x-direction. Each bin will have a
   % different number of data points. 
   % If we include 0, then there are three points for two pixel bins.
   xbins = linspace(min(X_ks), max(X_ks), npix_path_k+1); 
   U_k = [];
   for jj = 1:npix_path_k
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
       end
   end
   % -------------------
   % ---- visualize ------
   if abs(max(U_k) - min(U_k))> hole_depth*.5
       % Then we have an edge. 
        cs = 'r';
    else
        % have_edge = 0;
        cs = 'b';
    end
   plot(ax2, U_k, 'color', cs)
   
end

figure(100);
ax = gca();
imshow(I, [min(min(I)), max(max(I))])


%%
figure
% I = detrend_sampled_plane(I, pixelifsampled);
I = (I - max(max(I))).*pixelifsampled; 
imshow(I, [min(min(I)), max(max(I))]);
%%

% imshow_sane(I, ax2, width, width);
% Make the image square, to use smp.

size(I)
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

pixelifsampled_vec = PixelMatrixToVector(pixelifsampled);
I_vector = I_vector(find(pixelifsampled_vec>0.5));

A = @(x) IDCTfun(x,pixelifsampled_vec);
At = @(x) DCTfun(x,pixelifsampled_vec);

Ir_bp = idct(l1qc_logbarrier(At(I_vector), A, At, I_vector, 0.1));
Ir_bp = real(Ir_bp);
%%
% close all;
f5 = figure(6); clf
subplot(2,3,1)
ax3 = gca();
imshow_sane(I, ax3, width, width);
title('sample');

subplot(2,3,2)
ax4 = gca();
imshow_sane(PixelVectorToMatrix(Ir_bp,[n m]), ax4, width, width);
title('BP reconstruction');

subplot(2,3,3)
ax5 = gca();
imshow_sane(Ir_smp, ax5, width, width)
title('SMP reconstruction');

clc
% subplot(2,3,[4,5,6])

state_ticks = ExpMetaData.state_counts;
state_times = state_ticks*Ts;
time_total = sum(state_times);

s_time = sprintf('xy-move | z-down | z-settle | xy scan | z-up ||| total');
s_dat  = sprintf('%.2f      %.2f      %.2f      %.2f    %.2f     %.2f    %g',...
                    state_times,  sum(state_times));
s= sprintf('%s\n%s', s_time, s_dat);
s_row = '------------------------------------------------------------------';
for fld =fields(ExpMetaData)'
    if strcmp(fld{1},'state_counts')
        continue
    end
   s_i = sprintf('%s: %g', fld{1}, ExpMetaData.(fld{1}));
   if length(s_row) + length(s_i) < 79-5
       s_row = sprintf('%s  |  %s',s_row, s_i);
   else
       % Row is to long. Tack it onto the whole thing. 
       s = sprintf('%s\n%s', s, s_row);
       s_row = sprintf('%s', s_i);
   end

end

disp(s)
f5.CurrentAxes = ax3;
text(0,-1.2, s, 'Units', 'normalized')
%%
fig_root = 'C:\Users\arnold\Documents\labview\afm_imaging\matlab-code\figures';
cs_exp_fig_name = strrep(cs_exp_data_name, '.csv', '-fig')

fig_path = fullfile(fig_root, cs_exp_fig_name);
saveas(f5, fig_path)


