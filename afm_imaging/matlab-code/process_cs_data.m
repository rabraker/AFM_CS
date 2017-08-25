clc
clear
pix = 256;
width = 5;
microns_per_volt = 50/10;
pix_per_volt = (pix/width)*microns_per_volt;
Ts = 40e-6;

% dat = csvread('C:\Users\arnold\Documents\labview\afm_imaging/data/cs-data-out02.csv');
% dat = csvread('C:\Users\arnold\Documents\labview\afm_imaging/data/data-out_ontable_csimage.csv');
% dat = csvread('C:\Users\arnold\Documents\labview\afm_imaging\data\data-out_ontable_csimage15-250.csv');
% dat = csvread('C:\Users\arnold\Documents\labview\afm_imaging\data\data-out_ontable_csimage10-500.csv');
% dat = csvread('C:\Users\arnold\Documents\labview\afm_imaging\data\cs-data-10-250.csv');
% dat = csvread('C:\Users\arnold\Documents\labview\afm_imaging\data\data-out_cs-traj10-500-02.csv');  %% WORKS


data_root = 'C:\Users\arnold\Documents\labview\afm_imaging\data\';
cs_exp_data_name = 'cs-traj10-500_8-22-2017_06.csv';
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
% plot(uz2)

dat_meas(:,5) = detrend(dat_meas(:,5)); % detrend u_z as a long vector
plot(dat_meas(:,5))
%%
% Split each measurement entities into a cell array.

% Since we want to convert xy-coords to pixels, move all data to the
% positive orthant:
dat_meas(:,1) = dat_meas(:,1) - min(dat_meas(:,1));  %x dir
dat_meas(:,2) = dat_meas(:,2) - min(dat_meas(:,2));  %y dir


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
   
   % This gives poor results. Think it destroys step height info.
%    U_ks = detrend(dat_meas(inds, 5)); 
   
    
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

% close all
f2 = figure(2); clf
ylabel('u_z')
xlabel('pixel')
ax1 = gca();
hold on;

figure(100); clf; hold on;

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
    
    % We assume that each path contains data from the flat surface.
    % Subtracting off the max is a crude way of registering all the path
    % data to the same height. 
    u_pix = u_pix - max(u_pix);
%     if min(u_pix) > -0.1
        % Just to get a feel for things:
        figure(2)
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
%         figure(100);
%         imshow(I, [min(min(I)), max(max(I))])
%         keyboard
%     else
%         keyboard
%     end

end
%%



I = detrend_sampled_plane(I, pixelifsampled);
% f2 = figure(2); clf
% ax2 = gca();
% 
% imshow_sane(I, ax2, width, width);
% Make the image square, to use smp.
% I = I(1:end-1, 19:end);
size(I)
% pixelifsampled = pixelifsampled(1:end-1, 19:end);
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
f5 = figure(6); clf
subplot(2,3,1)
ax3 = gca();
imshow_sane(I.*PixelVectorToMatrix(pixelifsampled,[n m]), ax3, width, width);
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
cs_exp_fig_name = strrep(cs_exp_data_name, '.csv', '-fig');
fig_path = fullfile(fig_root, cs_exp_fig_name);
saveas(f5, fig_path)


