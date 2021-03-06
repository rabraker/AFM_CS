function [ img_data] = csdata2mat(cs_data_path, cs_meta_path, CsExpMetaIn, verbose)

% This path should be set to the root of your github working directory. 
if ispc
    root = 'C:\Users\arnold\Documents\labview';  % For windows
else
    root = '/home/arnold/gradschool/afm-cs';
end
addpath(fullfile(root, 'afm_imaging/matlab-code/functions'));

% TODO: These magic numbers to become inputs or otherwise dissapear!!
pix = CsExpMetaIn.npix; %256;
width = CsExpMetaIn.width %5;
microns_per_volt = 50/10;
pix_per_volt = (pix/width)*microns_per_volt;
Ts = 40e-6;
hole_depth = (20/7)*(1/1000)*(20);



dat = csvread(cs_data_path); 
load(cs_meta_path);  % Provides ExpMetaData

state_ticks = ExpMetaData.state_counts;
state_times = state_ticks*Ts;
time_total = sum(state_times);



% Drop all data corresponding to movement between points.
ind_meas = find(dat(:, 6) > 0);  % Index of the movements are all 0.
meta_ind = dat(:,6);
dat_meas = dat(ind_meas, :);

% Fit a line to the control data and subtract it. Helps remove some of
% the piezo drift. 
uz = dat_meas(:,5);
if verbose > 1
    figure(1); clf;
    plot(uz); hold on
end


dat_meas(:,5) = detrend(dat_meas(:,5)); % detrend u_z as a long vector
plot(dat_meas(:,5))
%
if verbose >1
    f2 = figure(2); clf
    ylabel('u_z')
    xlabel('pixel')
    ax2 = gca();
    hold on;
end
% Since we want to convert xy-coords to pixels, move all data to the
% positive orthant:
dat_meas(:,1) = dat_meas(:,1) - min(dat_meas(:,1));  %x dir
dat_meas(:,2) = dat_meas(:,2) - min(dat_meas(:,2));  %y dir

I = zeros(pix,pix);
pixelifsampled=zeros(pix, pix);
% bin all the data into pixels.
line_data = {}
line_data_avgd = {}
for k = 1:max(dat_meas(:,end))

   % Get indexes corresponding to measurement entity k
   inds = find(dat_meas(:,end) == k); 
   
   % Slice out the data for the current mu-path.
   U_ks = dat_meas(inds, 5);
   Y_ks = dat_meas(inds, 2)*pix_per_volt;
   X_ks = dat_meas(inds, 1)*pix_per_volt;
   line_data{k} = U_ks;
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
   line_data_avgd{k} = U_k;
   % -------------------
   % ---- visualize ------
   if abs(max(U_k) - min(U_k))> hole_depth*.5
       % Then we have an edge. 
        cs = 'r';
    else
        % have_edge = 0;
        cs = 'b';
   end
    if verbose > 1
       plot(ax2, U_k, 'color', cs)
    end
   
end

figure(100);
ax = gca();
if verbose
    imshow(I, [min(min(I)), max(max(I))])
end


% I = detrend_sampled_plane(I, pixelifsampled);
I = (I - max(max(I))).*pixelifsampled; 
if verbose
    figure
    imshow(I, [min(min(I)), max(max(I))]);
end

% imshow_sane(I, ax2, width, width);
% Make the image square, to use smp.

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

reconstruct_root = fullfile(root, 'reconstruction/BP');
reconstruct_path = genpath(reconstruct_root)
addpath(reconstruct_path)

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

s = metadata2text(ExpMetaData, Ts)
if verbose
    disp(s)
end
if verbose >2
% close all;
    f5 = figure(6); clf
    subplot(2,3,1)
    ax3 = gca();
    imshow_sane(I, ax3, width, width);
    title('sample');

    subplot(2,3,2)
    ax4 = gca();
    % imshow_sane(PixelVectorToMatrix(Ir_bp,[n m]), ax4, width, width);
    imshow_sane(bp_im, ax4, width, width);
    title('BP reconstruction');

    subplot(2,3,3)
    ax5 = gca();
    imshow_sane(Ir_smp, ax5, width, width)
    title('SMP reconstruction');
    
    f5.CurrentAxes = ax3;
    text(0,-1.2, s, 'Units', 'normalized')
end


savedata = 1;
if savedata
   img_data.cs_im = I;
   img_data.bp_im = bp_im;
   img_data.smp_im = Ir_smp;
   img_data.pixelifsampled = pixelifsampled;
   
   img_data.width = width;
   img_data.meta = ExpMetaData;
   img_data.Ts = Ts;
    
end
% % %%
% % fig_root = 'C:\Users\arnold\Documents\labview\afm_imaging\matlab-code\figures';
% % cs_exp_fig_name = strrep(cs_exp_data_name, '.csv', '-fig')
% % 
% % fig_path = fullfile(fig_root, cs_exp_fig_name);
% % saveas(f5, fig_path)



end

