% The middle of our grating has a 500 nm pitch size. The holes themselves 
% appear to be about 1/4 of this. Thus, i think that if we travel for
% 1/2 * 500 nm = 250 nm, we should be able to detect wheather or not we 
% started in a hole. The idea being that if no jump occures in the u_z 
% signal, we assume we are on the flat, and if we see a jump, we can 
% backout which part is the hole. 
% 
% This will essentially come down to edge detection. If an edge is present,
% then we know the lower part of the signal was in a hole. If no edge is
% present, the entire signal was on the flat part of the grating. 

% close all
% clear, clc
addpath('functions')
try; rmpath('functions/scanning_v0'); end
addpath('functions/scanning_v1')
addpath('functions/state_space_x')


%-------------- Location of System Model -------------------------

plants = CanonPlants.plants_ns14(9);
PLANT_init_x = plants.PLANT;

%-------------- Define scan Params -------------------------

TsTicks = 1600;
Ki_x = 0.01;

if 1
    width =5;  % microns
    npix = 512;  % image resolution.
    mu_length = 0.5;  % 1000 nm. length of the horizontal mu-path. 
    sub_sample_frac = 0.10;  % Percent of pixels to subsample. 
    raster_freq = 1;  % hz
end
N_prescan = 150;
% Unit conversions.
pix_per_micron = npix/width;
mu_pix = ceil(mu_length*pix_per_micron);

% ************************************************
% ************************************************
rng(1)

microns2volts = 1/5;  %10/50;

raster_period = 1/raster_freq;

microns_per_second = width/(raster_period/2);
volts_per_second = microns_per_second*microns2volts;
volts_per_sample = volts_per_second * AFM.Ts;  % Ts = seconds per sample

% Convert mu-path length in pixels to volts
mu_micron = (1/pix_per_micron) * mu_pix;
mu_volts = mu_micron * microns2volts;
mu_Nsamples = ceil(mu_volts / volts_per_sample);


% ************************************************

[pix_mask] = mu_path_mask(mu_pix, npix, npix, sub_sample_frac, false);

actual_sub_sample_frac = length(find(pix_mask == 1))/npix^2;

fprintf('Desired sub sample fraction: %f\n', sub_sample_frac)
fprintf('Actual  sub sample fraction: %f\n', actual_sub_sample_frac)


I = ones(npix,npix)-pix_mask;
figure(1)
imshow(I)

% '------------------correct by overscan for ramp ramp set ---------------
N = mu_Nsamples;

% N-1 because x0 is sample 1.
xdirControl = get_xdir_standard_control('const-sig');
x_N =  (N-1)*volts_per_sample;
N_extra = mu_overscan(xdirControl.Hyr, volts_per_sample, mu_Nsamples, 1, N_prescan)+10

mu_Nsamples_extra = N+N_extra
%%
% ----- Use this to generate the mica CS pattern -------------------
% nr = 20
% mu = [ones(1, 63), 0];
% row = repmat(mu, 1, 8);
% pix_mask_ = repmat(row, nr, 1);
% pix_mask = zeros(npix, npix);
% pix_mask(1:nr, :) = pix_mask_;
% figure(1);clf;imshow(pix_mask)  

mpt = MuPathTraj(pix_mask, width, mu_length, microns_per_second, AFM.Ts,...
  'overscan_samples', N_extra, 'pre_pad_samples', N_prescan);
if 0
  % Magic numbers for connection threshold radius.
  tcon_min = 0.04;
  tcon_est = 0.06;
  rad_min = tcon_min * volts_per_second;
  figure(3);clf; hold on
  ax = gca();  
  mpt_connect_opts = struct('rad_min', rad_min, 'volts_per_sample', volts_per_sample,...
    'ax', ax, 'ensure_forward', true, 'Npasses', N);
  
  mpt.mu_path_connect_rad(mpt_connect_opts);
  
  NN = 0;
  for k=1:length(mpt.mu_path_traj_s)
    NN = NN + length(mpt.mu_path_traj_s{k}.xt);
  end
  
  NN*AFM.Ts + length(mpt.mu_path_traj_s)*tcon_min
end

perc = floor(actual_sub_sample_frac*100);
fname = mpt.get_fname();
target_dir = sprintf('%dmicrons/parents', width);
data_root = fullfile(PATHS.exp(), 'imaging', 'cs-imaging', target_dir);
% data_root = fullfile(PATHS.exp(), 'z-bounce/parents')

fpath_json = fullfile(data_root, fname);
fprintf('File root:\n%s\n', data_root)
fprintf('File name:\n%s\n', fname)


%%
mpt.write_data_json(fpath_json)

%%
% mpt.write_data(fpath_csv)


xt = mpt.mu_path_traj_s{1}.xt;
yt = mpt.mu_path_traj_s{1}.yt;

for k=2:length(mpt.mu_path_traj_s)

  x0 = xt(end);
  y0 = yt(end);

  xt_k1 = mpt.mu_path_traj_s{k}.xt;
  yt_k1 = mpt.mu_path_traj_s{k}.yt;

  x_next = xt_k1(1);
  y_next = yt_k1(1);

  max_dist = max(abs(x0 - x_next), abs(y0-y_next));

  N = max_dist/volts_per_sample;
  x_con = linspace(x0, x_next, N)';
  y_con = linspace(y0, y_next, N)';

  xt = [xt; x_con; xt_k1];
  yt = [yt; y_con; yt_k1];
end

%%
length(xt)*AFM.Ts
figure;
plot(xt)
