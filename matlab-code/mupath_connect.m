Fscan = 1;

W = 5e-6;

V = Fscan * 2 * W;

d1 = 1e-6;

Tskip = d1/V


% Time to do zup-xymove,zeng, from slide 42 of comps
Tmove = 0.014; 
Tupdown = 0.066;
Tmu = Tmove + Tupdown;

% Even though Tmove depends on d1, it doesn't change that much. So assume it is
% constant. Then we can solve for d1:

d1_max = Tmu*V
d1_max_mic = d1_max*1e6
%%

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
clear, clc
% close all

addpath('functions')
addpath('functions/state_space_x')


%-------------- Location of System Model -------------------------

% plant_data = load(fullfile(PATHS.sysid(), 'x-axis_sines_infoFourierCoef_9-11-2018-01.mat'));
plants = CanonPlants.plants_ns14(9);
PLANT_init_x = plants.PLANT;

%-------------- Define scan Params -------------------------
TsTicks = 1600;
Ki_x = 0.01;

if 1
    width =5;  % microns
    npix =512;  % image resolution.
    mu_length = 0.5;  % 1000 nm. length of the horizontal mu-path. 
    sub_sample_frac = 0.125;  % Percent of pixels to subsample. 
end
% Unit conversions.
pix_per_micron = npix/width;
mu_pix = ceil(mu_length*pix_per_micron);

% ************************************************
% ************************************************

Ts = AFM.Ts;  % AFM sample rate. 

microns2volts = 1/5;  %10/50;
raster_freq = 1;  % hz
raster_period = 1/raster_freq;

microns_per_second = width/(raster_period/2);
pixels_per_second = pix_per_micron * microns_per_second;
volts_per_second = microns_per_second*microns2volts;
volts_per_sample = volts_per_second * Ts;  % Ts = seconds per sample

% Convert mu-path length in pixels to volts
mu_micron = (1/pix_per_micron) * mu_pix;
mu_volts = mu_micron * microns2volts;
mu_Nsamples = ceil(mu_volts / volts_per_sample);


% ************************************************
rng(1)
[pix_mask, XR, YR] = mu_path_mask(npix, mu_pix, sub_sample_frac, width, microns2volts);

actual_sub_sample_frac = length(find(pix_mask == 1))/npix^2;

fprintf('Desired sub sample fraction: %f\n', sub_sample_frac)
fprintf('Actual  sub sample fraction: %f\n', actual_sub_sample_frac)


I = ones(npix,npix)-pix_mask;
I = flipud(I); % why the fuck does it start from the top???
figure(1)
imshow(I)
%%


% '------------------correct by overscan for ramp ramp set ---------------
G = PLANT_init_x;
D = tf([0.01, 0], [1 -1], Ts);
H = feedback(D*G, 1);
% N-1 because x0 is sample 1.
% x_N =  (mu_Nsamples-1)*volts_per_sample;
N_extra = mu_overscan(G, D, volts_per_sample, mu_Nsamples, 0);

figure(3);clf; hold on
ax = gca();

mpt = MuPathTraj(pix_mask, width, mu_length, microns_per_second, Ts,...
  'overscan_samples', N_extra, 'pre_pad_samples', 250, 'ax', ax);



tcon_min = 0.04;
tcon_est = 0.06;
rad_min = tcon_min * volts_per_second;


figure(3); hold on
ax = gca();
N = 1;
mptc_0 = mpt.mu_path_traj_s;
mpt_connect_opts = struct('rad_min', rad_min, 'volts_per_sample', volts_per_sample,...
        'ax', ax, 'ensure_forward', true, 'Npasses', N);
mpt.mu_path_connect_rad(mpt_connect_opts);


%%
% compute estimated scanning time for the original and connected.
ttot_up_down_mov_og = tcon_est * length(mptc_0);
ttot_scan_og = scan_time(mptc_0);

ttot_og = ttot_scan_og + ttot_up_down_mov_og;

ttot_up_down_mov_con = tcon_est * length(mpt.mu_path_traj_s);
ttot_scan_con = scan_time(mpt.mu_path_traj_s);
ttot_con = ttot_scan_con + ttot_up_down_mov_con;


fprintf('Original Scan Time (estimated): %f\n', ttot_og);
fprintf('Connected Scan Time (estimated): %f\n', ttot_con);
fprintf('Estimated percent savings: %f\n', (1-ttot_con/ttot_og)*100);


perc = floor(actual_sub_sample_frac*100);
fname = sprintf('cs-traj-%dpix-%dperc-%dnm-%dmic-%.2dHz_250prescan.json',npix, perc, mu_length*1000,width, raster_freq);
target_dir = sprintf('%dmicrons/parents', width);
data_root = fullfile(PATHS.exp(), 'imaging', 'cs-imaging', target_dir);
fpath_json = fullfile(data_root, fname);

fprintf('File root:\n%s\n', data_root)
fprintf('File name:\n%s\n', fname)


mpt.write_data_json(fpath_json)
% % figure(5); clf
% % for k=1:length(mpt.mu_path_traj_s)
% %    hold on;
% %   plot(mpt.mu_path_traj_s{k}.met_idx)
% % end

function t_scan = scan_time(mptc_s)
  
  t_scan = 0;
  for k=1:length(mptc_s)
    t_scan = t_scan + length(mptc_s{k}.xt)*AFM.Ts;
  end
end
