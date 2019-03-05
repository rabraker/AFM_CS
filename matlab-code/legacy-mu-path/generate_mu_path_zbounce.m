%
% Generate "mu-path like" trajectories, where we don't actually move anywhere,
% but the point is to test the z-up-down scheme.


clear, clc

addpath('classes')
addpath('functions')

%-------------- Location of System Model -------------------------
dataroot      = fullfile(getMatPath(), 'AFM_SS',...
                'System_Identification', 'data', 'data_xaxis'); 
expName           = ['22-Jun-2016_exp01'];
modFitName    = [expName, '.mat'];
modFitPath    = fullfile(dataroot, modFitName);
load(modFitPath, 'modelFit')


% FitNum    = 'FitNum001';
PLANT_init_x = ltiFit(modFitPath, 'SS02').sys;
%-------------- Define scan Params -------------------------

Ts = 40e-6;
TsTicks = 1600;
Ki_x = 0.01;

if 1
    width = 5;  % microns
    pix =512;  % image resolution.
    mu_length = .5;  % 500 nm. length of the horizontal mu-path. 
    sub_sample_frac = 0.20;  % Percent of pixels to subsample. 
end
% Unit conversions.
pix_per_micron = pix/width;
mu_pix = ceil(mu_length*pix_per_micron);

% ************************************************
% ************************************************

Ts = 40e-6;  % AFM sample rate. 
Fs = 1/Ts;
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

ME = MeasEntityMu.factory([volts_per_sample*0, 0]);

N_mve = 1;
% ************************************************

XR = zeros(100,1);
YR = XR;

% '------------------correct by overscan for ramp ramp set ---------------
N = mu_Nsamples;
x_rate = volts_per_sample*0;


N_extra = 500;
mu_Nsamples_extra = N+N_extra;
meta_cell = repmat({mu_Nsamples_extra}, 1, length(XR));
MT = MasterTrajster(XR*0, YR*0, meta_cell, MoveEntityStatic.factory(N_mve), ME);



MT.visualize_sampling;

xlabel('x [v]')
ylabel('y [v]')
grid on

fname = 'cs-traj-z-bounce_500nmequiv.csv';

% create meta file name
% target_dir = sprintf('%dmicrons/parents', width)
if ispc
  data_root = 'Z:\afm-cs\z-bounce\parents';
else
  data_root = '/media/labserver/afm-cs/z-bounce/parents';
end

meta_in = strrep(fname, '.csv', '.mat');
meta_data_path = fullfile(data_root, meta_in);

%%


tip_velocity = width/(0.5*(1/raster_freq));

CsExpMetaIn.width = width;
CsExpMetaIn.nom_perc = sub_sample_frac;
CsExpMetaIn.mu_length = mu_length;
CsExpMetaIn.tip_velocity = tip_velocity;
CsExpMetaIn.npix = pix
CsExpMetaIn.pixelifsampled = 0; %pixelifsampled;
CsExpMetaIn.actual_perc = 0; %perc;

datafile = fullfile(data_root, fname)

if 1
    MT.write_csv(datafile)
    save(meta_data_path, 'CsExpMetaIn')
end
   



    

