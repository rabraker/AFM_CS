%
% Generate "mu-path like" trajectories, where we don't actually move anywhere,
% but the point is to test the z-up-down scheme.


clear, clc

addpath('classes')
addpath('functions')
clc
%-------------- Location of System Model -------------------------
samps_per_path = 1500;
pre_scan_samples = 250;
Npaths = 160
mpb = MuPathBounce(Npaths, samps_per_path, 'pre_pad_samples', pre_scan_samples)


% Unit conversions.

clc
fname = sprintf('cs-traj-z-bounce_%dsamples_Npath%d_prescan%d.json', samps_per_path, Npaths, pre_scan_samples);
fname = fullfile(PATHS.exp, 'z-bounce', 'parents', fname)

mpb.write_data_json(fname)



