%
% Generate "mu-path like" trajectories, where we don't actually move anywhere,
% but the point is to test the z-up-down scheme.


clear, clc

addpath('classes')
addpath('functions')
clc
%-------------- Location of System Model -------------------------
samps_per_path = 5000;
mpb = MuPathBounce(5, samps_per_path)


% Unit conversions.

clc
fname = sprintf('cs-traj-z-bounce_%dsamples.json', samps_per_path);
fname = fullfile(PATHS.exp, 'z-bounce', 'parents', fname)

mpb.write_data_json(fname)



