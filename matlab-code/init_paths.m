
function init_paths()
% init_paths()
% Initialize paths for for working on afm-cs
% ----------------------------------------------------------------------- %
% Run this script to initialize paths 
% ----------------------------------------------------------------------- %

addpath(fullfile(getCsRoot, 'matlab-code', 'functions'));
addpath(fullfile(getCsRoot, 'reconstruction', 'BP', 'l1magic', 'Optimization'));
addpath(fullfile(getCsRoot, 'reconstruction', 'BP'));
addpath(fullfile(getCsRoot, 'reconstruction', 'SMP_1D'));

addpath(fullfile(getCsRoot, 'matlab-code', 'models'))


end