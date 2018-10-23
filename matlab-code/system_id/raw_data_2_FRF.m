% After doing the swept sines experiement in LabView, you should run this
% script. It's overall function is to process the sines of the reference
% input, control output, and system output and generate Fourier
% coefficients (by numerically integrating each signal against the base Fourier
% kernal for that frequency). The ratio of these coeffiecients will be the
% different frequency responses (roughly. In actuallity, we do some sneaky
% stuff with all the averages to help with noise). 

% The script will automatically save the resulting FRF's into a .mat file
% with the same name as sinesOut_FileName, replacing only .csv with .mat.

if ispc
  rmpath('C:\Users\arnold\Documents\MATLAB\miscScripts\system_id\')
  addpath('C:\Users\arnold\Documents\labview\sysID\matlab\functions')
else
  rmpath(fullfile(getMatPath(), 'toolboxes', 'system_id'))
  addpath('/home/arnold/gradschool/sysID/matlab/functions')
end

% addpath('functions')

clear
clc


% Experiemental data out file. Although we need the data from the
% input-data csv file, labview saves that file into the header of the
% output file. Labview will also save the data from the z-axis-settings
% cluster into the header. The following variable is the filename displayed
% in 'output-data-file-name' indicator in play_sysID_Z_Axis.vi'.



FC_data_file = 'x-axis_sines_infoFourierCoef_10-21-2018-03.csv';
dataRoot = PATHS.sysid;
FC_path = fullfile(dataRoot, FC_data_file);

% For saving results:
frf_FileName = strrep(FC_data_file, '.csv', '.mat');
% save FRF data to .mat file?
save_data = true; 
% save the frf figure?
save_fig = false;


ssOpts = sweptSinesMeta('read', FC_path);

Ts = ssOpts.Ts;

[FC_s, E_s, freqs, ssOpts] = SweptSines.read_FC_data(FC_path, ssOpts);
freqs = freqs(:);
idx_uz = 1;
idx_stage = 2;



G_uz2stage = FC_s(:, idx_stage)./FC_s(:, idx_uz);
figure(100)
semilogx(freqs, E_s(:, idx_uz))

%%
F2 = figure(2); clf
F2.PaperPosition = [1.3376    2.3454    5.8247    6.3093];
h1 = frfBode(G_uz2stage, freqs, F2,  'Hz', '-r');
subplot(2,1,1)
title('Control signal to stage output')

% F2 = figure(2); clf
% F2.PaperPosition = [1.3376    2.3454    5.8247    6.3093];
% 
% rand_fname = fullfile(PATHS.sysid, 'rand_noise_frf_z_axis_10-22-2018.csv');
% [G_uz2stage, freqs] = load_randnoise_frf(rand_fname);
% G_uz2stage(1) = [];
% freqs(1) = [];
% frfBode(G_uz2stage, freqs, F2, 'Hz', '--k')
%%
gk = load('../g_k.mat')
gk = gk.g_k;

clc
ss_opts = frf2ss_opts('Ts', Ts, 'r', 500, 's', 300);
Nd2 = 2;
ns2 = 12;
f2ss = frf2ss(G_uz2stage, freqs*2*pi, Nd2, ss_opts); % 12
sys_z = f2ss.realize(ns2)*gk; % 12
frfBode(sys_z, freqs, F2, 'Hz', '-k')

k_estmax = find(freqs > 4000, 1, 'first');
subplot(2,1,1)
ylm = ylim;
subplot(2,1,2)
ylim([-180*3, 180])
plot([freqs(k_estmax), freqs(k_estmax)], ylm, 'k');
% plotPZ_freqs(sys_z, F2);
%%
LGopts = optimoptions(@lsqnonlin, 'Display', 'iter',...
    'FunctionTolerance', 1e-7, 'MaxIter', 5000,'MaxFunctionEvaluations', 5000,...
    'StepTolerance', 1e-8, 'Jacobian','off'); %, 'CheckGradients', false);

sos_fos = SosFos(sys_z, 'iodelay', sys_z.InputDelay);
LG = LogCostZPK(G_uz2stage(1:k_estmax), freqs(1:k_estmax)*2*pi, sos_fos);
LG.solve_lsq(3, LGopts)
[sys_stage_log, p] = LG.sos_fos.realize();
sys_stage_log.InputDelay = max(round(p, 0), 0);
fprintf('LG says delay = %.2f\n', p);

frfBode(sys_stage_log, freqs, F2,  'Hz', '--g');
plotPZ_freqs(sys_stage_log, F2);
%%
model_path = strrep(FC_path, '.csv', '.mat');

try 
  load(model_path)
end

modelFit.frf.G_uz2stage = G_uz2stage;
modelFit.frf.w_s       = freqs*2*pi;
modelFit.frf.freqs_Hz  = freqs;
modelFit.frf.Ts        = Ts;
modelFit.frf.freq_s    = freqs;
modelFit.G_zdir        = sys_stage_log;


if save_data
        save(model_path, 'modelFit');
end










