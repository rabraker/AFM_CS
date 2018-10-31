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

% clear
clc


% Experiemental data out file. Although we need the data from the
% input-data csv file, labview saves that file into the header of the
% output file. Labview will also save the data from the z-axis-settings
% cluster into the header. The following variable is the filename displayed
% in 'output-data-file-name' indicator in play_sysID_Z_Axis.vi'.


freq0 = [];
freq1 = linspace(1, 190, 40)';
freq2 = linspace(191, 235, 45)';
% freq3 = linspace(236, 3
freq3 = logspace(log10(235), log10(12500-1), 100)';

freqs = [freq0; freq1; freq2; freq3];
rand_fname = fullfile(PATHS.sysid, 'rand_noise_zaxis_10-30-2018_01.csv');

[G_uz2stage_rnd, freqs_rnd, Coh, params] = load_randnoise_frf(rand_fname);
G_uz2stage_rnd(1) = [];
freqs_rnd(1) = [];


% For saving results:
frf_FileName = strrep(rand_fname, '.csv', '.mat');
% save FRF data to .mat file?
save_data = true; 
% save the frf figure?
save_fig = false;


Ts = params.Ts;

F1 = figure(1); clf
F1.PaperPosition = [1.3376    2.3454    5.8247    6.3093];
h1 = frfBode(G_uz2stage_rnd, freqs_rnd, F1,  'Hz', '-r');
subplot(2,1,1)
title('Control signal to stage output')

G_uz2stage = interp1(freqs_rnd, G_uz2stage_rnd, freqs, 'spline');
frfBode(G_uz2stage, freqs, F1, 'Hz', '--k')
% F2 = figure(2); clf
% F2.PaperPosition = [1.3376    2.3454    5.8247    6.3093];


%%
F2 = figure(2); clf
F2.PaperPosition = [1.3376    2.3454    5.8247    6.3093];
h1 = frfBode(G_uz2stage, freqs, F2,  'Hz', '-r');
subplot(2,1,1)
title('Control signal to stage output')

if 1
  gk = load('../g_k.mat')
  gk = gk.g_k;
else
  gk = zpk([], [], 1, Ts);
end
clc
ss_opts = frf2ss_opts('Ts', Ts, 'r', 500, 's', 300);
Nd2 = 1;
ns2 = 14;
f2ss = frf2ss(G_uz2stage, freqs*2*pi, Nd2, ss_opts); % 12
sys_z = f2ss.realize(ns2)*gk; % 12

% sys_z = eject_nmpz(sys_z);

h1 = frfBode(sys_z, freqs, F2,  'Hz', '-k');

k_estmax = find(freqs > 6000, 1, 'first');
subplot(2,1,1)
ylm = ylim;
subplot(2,1,2)
ylim([-180*3, 180])
plot([freqs(k_estmax), freqs(k_estmax)], ylm, 'k');
% plotPZ_freqs(sys_z, F2);

clc
LGopts = optimoptions(@lsqnonlin, 'Display', 'iter',...
    'FunctionTolerance', 1e-9, 'MaxIter', 5000,'MaxFunctionEvaluations', 5000,...
    'StepTolerance', 1e-9, 'Jacobian','on', 'CheckGradients', false);

sos_fos = SosFos(sys_z, 'iodelay', sys_z.InputDelay);
LG = LogCostZPK(G_uz2stage(1:k_estmax), freqs(1:k_estmax)*2*pi, sos_fos);
W = ones(length(LG.omegas),1);
W(LG.omegas < 2*pi*300) = 15;
LG.solve_lsq(2, LGopts, W)
[sys_stage_log, p] = LG.sos_fos.realize();
sys_stage_log.InputDelay = 1; %max(round(p, 0), 0);
fprintf('LG says delay = %.2f\n', p);

stab = isstable(sys_stage_log);


if stab
  fprintf('System is stable.\n');
else
  fprintf('System is NOT STABLE.\n');
end

frfBode(sys_stage_log, freqs, F2,  'Hz', '--g');
plotPZ_freqs(sys_stage_log, F2);

size(zero(sys_stage_log))
sys_log = eject_nmpz(sys_stage_log);
size(zero(sys_log))
sys_log.InputDelay = 0;
sys_log.IODelay = 0;
frfBode(sys_log, freqs, F2,  'Hz', '-b');


%%
clc
F3 = figure(3); clf
F3.PaperPosition = [1.3376    2.3454    5.8247    6.3093];
h1 = frfBode(G_uz2stage, freqs, F3,  'Hz', '-r');
subplot(2,1,1)
title('Control signal to stage output')

opts = balredOptions('StateElimMethod', 'Truncate', 'FreqIntervals', [0, 2000]*2*pi)
sys_log2 = balred(sys_log, 5, opts)*gk;
sys_log2 = eject_nmpz(sys_log2);
frfBode(sys_log2, freqs, F3, 'Hz', '-k');
plotPZ_freqs(sys_log2, F3);


k_estmax = find(freqs > 2000, 1, 'first');
sos_fos = SosFos(sys_log2, 'iodelay', 1);
LG = LogCostZPK(G_uz2stage(1:k_estmax), freqs(1:k_estmax)*2*pi, sos_fos);
LG.solve_lsq(2, LGopts)
[sys_log3, p] = LG.sos_fos.realize();

sys_log3.InputDelay = max(0, round(p, 0));
fprintf('LG says delay = %.2f\n', p);

stab = isstable(sys_log3);
%%

if stab
  fprintf('System is stable.\n');
else
  fprintf('System is NOT STABLE.\n');
end

frfBode(sys_log3, freqs, F3,  'Hz', '--g');
plotPZ_freqs(sys_log3, F3);

%%
G = sys_log;
% pzplot(G)
p = pole(G);
% p(1:2) = [];
z = zero(G)
% z(end-1:end) = [];

gdrift = zpk(z(end-1:end), p(1:2), 1, Ts);
gdrift = gdrift/dcgain(gdrift);

figure(15); clf
step(gdrift)
hold on
% figure(16)
step(-G)
% G = zpk(z, p, 1, Ts);
% G = G*dcgain(sys_log)/dcgain(G);
% hold on, pzplot(G)
%%
clc
z = zero(G)
p = pole(G)

zz = z(end-4:end-3)
pp = p(end-1:end)
Dinv = zpk(pp, zz, 1, Ts);
Dinv = Dinv/dcgain(Dinv);

%%

KI = -0.025
D_I = zpk(0, 1, KI, Ts)

D = D_I * Dinv;
loop = D*G;
F10 = figure(10); clf
frfBode(loop, freqs, F10, 'Hz', 'k')
frfBode(G, freqs, F10, 'Hz', 'b')

H = -minreal(feedback(D, G))*gdrift;
isstable(H)
frfBode(H, freqs, F10, 'Hz', 'r')

figure(5)
margin(loop)

figure(4)
step(H)
%%

Loop = Dinv*D_I*sys_log;
H_yr = minreal(Loop/(1+Loop));
isstable(H_yr)
F10 = figure(10); clf

frfBode(Loop, freqs, F10, 'Hz', '-b')
frfBode(H_yr, freqs, F10, 'Hz', '-k')


%%
% model_path = strrep(FC_path, '.csv', '.mat');

try 
  load(frf_FileName)
end

modelFit.frf.G_uz2stage = G_uz2stage;
modelFit.frf.w_s       = freqs*2*pi;
modelFit.frf.freqs_Hz  = freqs;
modelFit.frf.Ts        = Ts;
modelFit.frf.freq_s    = freqs;
modelFit.G_zdir        = sys_log;
modelFit.G_reduce      = sys_log3;
modelFit.Dinv          = Dinv;
modelFit.gdrift        = gdrift;
if save_data
        save(frf_FileName, 'modelFit');
end










