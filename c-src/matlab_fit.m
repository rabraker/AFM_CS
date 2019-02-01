
% init_paths

freq0 = [];
freq1 = linspace(1, 190, 40)';
freq2 = linspace(191, 235, 45)';
% freq3 = linspace(236, 3
freq3 = logspace(log10(235), log10(12500-1), 100)';

freqs = [freq0; freq1; freq2; freq3];
rand_fname = fullfile(PATHS.sysid, 'rand_noise_zaxis_10-30-2018_01.csv');

[G_frf, freqs, Coh, params] = load_randnoise_frf(rand_fname);
G_frf(1) = [];
freqs(1) = [];

%%
f1 = figure(2001);

frfBode(G_frf, freqs, f1, 'Hz')


idx = find(freqs < 264 & freqs > 160);
G_res = G_frf(idx);
freqs_res = freqs(idx);

%%
clc


Ts = 40e-6;
wn1 = 210.6*2*pi;
wn2 = 215.1*2*pi;
z = 0.01;
go = tf([1, 2*z*wn2, wn2^2], [1, 2*z*wn1, wn1^2]);
go = c2d(go, Ts, 'matched');
go = zpk((10^(-3/20)) * go /dcgain(go))

sos_fos = SosFos(go, 'iodelay', 2);
g1 = sos_fos.realize()

%
LG = LogCostZPK(-G_res, freqs_res*2*pi, sos_fos);

LG.solve_lsq(1);

[g2, p]= LG.sos_fos.realize();
fprintf('LG says delay = %f\n', p)
g2.IODelay = floor(p);




f2 = figure(2002); clf
frfBode(-G_res, freqs_res, f2, 'hz', 'r')
frfBode(g1, freqs_res, f2, 'Hz', 'g')
frfBode(g2, freqs_res, f2, 'Hz', '--b')

% % compare
% fname = fullfile(PATHS.sysid, 'x-axis_sines_info_intsamps_zaxisFourierCoef_10-29-2018-01.mat');
% models = load(fname);
% 
% 
% Gres_old = 1/models.modelFit.Dinv;
% f3 = figure(3)
% frfBode(Gres_old, freqs_res, f3, 'Hz', '-k')
% frfBode(g2/dcgain(g2), freqs_res, f3, 'Hz', '--r')

%%
% write json data

opts.FileName = 'frf_test_data.json'
theta = [0.6786, -1.9959, 0.99891, -1.996142, 0.9989];
savejson('', struct('omegas', omegas(:)', 'resp_real', real(resp(:)'),...
  'resp_imag', imag(resp(:)'), 'theta0', theta), opts)


%%
omegas = freqs_res*2*pi;
resp = -G_res;
WH = WriteHeader('include/test_data')

WH.open();
WH.write_ifndef()

WH.write_test_data(omegas, 'omegas');
WH.write_test_data(real(resp), 'resp_real');
WH.write_test_data(imag(resp), 'resp_imag');

WH.close_ifdef()

WH.close();

