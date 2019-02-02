
% init_paths

freq0 = [];
freq1 = linspace(1, 190, 40)';
freq2 = linspace(191, 235, 45)';
% freq3 = linspace(236, 3
freq3 = logspace(log10(235), log10(12500-1), 100)';

freqs = [freq0; freq1; freq2; freq3];
% ss_fname = fullfile(PATHS.sysid, 'z-axis_sines_info_quick_firstResFourierCoef_11-26-2018-03.json');
ss_fname = ' Z:\afm-cs\sysID\ALL-axis_sines_info_intsamps_quickFourierCoef_1-17-2019z-drive-03.json';
ss_data = SweptSinesOnline(ss_fname);

% [G_frf, freqs, Coh, params] = load_randnoise_frf(rand_fname);
% G_frf(1) = [];
% freqs(1) = [];
G_frf = ss_data.FC_s(:,4)./ss_data.FC_s(:,1);
freqs = ss_data.freq_s;
%%
f1 = figure(2001);

frfBode(G_frf, freqs, f1, 'Hz')


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
LG = LogCostZPK(-G_frf, freqs*2*pi, sos_fos);

LG.solve_lsq(2);

[g2, p]= LG.sos_fos.realize();
fprintf('LG says delay = %f\n', p)
g2.IODelay = floor(p);




f2 = figure(2002); clf
frfBode(-G_frf, freqs, f2, 'hz', 'r')
% frfBode(g1, freqs, f2, 'Hz', 'g')
frfBode(g2, freqs, f2, 'Hz', '--b')
% frfBode(11models.modelFit.Dinv, freqs, f2, 'Hz', ':k')
%%
g_1st_res = g2
mat_file = strrep(ss_fname, '.json', '.mat');
save(mat_file, 'ss_data', 'g_1st_res');
%%

