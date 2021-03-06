% After doing the swept sines experiement in LabView, you should run this
% script. It's overall function is to process the sines of the reference
% input, control output, and system output and generate Fourier
% coefficients (by numerically integrating each signal against the base Fourier
% kernal for that frequency). The ratio of these coeffiecients will be the
% different frequency responses (roughly. In actuallity, we do some sneaky
% stuff with all the averages to help with noise). 

% The script will automatically save the resulting FRF's into a .mat file
% with the same name as sinesOut_FileName, replacing only .csv with .mat.


init_paths();

% save FRF data to .mat file?
save_data = true; 
% save the frf figure?
save_fig = false;

% Experiemental data out file. 
%%
FC_data_file_ux = 'x-axis_sines_info_intsamps_quickFourierCoef_11-11-2018xdrive-01.csv';
FC_data_file_uy = 'x-axis_sines_info_intsamps_quickFourierCoef_11-11-2018ydrive-01.csv';
dataRoot = PATHS.sysid;
FC_path_ux = fullfile(dataRoot, FC_data_file_ux);
FC_path_uy = fullfile(dataRoot, FC_data_file_uy);


ss_ux = SweptSinesOnline(FC_path_ux);
ss_uy = SweptSinesOnline(FC_path_uy);

freqs = ss_ux.freq_s_adjusted(:);
assert( length(ss_ux.freq_s_adjusted(:))==  length(ss_uy.freq_s_adjusted(:)))

idx_u = 1;
idx_xdir = 2;
idx_ydir = 3;
idx_zdir = 4;

G_ux2xdir = ss_ux.FC_s(:, idx_xdir)./ss_ux.FC_s(:, idx_u);
G_ux2ydir = ss_ux.FC_s(:, idx_ydir)./ss_ux.FC_s(:, idx_u);
G_ux2zdir = ss_ux.FC_s(:, idx_zdir)./ss_ux.FC_s(:, idx_u);
%
G_uy2xdir = ss_uy.FC_s(:, idx_xdir)./ss_uy.FC_s(:, idx_u);
G_uy2ydir = ss_uy.FC_s(:, idx_ydir)./ss_uy.FC_s(:, idx_u);
G_uy2zdir = ss_uy.FC_s(:, idx_zdir)./ss_uy.FC_s(:, idx_u);

F2 = figure(2); clf
F2.PaperPosition = [1.3376    2.3454    5.8247    6.3093];

h1 = frfBode(G_ux2xdir, freqs, F2,  'Hz', '-k');
h2 = frfBode(G_ux2ydir, freqs, F2,  'Hz', '-g');
h3 = frfBode(G_ux2zdir, freqs, F2,  'Hz', '-r');

h1.DisplayName = 'ux to xdir';
h2.DisplayName = 'ux to ydir';
h3.DisplayName = 'ux to zdir';

subplot(2,1,1)
title('xdir Control signal to stage output')

F3 = figure(3); clf
F2.PaperPosition = [1.3376    2.3454    5.8247    6.3093];

h1y = frfBode(G_uy2xdir, freqs, F3,  'Hz', '-k');
h2y = frfBode(G_uy2ydir, freqs, F3,  'Hz', '-g');
h3y = frfBode(G_uy2zdir, freqs, F3,  'Hz', '-r');

h1y.DisplayName =  'uy to xdir';
h2y.DisplayName =  'uy to ydir';
h3y.DisplayName =  'uy to zdir';

subplot(2,1,1)
title('ydir Control signal to stage output')

legend([h1, h2, h3])
legend([h1y, h2y, h3y])
%%
H = [];
opts = frf2ss_opts('Ts', Ts, 'nd', 9)
H(1,1,:) = G_ux2xdir;
H(2,1,:) = G_ux2ydir;
H(1,2,:) = G_uy2xdir;
H(2,2,:) = G_uy2ydir;

H(3,1,:) = G_ux2zdir;
H(3,2,:) = G_uy2zdir;

opts.r = 350;

fd = frf2ss(H, freqs*2*pi, 10, opts);
figure(4); clf
bar(fd.sigmas)
%%
for k = 1:100
sys = fd.realize(k);
if max(abs(pole(sys))) < 1
  fprintf('stable for ns=%d\n', k)
  kstab = k;
end

end
%%

clc
try
  delete(h_ux_fit)
  delete(h_uy_fit)
end

clrs = {'--b', '--m', '--g'};
clrs = get(gca, 'colororder');
h_ux_fit = gobjects(3,1);
h_uy_fit = gobjects(3,1);
names_x = {'ux to xdir, fit', 'ux to ydir fit', 'ux to zdir fit'};
names_y = {'uy to xdir, fit', 'uy to ydir fit', 'uy to zdir fit'};


for k=1:3
  h_ux_fit(k) = frfBode(sys(k,1), freqs, F2,  'Hz', 'Color', clrs(k,:), 'LineStyle', '--');
  h_ux_fit(k).DisplayName = names_x{k};
end
legend([ h1, h2, h3, h_ux_fit'])
for k=1:3
  h_uy_fit(k) = frfBode(sys(k,2), freqs, F3,  'Hz', 'Color', clrs(k, :), 'LineStyle', '--');
  h_uy_fit(k).DisplayName = names_y{k};
end
legend([ h1y, h2y, h3y, h_uy_fit'])

%%
clc
fc_time_file = 'x-axis_sines_info_intsamps_quick_out_11-11-2018xdrive-01.csv';
time_path = fullfile(dataRoot, fc_time_file);
ss_ux_off = SweptSinesOffline(time_path);
[FC_s, E_s] = SweptSines.fourierCoefs_all_periods(ss_ux_off.xk, ss_ux_off);
ss_ux_off.FC_s = FC_s;

% F2 = figure(2); clf
Gux2xdir2 = FC_s(:,2)./FC_s(:,1);
Gux2ydir2 = FC_s(:,3)./FC_s(:,1); 
Gux2zdir2 = FC_s(:,4)./FC_s(:,1); 

names = {'ux', 'x-dir', 'y-dir', 'z-dir'};

F2 = figure(2); clf
clrs = {'k', 'g', 'r'};
F3 = figure(3);clf;
Gz_list = {Gux2xdir2, Gux2ydir2, Gux2zdir2};
freqs = ss_ux_off.freq_s_adjusted(:);
frf_bode_dataview(Gz_list, ss_ux_off.FC_s*2,freqs, Ts, F2, ss_ux_off.xk, F3, names, 'Hz', clrs);
%%
save('good_frf_xdir.mat', 'Gz_list', 'freqs', 'Ts', 'ss_ux_off', 'names', 'clrs')

%%
% S_ij = SweptSines.getCrossSpectra(ss_ux_off.xk, ss_ux_off);
% Gux2xdir2 = squeeze(S_ij(2,1,:)./S_ij(1,1,:));
% Gux2ydir2 = squeeze(S_ij(3,1,:)./S_ij(1,1,:));
% % Gux2zdir2 = squeeze(S_ij(4,1,:)./S_ij(1,1,:));
% frfBode(Gux2xdir2, ss_ux_off.freq_s_adjusted', F2, 'Hz', '--b');
% frfBode(Gux2ydir2, ss_ux_off.freq_s_adjusted', F2, 'Hz', '--m');
% frfBode(Gux2zdir2, ss_ux_off.freq_s_adjusted', F2, 'Hz', '--y');
%%
model_path = 'xy-axis_sines_info_intsamps_quickFourierCoef_11-11-2018xydrive-01.mat';
names = {'ux2x', 'uy2x';
         'ux2y', 'uy2y';
         'ux2z', 'uy2z'}

try 
  load(model_path)
end

modelFit.frf.all_axes  = H;
modelFit.frf.w_s       = freqs*2*pi;
modelFit.frf.freqs_Hz  = freqs;
modelFit.frf.Ts        = Ts;
modelFit.frf.freq_s    = freqs;
modelFit.frf.names     = names;
modelFit.G_all_axes    = sys;
save_data = true

if save_data
        save(model_path, 'modelFit');
end










