clc
clear

% This path should be set to the root of your github working directory. 
if ~ispc
  addpath /home/arnold/gradschool/sysID/matlab/functions
end

addpath(fullfile(getMatPath(), 'afm-cs', 'matlab-code', 'functions'))
addpath(genpath(fullfile(getMatPath(), 'afm-cs', 'reconstruction', 'BP')))
addpath(fullfile(getMatPath(), 'afm-cs', 'reconstruction', 'SMP_1D'))
Ts = 40e-6;


% cs_exp_data_name_s{1} = 'cs-traj-z-bounce_out_10-17-2018-03.csv';
cs_exp_data_name_s{1} = 'cs-traj-512pix-10perc-500nm-5mic-01Hz_out_11-4-2018-01.csv';
sub_dir = '5microns/11-04-2018';
data_root = fullfile(PATHS.exp(), 'imaging', 'cs-imaging', sub_dir);
chan_map = ChannelMap([1:5]);
% -------------------
fname = fullfile(PATHS.sysid, 'x-axis_sines_info_intsamps_zaxisFourierCoef_10-29-2018-01.mat');
models = load(fname);
G = -models.modelFit.G_zdir;
p = pole(G);
z = zero(G);

gg = zpk(z(end-1:end), p(1:2), 1, G.Ts);
% gpz = zpk(p(end-1:end), z(9:10), 1, G.Ts);
% D = zpk(0, 1, 0.025, G.Ts);
% LPF = zpk([], [0.85, 0.85], 1, Ts);
% LPF = LPF/dcgain(LPF);
% H = minreal(ss( (1+D*G)/D)*LPF);
% gg = zpk([], [], 1, Ts);
% gg = gg/dcgain(gg);
% ---------------------
cs_paths = cs_exp_paths(data_root, cs_exp_data_name_s{1});

hole_depth = (20);

cs_exp = CsExp(cs_paths, chan_map, Ts, hole_depth, gg);
fprintf('Total Imaging time: %.2f\n', cs_exp.time_total)
%%
figbase = 100;

Fig_uz = figure(20+figbase); clf

ax1 = gca();
Fig_ze = figure(30+figbase); clf

ax2 = gca();
cs_exp.plot_all_cycles(ax1, ax2);
%%
clc
G = models.modelFit.G_zdir;
[z, p,k] = zpkdata(G, 'v')

zdrift = z(end-1:end);
pdrift = p(1:2);

gdrift = zpk(zdrift, pdrift, 1, G.Ts);
gvib = ss(minreal(G/gdrift));

figure(10); clf
ax1 = gca()
figure(11); clf
ax2 = gca()

cs_exp.fit_gdrift_per_cycle(ax1, ax2, gdrift, gvib);


%%


verbose = false;
if verbose
  figs{1}= figure(1000); clf; hold on, grid on;
  figs{2}= figure(2000); clf; hold on, grid on;
  figs{3}= figure(3000); clf; hold on, grid on;
else
  figs = [];
end
cs_exp.process_cs_data(verbose, figs);

%
pixelifsampled = cs_exp.pix_mask;
I = cs_exp.Img_raw;
fprintf('finished processing raw CS data...\n');


bp = true;
% ********* SMP *************
clear CsExp
cs_exp.solve_smp1d();

ht = cs_exp.feature_height;
figure(10+figbase)
ax = gca();
figure(11+figbase)
axx = gca();
imshow_dataview(cs_exp.Img_raw, [-ht, ht], ax, axx)

figure(12+figbase)
ax = gca();
figure(13+figbase)
axx = gca();
imshow_dataview(cs_exp.Img_smp1d - mean(cs_exp.Img_smp1d(:)), [-ht, ht], ax, axx)


ht = cs_exp.feature_height;
f5 = figure(9 + 10*i); clf
subplot(2,3,1)
ax3 = gca();

subplot(2,3,2)
ax4 = gca();

subplot(2,3,3)
ax5 = gca();
    
subplot(2,3,[4,5,6])
ax6 = gca();
ax6.Visible = 'off';

cs_exp.Img_raw = cs_exp.Img_raw - mean(cs_exp.Img_raw(:));
cs_exp.Img_smp1d = cs_exp.Img_smp1d - mean(cs_exp.Img_smp1d(:));
imshow_sane(cs_exp.Img_raw, ax3, cs_exp.width, cs_exp.width, [-ht, ht]);
title(ax3, 'sample');

imshow_sane(cs_exp.Img_smp1d, ax5, cs_exp.width, cs_exp.width, [-ht, ht])
title(ax5, 'SMP reconstruction');
drawnow
%%
if bp
    cs_exp.solve_basis_pursuit();

    % imshow_sane(PixelVectorToMatrix(Ir_bp,[n m]), ax4, cs_exp.width, cs_exp.width);
    
    
    imshow_sane(cs_exp.Img_bp, ax4, cs_exp.width, cs_exp.width, [-ht, ht]);

    title(ax4, 'BP reconstruction');
end

figure(14+figbase)
ax = gca();
figure(15+figbase)
axx = gca();
imshow_dataview(cs_exp.Img_smp1d - mean(cs_exp.Img_smp1d(:)), [-ht, ht], ax, axx)

%%
axes(ax6)
s = metadata2text(cs_exp.meta_exp, Ts)
s = sprintf('%s\nperc=%.3f', s, sum(sum(pixelifsampled))/cs_exp.npix/cs_exp.npix);
s2 = sprintf('%s', cs_paths.data_path)
disp(s)
ax6.Visible = 'off';

t1 = text(0,.5, s, 'Units', 'normalized');
t2 = text(0, t1.Extent(2)-.1, s2, 'Units', 'normalized', 'interpreter', 'none');


figure(12)
ax = gca();
figure(13)
axx = gca();
imshow_dataview(cs_exp.Img_raw, [-ht, ht], ax, axx)

% text(0,-1.2, s, 'Units', 'normalized')
%%
savedata = 1;
if savedata
   img_data.cs_im = I;
   img_data.bp_im = bp_im;

   img_data.smp_im = Ir_smp;
   img_data.pixelifsampled = pixelifsampled;
   
   img_data.width = width;
   img_data.meta = ExpMetaData;
   img_data.Ts = Ts;
   img_data_file_name = strrep(cs_exp_data_name, '.csv', '_img-data.mat');
   img_data_path = fullfile(data_root, img_data_file_name);
   
   save(img_data_path, 'img_data')
    
end
%
% fig_root = 'C:\Users\arnold\Documents\labview\afm_imaging\matlab-code\figures';
cs_exp_fig_name = strrep(cs_data_path, '.csv', '-fig.fig')

fig_path = fullfile(cs_exp_fig_name);
saveas(f5, fig_path)

