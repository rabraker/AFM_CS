clc
clear

% initialize paths.
init_paths();

Ts = 40e-6;

size = '5microns';
% cs_exp_data_name_s{1} = 'cs-traj-512pix-9perc-500nm-5mic-01Hz_250prescan_out_11-26-2018newDinv-03.csv';

% data_root = PATHS.cs_image_data(size, '11-26-2018');

% This one is pretty decent. move   |  lower  |  settle  | scan   | up 
% 4.885 | 1.376  | 10.100   | 28.671  | 0.643 |
% Total Imaging time: 45.68
% cs_exp_data_name_s{1} = 'cs-traj-512pix-10perc-500nm-5mic-01Hz_v2_out_11-17-2018-02.csv';
% data_root = PATHS.cs_image_data(size, '11-17-2018');


%move   |  lower  |  settle  | scan   | up 
% 4.035 | 0.971  | 9.658   | 22.607  | 0.428 |
% Total Imaging time: 37.70
% cs_exp_data_name_s{1} = 'cs-traj-512pix-9perc-500nm-5mic-01Hz_250prescan_out_11-24-2018-03.csv';
data_root = PATHS.cs_image_data(size, '11-24-2018');
cs_exp_data_name_s{1} = 'cs-traj-256pix-10perc-500nm-5mic-01Hz_250prescan.json_2-14-2019-03.csv';
data_root = PATHS.cs_image_data(size, '2-14-2019');


chan_map = ChannelMap([1:5]);
% -------------------
fname = fullfile(PATHS.sysid, 'x-axis_sines_info_intsamps_zaxisFourierCoef_10-29-2018-01.mat');
models = load(fname);
G = -models.modelFit.G_zdir;
p = pole(G);
z = zero(G);
%%
close all
gg = zpk(z(end-1:end), p(1:2), 1, G.Ts);
% gpz = zpk(p(end-1:end), z(9:10), 1, G.Ts);
% D = zpk(0, 1, 0.025, G.Ts);
% LPF = zpk([], [0.85, 0.85], 1, Ts);
% LPF = LPF/dcgain(LPF);
% H = minreal(ss( (1+D*G)/D)*LPF);
% gg = zpk([], [], 1, Ts);
% gg = gg/dcgain(gg);
% ---------------------
% gg = [];
cs_paths = get_cs_paths(data_root, cs_exp_data_name_s{1});

hole_depth = (20);

cs_exp = CsExp(cs_paths, 'feature_height', hole_depth, 'gg', gg);

% cs_exp.print_state_times();
% fprintf('Total Imaging time: %.2f\n', cs_exp.time_total)

% % something screwed up. Was supposed to only have 250 pre-scan samples, but
% % they are in the settle.
% for k=1:length(cs_exp.idx_state_s.tsettle)
%   idx_k = cs_exp.idx_state_s.tsettle{k}(end-200:end);
%   cs_exp.idx_state_s.tsettle{k}(end-200:end) = [];
%   cs_exp.idx_state_s.scan{k} = [idx_k, cs_exp.idx_state_s.scan{k}];
% end
%
figbase = 20;

[~, axs] = make_traj_figs(figbase);
cs_exp.plot_all_cycles(axs{:});

%%

verbose = false;
if verbose
  fig_inc = 10;
  figs{1}= figure(1000+fig_inc); clf; hold on, grid on;
  figs{2}= figure(2000+fig_inc); clf; hold on, grid on;
  figs{3}= figure(3000+fig_inc); clf; hold on, grid on;
else
  figs = [];
end

cs_exp.process_cs_data(verbose, figs);
pixelifsampled = cs_exp.pix_mask;
I = cs_exp.Img_raw;
fprintf('finished processing raw CS data...\n');

fprintf('nperc=%.3f\n', sum(cs_exp.pix_mask(:))/cs_exp.npix^2);

ht = cs_exp.feature_height;
figure(10+figbase)
ax = gca();
figure(11+figbase)
axx = gca();
imshow_dataview(cs_exp.Img_raw-mean(cs_exp.Img_raw(:)), [-ht, ht], ax, axx)


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


if bp
    cs_exp.solve_basis_pursuit();
    imshow_sane(cs_exp.Img_bp, ax4, cs_exp.width, cs_exp.width, [-ht, ht]);
    title(ax4, 'BP reconstruction');
end

figure(14+figbase)
ax = gca();
figure(15+figbase)
axx = gca();
imshow_dataview(cs_exp.Img_bp - mean(cs_exp.Img_bp(:)), [-ht, ht], ax, axx)


 
% for k=1:length(cs_exp.idx_state_s.scan)
%   idx_k = cs_exp.idx_state_s.scan{k};
%   xmin_k = max(self.x(idx_k));
%   ymin_k = max(self.y(idx_k));
%   xmin_meas = max([xmin_meas, xmin_k]); % brackets necessary
%   ymin_meas = min([ymin_meas, ymin_k]); % brackets necessary
% end




axes(ax6)
s = metadata2text(cs_exp.meta_exp, Ts)
s = sprintf('%s\nperc=%.3f', s, sum(sum(pixelifsampled))/cs_exp.npix/cs_exp.npix);
s2 = sprintf('%s', cs_paths.data_path)
disp(s)
ax6.Visible = 'off';

t1 = text(0,.5, s, 'Units', 'normalized');

% pixmat2 = cs_exp.Img_bp;
% save('Z:\afm-cs\tuesday-figs\11-5-2016\cs_img_slowdown.mat', 'pixmat2')


%% Compare to ACC
acc_10_bp = get_acc_cs_img(10)
figure(214)
ax1 = gca();

figure(2)
ax2 = gca();

acc_10_bp.thresh

imshow_dataview(acc_10_bp.bp_im_normalized, [-ht, ht], ax1, ax2)
%%
figure(230)
imshow(1-cs_exp.pix_mask)
%%
figure(240)
imshow(1-cs_exp.meta_in.pix_mask)
figure(250)
imshow(cs_exp.meta_in.pix_mask - cs_exp.pix_mask, [-1, 1])
%%
sum(cs_exp.meta_in.pix_mask(:))/512^2
sum(cs_exp.pix_mask(:))/512^2

%%
% figure(12)
% ax = gca();
% figure(13)
% axx = gca();
% imshow_dataview(cs_exp.Img_raw, [-ht, ht], ax, axx)

% text(0,-1.2, s, 'Units', 'normalized')
%%
% savedata = 1;
% if savedata
%    img_data.cs_im = I;
%    img_data.bp_im = bp_im;
% 
%    img_data.smp_im = Ir_smp;
%    img_data.pixelifsampled = pixelifsampled;
%    
%    img_data.width = width;
%    img_data.meta = ExpMetaData;
%    img_data.Ts = Ts;
%    img_data_file_name = strrep(cs_exp_data_name, '.csv', '_img-data.mat');
%    img_data_path = fullfile(data_root, img_data_file_name);
%    
%    save(img_data_path, 'img_data')
%     
% end
% %
% % fig_root = 'C:\Users\arnold\Documents\labview\afm_imaging\matlab-code\figures';
% cs_exp_fig_name = strrep(cs_data_path, '.csv', '-fig.fig')
% 
% fig_path = fullfile(cs_exp_fig_name);
% saveas(f5, fig_path)

function [figs, axs] = make_traj_figs(figbase)
  Fig_uz = figure(20+figbase); clf
  
  ax1 = gca();
  Fig_ze = figure(30+figbase); clf
  ax2 = gca();
  
  Fig_x = figure(40+figbase); clf
  ax3 = gca();
  Fig_y = figure(50+figbase); clf
  ax4 = gca();
  
  figs = {Fig_uz, Fig_ze, Fig_x, Fig_y};
  axs = {ax1, ax2, ax3, ax4};
end
