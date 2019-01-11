

% Holes are on a  500 nm pitch, and we scan 5 microns, so 10 holes.
% Thus, assuming the holes are half the pitch, each hole is 250 nm and the
% space between is 250 nm. Thus, 1/20th of the image width.
%
init_paths();
npix = 512;
x_start = 26;
y_start = 26;
img_mat = make_CS20NG(x_start, y_start);

figure(14); clf
imshow(img_mat, [0, 1])

%%

% Now, load some cs data. Use the pix_mask, and reconstruct.


img_size = '5microns';
hole_depth = (20);
clear CsExp;
clear CsExp.load_raw_data
clc

% -------------------
fname = fullfile(PATHS.sysid, 'x-axis_sines_info_intsamps_zaxisFourierCoef_10-29-2018-01.mat');
models = load(fname);
G = -models.modelFit.G_zdir;
p = pole(G);
z = zero(G);
gg = zpk(z(end-1:end), p(1:2), 1, G.Ts);


% chan_map = ChannelMap([1:5]);
cs_exp_data_name_s =...
  {'cs-traj-512pix-9perc-500nm-5mic-01Hz_250prescan_out_11-24-2018-01.csv',...
  'cs-traj-512pix-9perc-500nm-5mic-01Hz_250prescan_out_11-24-2018-02.csv',...
  'cs-traj-512pix-9perc-500nm-5mic-01Hz_250prescan_out_11-24-2018-03.csv'};

data_root = PATHS.cs_image_data(img_size, '11-24-2018');

cs_paths1 = get_cs_paths(data_root, cs_exp_data_name_s{1});
cs_exp1 = CsExp(cs_paths1, 'feature_height', hole_depth, 'gg', gg, 'load_full', true);

clear CsExp;
cs_exp1.print_state_times();
fprintf('Total Imaging time: %.2f\n', cs_exp1.time_total)
cs_exp1.process_cs_data();
%%
cs_exp1.solve_bp();
%%

cs_paths2 = get_cs_paths(data_root, cs_exp_data_name_s{2});
cs_exp2 = CsExp(cs_paths2,  'feature_height', hole_depth, 'gg', gg, 'reload_raw', false);
cs_exp2.print_state_times();
fprintf('Total Imaging time: %.2f\n', cs_exp2.time_total)
cs_exp2.process_cs_data();
cs_exp2.solve_basis_pursuit();


cs_paths3 = get_cs_paths(data_root, cs_exp_data_name_s{3});
cs_exp3 = CsExp(cs_paths3,  'feature_height', hole_depth, 'gg', gg, 'reload_raw', false);
cs_exp3.print_state_times();
fprintf('Total Imaging time: %.2f\n', cs_exp3.time_total)
cs_exp3.process_cs_data();
cs_exp3.solve_basis_pursuit();


cs_exp1.save()
cs_exp2.save()
cs_exp3.save()

thresh = (20/7)*(1/1000)*20;
cs_sim_path = '~/matlab/afm-cs/matlab-code/notes/data/cs_sim_CS20NG.mat';

if 0
  cs_sim = CsSim(img_mat*thresh, cs_exp1.pix_mask);
  cs_sim.solve_bp();
%   cs_sim.Img_bp = cs_sim.Img_bp - mean(cs_sim.Img_bp(:));
  save(cs_sim_path, 'cs_sim')
else
  load(cs_sim_path);
  
end

raster_root = PATHS.raster_image_data(img_size, '11-25-2018');
dat_name = 'raster_scan_512pix_5mic_01Hz_out_11-25-2018bungee_extfan-01.csv';
raster_paths = get_raster_paths(raster_root, dat_name);
rast_exp = RasterExp(raster_paths, 'reload_raw', false, 'load_full', true)

rast_exp.bin_raster_really_slow(@detrend);

rast_exp.save()
%%
f2 = figure(2); clf
ImshowDataView.setup(f2);
th_mx_sim = max(cs_sim.Img_bp(:));
th_mn_sim = min(cs_sim.Img_bp(:));

ax1 = subplot(3,1,1:2);
ax2 = subplot(3,1,3);
cb_sim =  @(event_obj)cs_sim.dataview_callback(event_obj, ax1, ax2);
ImshowDataView.imshow(cs_sim.Img_bp,  [th_mn_sim, th_mx_sim], ax1, ax2, cb_sim)
%%
tr_mx = max(cs_exp1.Img_bp(:));
tr_mn = min(cs_exp1.Img_bp(:))

f3 = figure(3); clf
ImshowDataView.setup(f3);
ax3 = subplot(3,1,1:2);
ax4 = subplot(3,1,3);
cb_exp1 =  @(event_obj)cs_exp1.dataview_callback(event_obj, ax3, ax4);
ImshowDataView.imshow(cs_exp1.Img_bp, [tr_mn, tr_mx], ax3, ax4, cb_exp1)

stit = sprintf('Experimental reconstruction. Total time = %.1f', cs_exp1.time_total);
title(ax3, stit)


%%
% disableDefaultInteractivity()
% F11 = figure(11); clf
% set(F11, 'Position', [-880 41 877 956]);
% ypad = 0.01;
% xpad1 = 0.05;
% % Simulation axes
% wd = 0.45;
% wd_im = 0.4;
% 
% bt2 = 0.14;
% bt1 = 0.04;
% bt3 = bt2 + wd_im - 0.0;
% 
% 
% lft_im_a = 0.085;
% lft_im_b = 0.58;
% lft1a = 0.065;
% lft1b = lft1a + xpad1 + wd;
% 
% ht1 = 0.15;
% ax0a = axes('Position', [lft_im_a, bt3, wd_im, 0.5154]);
% ax0b = axes('Position', [lft_im_b, bt3, wd_im, 0.5154]);
% 
% ax1 = axes('Position', [lft_im_a, bt2, wd_im, 0.5154]);
% ax2 = axes('Position', [lft1a bt1 wd, ht1]);
% 
% % Experimental exes
% ax3 = axes('Position', [lft_im_b, bt2, wd_im ,0.5154]);
% ax4 = axes('Position', [lft1b, bt1, wd, ht1]); 
%%
[F11, axs1] = make_fig_4img_2time(11);


clc
% Simulation master
ImshowDataView.imshow(img_mat, [0,1], axs1.top_lft_img, axs1.left_time);
 
% SImulartion CS
cb_sim =  @(event_obj)cs_sim.dataview_callback(event_obj, axs1.bot_lft_img, axs1.right_time);
ImshowDataView.imshow(cs_sim.Img_bp,   [th_mn_sim, th_mx_sim], axs1.bot_lft_img, axs1.right_time, cb_sim)
title(axs1.bot_lft_img, 'Simulated reconstruction') 
% Experiment.
%%
cb_exp2 =  @(event_obj)cs_exp2.dataview_callback(event_obj, axs1.bot_rt_img, axs1.left_time);
ImshowDataView.imshow(cs_exp2.Img_bp, [tr_mn*.8, tr_mx], axs1.bot_rt_img, axs1.left_time, cb_exp2)
stit = sprintf('Experimental reconstruction. Total time = %.1f', cs_exp2.time_total);
title(axs1.bot_rt_img, stit)

xlabel(axs1.left_time, 'pixel')
xlabel(axs1.right_time, 'pixel')
title(axs1.top_lft_img, 'Pure simulation')
 ImshowDataView.setup(gcf)


ImshowDataView.imshow(rast_exp.pix_mat, [-thresh, thresh], axs1.top_rt_img, axs1.left_time);
title(axs1.top_rt_img, '1 Hz raster scan')


% F11 = mkfig(7, 4, 6, true); clf
% ax1 = axes('Position', [0.11438 0.4513 0.775 0.51537]);
% ax2 = axes('Position', [0.13 0.11 0.775 0.27194]);
% ax2.Color = fig_color;
%
% save_fig(F11, 'notes/figures/cs_compare_sim_sim', false)
% %%
% % Now, lets compare, side-by-side, the what happens when we use the
% % TV-smoothing.
%%
%First, we do the simulation:


[fig12, axs2] = make_fig_4img_2time(12);
imagesc(axs2.top_lft_img, cs_sim.Img_bp)
colormap('gray')

imagesc(axs2.top_rt_img, cs_exp2.Img_bp)
colormap('gray')
%%
lambda_s = .01:.1:15;
ssim_tv_s = lambda_s*0;
psnr_tv_s = lambda_s*0;
lambda = 15;
mx_sim = max(abs(cs_sim.Img_bp(:)));
for k=1:length(lambda_s)
  lambda = lambda_s(k);
cs_sim_TV = SplitBregmanROF(cs_sim.Img_bp/mx_sim, lambda, .001)*mx_sim;
% cs_sim_TV = SB_ATV_1dy(cs_sim.Img_bp/mx_sim, lambda)*mx_sim;
% cs_sim_TV = SB_ATV_1dy(cs_sim_TV', .001)';
[psnr_tv, ssim_tv] = metrics(cs_sim.Img_original, cs_sim_TV);
psnr_tv_s(k) = psnr_tv;
ssim_tv_s(k) = ssim_tv;
fprintf('sim (TV): lambda = %f  psnr = %f   ssim=%f\n', lambda, psnr_tv, ssim_tv)

end
figure(100)
yyaxis left
plot(lambda_s, psnr_tv_s)
yyaxis right
plot(lambda_s, ssim_tv_s)
%%


%
lambda = 3;
cs_sim_TV = SplitBregmanROF(cs_sim.Img_bp/mx_sim, lambda, .001)*mx_sim;
[psnr_tv, ssim_tv] = metrics(cs_sim.Img_original, cs_sim_TV);
fprintf('sim (TV): lambda = %f  psnr = %f   ssim=%f\n', lambda, psnr_tv, ssim_tv)

mx_exp2 = max(abs(cs_exp2.Img_bp(:)));
cs_exp_TV = SplitBregmanROF(cs_exp2.Img_bp/mx_exp2, lambda, .001)*mx_exp2;
% cs_exp_TV = SB_ATV_1dy(cs_exp2.Img_bp/mx_exp2, lambda)*mx_exp2;

figure(fig12)
imagesc(axs2.bot_lft_img, cs_sim_TV)
colormap('gray')

imagesc(axs2.bot_rt_img, cs_exp_TV)
colormap('gray')


row_idx = 240;
sim_og_row = cs_sim.Img_bp(row_idx, :);
sim_TV_row = cs_sim_TV(row_idx, :);

exp2_og_row = cs_exp2.Img_bp(row_idx, :);
exp2_TV_row = cs_exp_TV(row_idx, :);

cla(axs2.left_time)
hold(axs2.left_time, 'on')
grid(axs2.left_time, 'on')
plot(axs2.left_time, sim_og_row)
plot(axs2.left_time, sim_TV_row)
plot(axs2.left_time, cs_sim.Img_original(row_idx, :), '--')
ylim(axs2.left_time, 1.05*[min(cs_sim.Img_bp(:)), max(cs_sim.Img_bp(:))] )

cla(axs2.right_time)

hold(axs2.right_time, 'on')
grid(axs2.right_time, 'on')
plot(axs2.right_time, exp2_og_row)
plot(axs2.right_time, exp2_TV_row)
ylim(axs2.right_time, 1.05*[min(cs_exp2.Img_bp(:)), max(cs_exp2.Img_bp(:))] )
%%


function [psnr_ ,ssim_] = metrics(img1, img2)
  dn = max( max(abs(img1(:))), max(abs(img2(:))));
  ssim_ = ssim(img1, img2, 'DynamicRange', dn);
  psnr_ = psnr(img1, img2, dn);


end

function [fig, axs] = make_fig_4img_2time(fig_num)
  fig = figure(fig_num); clf
  set(fig, 'Position', [-880 41 877 956]);
  ypad = 0.03;
  xpad1 = 0.05;
  % Simulation axes
  wd = 0.45;
  wd_im = 0.4;
  
  ht1 = 0.15;
  ht_img = 0.35;
  
  bt1 = 0.04;
  bt2 = bt1 + ht1+ypad;
  
  bt3 = bt2 + wd_im - 0.0;
    
  lft_im_a = 0.085;
  lft_im_b = 0.58;
  lft1a = 0.065;
  lft1b = lft1a + xpad1 + wd;
  
  % masters: sim and raster
  axs.top_lft_img = axes('Position', [lft_im_a, bt3, wd_im, ht_img]);
  axs.top_rt_img = axes('Position', [lft_im_b, bt3, wd_im, ht_img]);
  
  % Simulation
  axs.bot_lft_img = axes('Position', [lft_im_a, bt2, wd_im, ht_img]);
  axs.left_time = axes('Position', [lft1a, bt1, wd, ht1]);

  % Experimental exes
  axs.bot_rt_img = axes('Position', [lft_im_b, bt2, wd_im, ht_img]);
  axs.right_time = axes('Position', [lft1b bt1 wd, ht1]);
  
end