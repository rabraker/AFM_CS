clc
clear
%%
close all
addpath functions/scanning_v1/
addpath ~/matlab/dependencies/l1c/
addpath ~/matlab/dependencies/l1c/lib
addpath cs_simulations/splitBregmanROF_mex/
% initialize paths.
init_paths();
figbase = 20;


size_dir = '5microns';
fname = fullfile(PATHS.sysid, 'x-axis_sines_info_intsamps_zaxisFourierCoef_10-29-2018-01.mat');
models = load(fname);
G = -models.modelFit.G_zdir;
p = pole(G);
z = zero(G);
gg = zpk(z(end-1:end), p(1:2), 1, G.Ts);
hole_depth = (20);

chan_map = ChannelMap([1:5]);
exp_date = '4-26-2019'
% ----------------------- Load and Process CS-data -----------------------------
dat_root = PATHS.raster_image_data(size_dir, exp_date);

% ----------------------- Load and Process raster-data -------------------------
raster_files = {...
'raster_scan_512pix_5mic_01Hz_out_4-26-2019-01.csv',...
'raster_scan_512pix_5mic_04Hz_out_4-26-2019-01.csv',...
'raster_scan_512pix_5mic_05Hz_out_4-26-2019-01.csv',...
'raster_scan_512pix_5mic_08Hz_out_4-26-2019-01.csv',...
'raster_scan_512pix_5mic_10Hz_out_4-26-2019-01.csv',...
};
cs_files = {...
'cs-traj-512pix-10perc-500nm-5mic-01Hz-150prescan-notconnect_out_4-26-2019-01.csv',...
'cs-traj-512pix-10perc-500nm-5mic-02Hz-150prescan-notconnect_out_4-26-2019-01.csv',...
'cs-traj-512pix-10perc-500nm-5mic-04Hz-150prescan-notconnect_out_4-26-2019-01.csv',...
'cs-traj-512pix-10perc-500nm-5mic-05Hz-150prescan-notconnect_out_4-26-2019-01.csv',...
'cs-traj-512pix-10perc-500nm-5mic-08Hz-150prescan-notconnect_out_4-26-2019-01.csv',...
'cs-traj-512pix-15perc-500nm-5mic-01Hz-150prescan-notconnect_out_4-26-2019-01.csv',...
'cs-traj-512pix-15perc-500nm-5mic-02Hz-150prescan-notconnect_out_4-26-2019-01.csv',...
'cs-traj-512pix-15perc-500nm-5mic-04Hz-150prescan-notconnect_out_4-26-2019-01.csv',...
'cs-traj-512pix-15perc-500nm-5mic-05Hz-150prescan-notconnect_out_4-26-2019-01.csv',...
'cs-traj-512pix-15perc-500nm-5mic-08Hz-150prescan-notconnect_out_4-26-2019-01.csv',...
};


%%
rast_exps = {};
for k=1:length(raster_files)
  raster_paths = get_raster_paths(dat_root, raster_files{k});
  rast_exps{k} = RasterExp(raster_paths, 'load_full', true);
  if rast_exps{k}.time_total == 0
      rast_exps{k}.time_total = rast_exps{k}.samps_per_period*rast_exps{k}.npix*AFM.Ts;
  end
end
%

% dat = loadjson('/media/labserver/afm-cs/imaging/raster/5microns/parents/raster_scan_512pix_5mic_10Hz.json')


use_ze = false;
x1s = [63, 60, 57, 54, 52];
x2s = [488, 475, 473, 469, 417];
figbase = 10;
for k=1:length(rast_exps)
  rast_exps{k}.bin_raster_really_slow(@detrend, use_ze);
  
  pixmats_raw{k} = rast_exps{k}.pix_mat(1:end, 1:end);
%   rast_exps{k}.pix_mat_pinned = pixmats_raw{k};
  pixmat_ = pin_along_column(rast_exps{k}.pix_mat, x1s(k), x2s(k));
%   pixmat_ = pixmats_raw{k};
  rast_exps{k}.pix_mat_pinned = pixmat_ - mean(pixmat_(:));
  rast_exps{k}.pin_idx_s = [x1s(k), x2s(k)];
  stit = sprintf('(raster) %.2f Hz', rast_exps{k}.meta_in.raster_freq);
  plot_raster_data(rast_exps{k}.pix_mat_pinned, figbase*k, stit)
  % stit = sprintf('(raw) scan %d', k);
  % plot_raster_data(pixmats_raw{k}, (figbase-5)*k, stit)

end

for k=1:length(rast_exps)
%   rast_exps{k}.save()
end
%%
use_ze = false;
data_root = PATHS.cs_image_data(size_dir, exp_date);
cs_exps = {};
for k=1:length(cs_files)
  cs_paths = get_cs_paths(data_root, cs_files{k});
  % cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth, 'gg', gg);
  gg = @(u, idx_state_s) log_creep_detrend(u, idx_state_s);
  cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth, 'gg', gg,...
      'load_full', true);
  % cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth);
  cs_exps{k}.print_state_times();
  cs_exps{k}.sub_sample_frac()
end

%
bp = true;
recalc = false;
use_dct2 = true;
thresh = (20/7)*(1/1000)*20;
opts = l1qc_dct_opts('l1_tol', 0.00001, 'epsilon', 0.1);


% addpath ~/gradschool/publications/yufans/matlab-image-reconstruction-algorithms/'Mu path reconstruction with vertical penalty'/;

%%
for k=1:1%length(cs_exps)
  cs_exps{k}.process_cs_data(false, [], use_ze);
  fprintf('finished processing raw CS data...\n');
  fprintf('nperc=%.3f\n', sum(cs_exps{k}.pix_mask(:))/cs_exps{k}.npix^2);
  ht = cs_exps{k}.feature_height;
  if bp
   cs_exps{k}.solve_bp(recalc, use_dct2, opts);

   mu = 100;
   lamy = 0.1/10;

   gamma = mu/10;
   weight = 1.5; %0.001;%mu*5;
   sigma = 0.001;

   pix_idx = find(CsTools.pixmat2vec(cs_exps{k}.pix_mask) > 0.5);
   im = cs_exps{k}.Img_raw;

   [Ir] = bpvv_bregman(im,pix_idx,weight,mu,lamy, gamma, 10, sigma);
%%
    fprintf('Finished solving bp problem #%d\n', k);
    stit = sprintf('(CS) %.2f Hz equiv, \\%% %.2f sampling\nTotal time: %.2f',...
      cs_exps{k}.meta_in.tip_velocity/(cs_exps{k}.meta_in.width * 2),...
      cs_exps{k}.meta_in.actual_sub_samble_perc,...
      length(cs_exps{k}.x)*cs_exps{k}.Ts);
    
    f10=mkfig(1001 + 2*k, 6, 7.5); clf
    ax4 = subplot(3,1,[1,2]);
    ax4_2 =subplot(3,1,3);
    
    ImshowDataView.setup(f10);
    cb_exp =  @(event_obj)cs_exps{k}.dataview_callback(event_obj, ax4, ax4_2);
    bar_bp = mean(cs_exps{k}.Img_bp(:));
    cs_exps{k}.Img_bp = cs_exps{k}.Img_bp - bar_bp;
    cs_exps{k}.Img_raw = cs_exps{k}.Img_raw - bar_bp;
    
    im_tmp = cs_exps{k}.Img_bp - mean(cs_exps{k}.Img_bp(:));
    im_tmp = SplitBregmanROF(im_tmp, 100, 0.001);
    ImshowDataView.imshow(im_tmp,...
      [-thresh, thresh], ax4, ax4_2, cb_exp)
    title(ax4, stit)
  end
% keyboard
%   cs_exps{k}.save();
end
%%
for k=1:length(cs_exps)
  cs_exps{k}.save()
end
%%
Fig1 = mkfig(3000, 7, 9); clf;
ha1 = tight_subplot(4, 3, [0.025, 0.015], [.01, .02], [.02, .02], true);

Fig_err = mkfig(3001, 7, 9); clf
ha_err = tight_subplot(4, 3, [0.025, 0.015], [.01, .02], [.05, .02], true);
%
Fig_rows = mkfig(3002, 7, 4.5); clf
ha_row = tight_subplot(2, 1, [0.1, 0.015], [.1, .05], [.085, .02], false);
xlabel(ha_row(2), 'x-direction pixel', 'FontSize', 14)
ylabel(ha_row(1), 'height [nm]', 'FontSize', 14)
title(ha_row(1), 'raster', 'FontSize', 14)
title(ha_row(2), 'CS', 'FontSize', 14)
ylabel(ha_row(2), 'height [nm]', 'FontSize', 14)

grid(ha_row(1), 'on')
grid(ha_row(2), 'on')

mu = 100;
Img_filts = {};
mxs = [];
thresh = (20/7)*(1/1000)*20;
DRng = 2*thresh;

for k=1:length(cs_exps)
  Img_filts{k} = cs_exps{k}.Img_bp - min(cs_exps{k}.Img_bp(:));
  if ~isinf(mu)
    Img_filts{k} = SplitBregmanROFAn(Img_filts{k}, mu, 0.001);
  end
%   figure(101+k)
  mx = max(Img_filts{k}(:)) - min(Img_filts{k}(:));
  mxs = [mxs; mx];
%   mesh(Img_filts{k}), colormap('gray')
end

slice = 30:512-30;
master_idx = 1;
im_master = rast_exps{master_idx}.pix_mat_pinned - mean(rast_exps{master_idx}.pix_mat_pinned(:));
if ~isinf(mu)
  %im_master = SplitBregmanROF(im_master, mu, 0.001);
end
fprintf('---------------------------------------------------\n\n');
scan_metrics = {};
row_idx = 172;


excludes = [];
do_rows = [1, 3, 5];
err_idx = 1;
plt_idx = 1;
row_plt_idx = 1;
for k=1:length(rast_exps)
  imk = rast_exps{k}.pix_mat_pinned - mean(rast_exps{k}.pix_mat_pinned(:));
  if ~isinf(mu)
    %imk = SplitBregmanROF(imk, mu, 0.001);
  end
  imk_slice = imk(slice, slice);
  im1_ontok_fit = norm_align(imk_slice, im_master);
%   [im1_ontok_fit, imk_slice] = align_by_metric(im_master, imk, [], 'psnr');
  [psn_1k, ssm_1k] = ssim_psnr_norm(im1_ontok_fit, imk_slice, DRng);
  damage = rast_exps{k}.damage_metric();
  quality = rast_exps{k}.quality_metric();
  
  csm = ScanMetrics('ssim', ssm_1k, 'psnr', psn_1k,...
    'quality', quality, 'damage', damage,...
    'rate', rast_exps{k}.meta_in.raster_freq,...
    'time', rast_exps{k}.time_total,...
    'coverage', 100,...
    'type', 'raster');
    
  scan_metrics{end+1} = csm;
  
  if k ~= 1
      im_err = im1_ontok_fit - imk_slice;
      imagesc(ha_err(err_idx), im_err, [-thresh, thresh]);
      colormap(ha_err(err_idx), 'gray');
      stit_err = sprintf('raster (%.1f Hz)', csm.rate);
      title(ha_err(err_idx), stit_err);
      err_idx = err_idx + 1;
  end
  
  
  
  if k~=2
      imagesc(ha1(plt_idx), imk, [-thresh, thresh]);
      colormap(ha1(plt_idx), 'gray');
      stit = sprintf('raster (%.1f Hz)', csm.rate);
      title(ha1(plt_idx), stit);
      
      
      if ~isempty(intersect(k, do_rows))
          hold(ha1(plt_idx), 'on')
          hold(ha_row(1), 'on')
          
          plot(ha1(plt_idx), [1, 512], [row_idx, row_idx], 'r');
          hl_row = plot(ha_row(1), imk(row_idx, :)*AFM.volts2nm_z());
          hl_row.DisplayName = sprintf('%.1f Hz', csm.rate);
          hold(ha1(plt_idx), 'on')
          hold(ha_row(1), 'on')
          
      end
      plt_idx = plt_idx + 1;
  end

  fprintf('(raster %.1f Hz) psnr: %.4f, ssm: %.4f, damage: %.4f, time: %.4f\n',...
    rast_exps{k}.meta_in.raster_freq, psn_1k, ssm_1k, damage, csm.time);

end

remove_ticks(ha1)
remove_ticks(ha_err)
%

set(ha_row(1), 'XLim', [1, 512])
axes(ha_row(1))
leg1 = legend();
set(leg1, 'NumColumns', 3, 'FontSize', 11, 'Position', [0.5574 0.8968 0.4223 0.0500])


j = length(rast_exps);
fprintf('---------------------------------------------------\n');
cs_stats = {};
do_rows = [1, 4, 6, 10];
for k=1:length(cs_exps)

  imk = cs_exps{k}.Img_bp - mean(cs_exps{k}.Img_bp(:));
  if ~isinf(mu)
    imk = SplitBregmanROF(imk, mu, 0.001);
  end
  imk_slice = imk(slice, slice);
  im1_ontok_fit = norm_align(imk_slice, im_master);
  [psn_1k, ssm_1k] = ssim_psnr_norm(im1_ontok_fit, imk_slice, DRng);

  damage = cs_exps{k}.damage_metric();
  quality = cs_exps{k}.quality_metric();

  csm = ScanMetrics('ssim', ssm_1k, 'psnr', psn_1k,...
    'quality', quality, 'damage', damage,...
    'rate', cs_exps{k}.meta_in.tip_velocity/(cs_exps{2}.meta_in.width * 2),...
    'time', length(cs_exps{k}.x)*cs_exps{k}.Ts,...
    'coverage', cs_exps{k}.meta_in.actual_sub_samble_perc,...
    'type', 'CS');
  
  fprintf('(CS %.1f Hz, %% %.1f, %.1f) psnr: %.4f, ssm: %.4f, damage: %.4f, time: %.4f\n',...
    csm.rate, csm.coverage, cs_exps{k}.sub_sample_frac*100, psn_1k, ssm_1k, damage, csm.time);

    scan_metrics{end+1} = csm;
  [t_cycle_avg, t_connect] = cs_exps{k}.estimate_mpt_connect_savings();
  cs_stats{k, 1} = csm;
  cs_stats{k, 2} = t_cycle_avg;
  cs_stats{k, 3} = t_connect;
  

  if csm.rate ~= 4
      im_err = im1_ontok_fit - imk_slice;
      imagesc(ha_err(err_idx), im_err, [-thresh, thresh]);
      colormap(ha_err(err_idx), 'gray');
      title(ha_err(err_idx), sprintf('CS (%.1f \\%%, %.1f Hz)', csm.coverage, csm.rate));
      
      err_idx = err_idx + 1;
  end
  
  if csm.rate ~= 4
      imagesc(ha1(plt_idx), imk, [-thresh, thresh]);
      colormap(ha1(plt_idx), 'gray');
      title(ha1(plt_idx), sprintf('CS (%.1f \\%%, %.1f Hz)', csm.coverage, csm.rate));

      
      if ~isempty(intersect(k, do_rows))
          hold(ha1(plt_idx), 'on')
          hold(ha_row(2), 'on')
          
          plot(ha1(plt_idx), [1, 512], [row_idx, row_idx], 'r');
          hl_row = plot(ha_row(2), imk(row_idx, :)*AFM.volts2nm_z());
          hl_row.DisplayName = sprintf('%.1f Hz, %.1f \\%%', csm.rate, csm.coverage);
      end
      plt_idx = plt_idx + 1;
 end
  
  j=j+1;
end

remove_ticks(ha1)
remove_ticks(ha_err)

set(ha_row(2), 'XLim', [1, 512])
axes(ha_row(2))
leg2 = legend();
set(leg2, 'NumColumns', 4, 'FontSize', 11, 'Position', [0.1032 0.4253 0.8589 0.0500])

%%
csm_r = {};
csm_cs10 = {};
csm_cs15 = {};

tr = []
Hzr = [];
dmg_r = [];
psnr_r = [];
ssim_r = [];

t_cs10 = [];
Hz_cs10 = [];
dmg_cs10 = [];
psnr_cs10 = [];
ssim_cs10 = [];

t_cs15 = [];
Hz_cs15 = [];
dmg_cs15 = [];
psnr_cs15 = [];
ssim_cs15 = [];
for k = 1:length(scan_metrics)
    if strcmp(scan_metrics{k}.type, 'raster')
        csm_r{end+1} = scan_metrics{k};
        tr(end+1) = scan_metrics{k}.time;
        Hzr(end+1) = scan_metrics{k}.rate;
        dmg_r(end+1) = scan_metrics{k}.damage;
        ssim_r(end+1) = scan_metrics{k}.ssim;
        psnr_r(end+1) = scan_metrics{k}.psnr;
    elseif scan_metrics{k}.coverage < 14
        t_cs10(end+1) = scan_metrics{k}.time;
        Hz_cs10(end+1) = scan_metrics{k}.rate;
        dmg_cs10(end+1) = scan_metrics{k}.damage;
        ssim_cs10(end+1) = scan_metrics{k}.ssim;
        psnr_cs10(end+1) = scan_metrics{k}.psnr;
        
        csm_cs10{end+k} = scan_metrics{k};
    else
        t_cs15(end+1) = scan_metrics{k}.time;
        Hz_cs15(end+1) = scan_metrics{k}.rate;
        dmg_cs15(end+1) = scan_metrics{k}.damage;
        ssim_cs15(end+1) = scan_metrics{k}.ssim;
        psnr_cs15(end+1) = scan_metrics{k}.psnr;
        
        csm_cs15{end+k} = scan_metrics{k};
    end
end

% save('damage_data.mat', 'Hzr', 'tr', 'dmg_r', 't_cs10', 'Hz_cs10', 'dmg_cs10',...
%     't_cs15', 'Hz_cs15', 'dmg_cs15')
F = mkfig(3, 7, 3.5); clf
ha = tight_subplot(1, 2, .1, [0.12, 0.05], [0.1, 0.1]);

h3 = plot(ha(1), Hz_cs15, dmg_cs15, '*b');

xlabel(ha(1), 'scan rate [Hz]')
ylabel(ha(1), 'RDI')

hold(ha(1), 'on')
grid(ha(1), 'on')

h2 = plot(ha(1), Hz_cs10, dmg_cs10, 'or');
h1 = plot(ha(1), Hzr, dmg_r, 'x', 'MarkerSize', 8);


h4 = semilogx(ha(2), tr, dmg_r, 'x', 'MarkerSize', 8);
xlabel(ha(2), 'total time [s]')
ylabel(ha(2), 'RDI')

hold(ha(2), 'on')
grid(ha(2), 'on')

h4 = semilogx(ha(2), t_cs10, dmg_cs10, 'or');
h5 = semilogx(ha(2), t_cs15, dmg_cs15, '*b');


h1.DisplayName = 'raster';
h2.DisplayName = 'CS: 10\%';
h3.DisplayName = 'CS: 15\%';
legend([h1, h2, h3], 'location', 'northwest')

offsets = [50, 10, 10, 5, 5];
for k=1:5
    st = sprintf('%d Hz', Hzr(k));
    t1 = text(tr(k)+offsets(k), dmg_r(k), st, 'HorizontalAlignment', 'left')
end

offsets = [10, 10, 10, 5, 5];
for k=1:5
    st = sprintf('%d Hz', Hz_cs15(k));
    t1 = text(t_cs15(k)+offsets(k), dmg_cs15(k), st, 'HorizontalAlignment', 'left')
end


%%
save_fig(F, fullfile(PATHS.cs_final_fig(), 'cs_rast_damage'), false)

save_fig(F, fullfile(PATHS.defense_fig(), 'cs_rast_damage'), true)
%%
Hzr
dmg_r
A = [Hzr'*0+1, Hzr'];
bm = A\dmg_r'
m_r = bm(2)

A = [Hz_cs15'*0+1, Hz_cs15'];
bm = A\dmg_cs15'

m_cs = bm(2)


m_cs/m_r
pi/8

 % ---------------------------------------------------------------------- %
 %%
 
%From 1-Hz Raster-batch
ps = [23.25, 22.98, 24.99, 23.65, 22.51];
ssm = [0.7, 0.7, 0.73, 0.72, 0.71];

ps_mu = mean(ps)
ssm_mu = mean(ssm)

ps_std= std(ps)
ssm_std = std(ssm)


F = mkfig(4, 7, 3.5); clf
ha = tight_subplot(1, 2, .1, [0.12, 0.05], [0.1, 0.02], false);


xlabel(ha(1), 'total time [s]')
ylabel(ha(1), 'SSIM')
hold(ha(1), 'on')
grid(ha(1), 'on')

set(ha(1), 'XScale', 'log')

ylabel(ha(2), 'PSNR')
xlabel(ha(2), 'total time [s]')

hold(ha(2), 'on')
grid(ha(2), 'on')
set(ha(2), 'XScale', 'log')

h4 = semilogx(ha(1), tr(2:end), ssim_r(2:end), 'x', 'MarkerSize', 8);
h5 = semilogx(ha(1), t_cs10, ssim_cs10, 'or');
h6 = semilogx(ha(1), t_cs15, ssim_cs15, '*b');

low = [ssm_mu-ssm_std, ssm_mu-ssm_std];
up = [ssm_mu+ssm_std, ssm_mu+ssm_std]
hci = ciplot(low, up, [10, 200], 'b', ha(1));
alpha(hci, '0.25')

h7 = semilogx(ha(2), tr, psnr_r, 'x', 'MarkerSize', 8);

hold(ha(2), 'on')
grid(ha(2), 'on')

h8 = semilogx(ha(2), t_cs10, psnr_cs10, 'or');
h9 = semilogx(ha(2), t_cs15, psnr_cs15, '*b');

low = [ps_mu-ps_std, ps_mu-ps_std];
up = [ps_mu+ps_std, ps_mu+ps_std]
hci = ciplot(low, up, [10, 200], 'b', ha(2));
alpha(hci, '0.25')


h4.DisplayName = 'raster';
h5.DisplayName = 'CS: 10\%';
h6.DisplayName = 'CS: 15\%';
leg = legend([h4, h5, h6], 'Position', [0.3510 0.8104 0.1390 0.1370]);

save_fig(F, fullfile(PATHS.cs_final_fig(), 'cs_rast_time_vs_ssim_psnr'), false)

%%


F4 = mkfig(5, 7, 3.5); clf
ha4 = tight_subplot(1, 2, .1, [0.12, 0.05], [0.1, 0.1])


xlabel(ha4(1), 'scan rate [Hz]')
ylabel(ha4(1), 'SSIM')
hold(ha4(1), 'on')
grid(ha4(1), 'on')


xlabel(ha4(2), 'scan rate [Hz]')
ylabel(ha4(2), 'PSNR')
hold(ha4(2), 'on')
grid(ha4(2), 'on')

h10 = plot(ha4(1), Hz_cs15, ssim_cs15, '*b')
h11 = plot(ha4(1), Hz_cs10, ssim_cs10, 'or')
h12 = plot(ha4(1), Hzr(2:end), ssim_r(2:end), 'x', 'MarkerSize', 8)

h13 = plot(ha4(2), Hz_cs15, psnr_cs15, '*b')
h14 = plot(ha4(2), Hz_cs10, psnr_cs10, 'or')
h15 = plot(ha4(2), Hzr(2:end), psnr_r(2:end), 'x', 'MarkerSize', 8)



h1.DisplayName = 'raster';
h2.DisplayName = 'CS: 10\%';
h3.DisplayName = 'CS: 15\%';
legend([h1, h2, h3], 'location', 'northwest')
%%


% set(ha1(end), 'Visible', 'off')
% set(ha_err(end), 'Visible', 'off')


save_fig(Fig1, fullfile(PATHS.cs_final_fig(), 'cs_raster_images_4-26-2019'), false)
save_fig(Fig_err, fullfile(PATHS.cs_final_fig(), 'cs_raster_images_err_4-26-2019'),false)
save_fig(Fig_rows, fullfile(PATHS.cs_final_fig(), 'cs_raster_pixel_rows_4-26-2019'), false)

%%



%%
% skip the 0.5Hz scan
sm = scan_metrics(1:end)
S1 = scan_metrics_table(sm)

S2 = state_times_table(cs_exps)



fid = fopen(fullfile(PATHS.cs_final_table(), 'cs_raster_table_4-26-2019_muInf_dct2.tex'), 'w+');
fprintf(fid, '%s', S1);
fclose(fid);

fid = fopen(fullfile(PATHS.cs_final_table(), 'cs_state_times_table_4-26-2019_muInf_dct2.tex'), 'w+');
fprintf(fid, '%s', S2);
fclose(fid);


function S = scan_metrics_table(csm_s)
  S = 'type &  rate (Hz) & PSNR & SSIM & time [s] & damage\\';
  S = sprintf('%s\n\\toprule\n', S);
  
  for k=1:length(csm_s)
    csm = csm_s{k};
    if strcmp(csm.type, 'CS')
      description = sprintf('CS (%.2f~\\%%)', csm.coverage);
    else
      description = 'raster';
    end
    S = sprintf('%s%s & %.2f & %.2f & %.2f & %.1f & %.2f\\\\\n', S,...
      description, csm.rate, csm.psnr, csm.ssim, csm.time, csm.damage);
  end
end

function S = state_times_table(cs_exps)
  S = 'description & move & engage & tip-settle & scan & tip-up & total\\';
  S = sprintf('%s\n\\toprule\n', S);
  for k=1:length(cs_exps)
    cst = cs_exps{k}.get_state_times();
    description = sprintf('%.1f Hz, %.2f~\\%%',...
      cs_exps{k}.meta_in.tip_velocity/(cs_exps{2}.meta_in.width * 2),...
      cs_exps{k}.meta_in.actual_sub_samble_perc);
  
    S = sprintf('%s%s &%.2f & %.2f & %.2f & %.2f & %.2f & %.2f\\\\\n', S,...
      description, cst.move, cst.tdown, cst.tsettle, cst.scan,...
      cst.tup, cst.total);
  end
end




%%
function plot_raster_data(pixmat2, figbase, stit, plot_mesh)
  if nargin < 4
    plot_mesh = false;
  end
  thresh = (20/7)*(1/1000)*20;
  pixmat2 = pixmat2 - mean(pixmat2(:));
%   F10 = figure(figbase+2); clf
  F10 = mkfig(figbase+2, 6, 7.5, false);
  ax1 = subplot(3,1,[1,2]);
  ax2 =subplot(3,1,3);
    

  lo = min(min(pixmat2));
  hi = max(max(pixmat2));
  
  
  imshow_dataview(flipud(pixmat2 - mean(pixmat2(:))), [-thresh, thresh], ax1, ax2)
  try
    grid(ax1, 'on')
  catch
    keyboard
  end
  
  colormap(ax1, 'gray')
  grid(ax2, 'on')
  ax1.GridAlpha = 1;
  ax2.GridAlpha = 1;
  title(ax1, stit)
  title(ax2, stit)
  if plot_mesh
    figure(figbase+4)
    mesh(pixmat2)
    xlabel('x')
    ylabel('y')
    title(stit)
  end
end
