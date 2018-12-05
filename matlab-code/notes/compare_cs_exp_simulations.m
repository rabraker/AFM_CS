
npix = 512;
img_mat = ones(npix,npix);

% Holes are on a  500 nm pitch, and we scan 5 microns, so 10 holes.
% Thus, assuming the holes are half the pitch, each hole is 250 nm and the
% space between is 250 nm. Thus, 1/20th of the image width.
%
hole_width = npix/20
pitch = npix/10;

m_pitch = ceil(pitch);
n_hole = ceil(hole_width);

% First lets, create a prototype hole. Then we will insert it into the master
% image periodically.
sq = ones(n_hole,n_hole);
% circle radius
r = n_hole/2;
% center
xc = n_hole/2;
yc = n_hole/2;


for x_pt=1:n_hole
  for y_pt=1:n_hole
    
    rad_pt = sqrt( (x_pt - xc)^2 + (y_pt - yc)^2 );
    
    if rad_pt <= r
      sq(y_pt, x_pt) = 0; % make it black
    end
  end
end


figure(14);

imshow(sq, [0, 1])


x_start = 20;
y_start = 20;

figure(14); clf
imshow(img_mat, [0, 1])
x_lft = x_start;
while x_lft+n_hole < npix +n_hole
  y_tp = y_start;
  while y_tp+n_hole <= npix + n_hole
    img_mat(y_tp:y_tp+n_hole-1, x_lft:x_lft+n_hole-1) = sq;
    
    y_tp = y_tp + m_pitch;
%     imshow(img_mat, [0, 1])
%     drawnow();
%     pause
  end
  
  x_lft = x_lft + m_pitch;
  
end
img_mat = img_mat(1:npix, 1:npix);
figure(14); clf
imshow(img_mat, [0, 1])

%%
% Now, load some cs data. Use the pix_mask, and reconstruct.

init_paths();
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
cs_exp1.solve_basis_pursuit();

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
%%
cs_exp3.solve_basis_pursuit(true);
%%

cs_exp1.save()
cs_exp2.save()
cs_exp3.save()

thresh = (20/7)*(1/1000)*20;
cs_sim_path = '~/matlab/afm-cs/matlab-code/notes/data/cs_sim_CS20NG.mat';

if 0
  cs_sim = CsSim(img_mat*thresh, cs_exp1.pix_mask);
  cs_sim.solve_bp();
  cs_sim.Img_bp = cs_sim.Img_bp - mean(cs_sim.Img_bp(:));
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
%
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
F11 = figure(11); clf
set(F11, 'Position', [-880 41 877 956]);
ypad = 0.01;
xpad1 = 0.05;
% Simulation axes
wd = 0.45;
wd_im = 0.4;

bt2 = 0.14;
bt1 = 0.04;
bt3 = bt2 + wd_im - 0.0;


lft_im_a = 0.085;
lft_im_b = 0.58;
lft1a = 0.065;
lft1b = lft1a + xpad1 + wd;

ht1 = 0.15;
ax0a = axes('Position', [lft_im_a, bt3, wd_im, 0.5154]);
ax0b = axes('Position', [lft_im_b, bt3, wd_im, 0.5154]);

ax1 = axes('Position', [lft_im_a, bt2, wd_im, 0.5154]);
ax2 = axes('Position', [lft1a bt1 wd, ht1]);

% Experimental exes
ax3 = axes('Position', [lft_im_b, bt2, wd_im ,0.5154]);
ax4 = axes('Position', [lft1b, bt1, wd, ht1]); 


% Simulation master
ImshowDataView.imshow(img_mat, [0,1], ax0a, ax2);
 
% SImulartion CS
cb_sim =  @(event_obj)cs_sim.dataview_callback(event_obj, ax1, ax2);
ImshowDataView.imshow(cs_sim.Img_bp,   [th_mn_sim, th_mx_sim], ax1, ax2, cb_sim)
title(ax1, 'Simulated reconstruction') 
% Experiment.
cb_exp2 =  @(event_obj)cs_exp2.dataview_callback(event_obj, ax3, ax4);
ImshowDataView.imshow(cs_exp2.Img_bp, [tr_mn*.8, tr_mx], ax3, ax4, cb_exp2)
stit = sprintf('Experimental reconstruction. Total time = %.1f', cs_exp2.time_total);
title(ax3, stit)

xlabel(ax2, 'pixel')
title(ax0a, 'Pure simulation')
 ImshowDataView.setup(gcf)


ImshowDataView.imshow(rast_exp.pix_mat, [-thresh, thresh], ax0b, ax4);
title(ax0b, '1 Hz raster scan')


% F11 = mkfig(7, 4, 6, true); clf
% ax1 = axes('Position', [0.11438 0.4513 0.775 0.51537]);
% ax2 = axes('Position', [0.13 0.11 0.775 0.27194]);
% ax2.Color = fig_color;
%%
save_fig(F11, 'notes/figures/cs_compare_sim_sim', false)

