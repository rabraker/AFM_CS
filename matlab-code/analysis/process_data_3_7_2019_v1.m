clc
clear
%
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
exp_date = '3-7-2019'
% ----------------------- Load and Process CS-data -----------------------------
dat_root = PATHS.raster_image_data(size_dir, exp_date);

% ----------------------- Load and Process raster-data -------------------------
raster_files = {...
'raster_scan_512pix_5mic_5.00e-01Hz_out_3-7-2019-01.csv',...
'raster_scan_512pix_5mic_01Hz_out_3-7-2019-01.csv',...
'raster_scan_512pix_5mic_02Hz_out_3-7-2019-01.csv',...
'raster_scan_512pix_5mic_04Hz_out_3-7-2019-01.csv',...
};

cs_files = {...
'cs-traj-512pix-10perc-500nm-5mic-01Hz-250prescan-isconnect_out_3-7-2019-01.csv',...
'cs-traj-512pix-10perc-500nm-5mic-02Hz-250prescan-isconnect_out_3-7-2019-01.csv',...
'cs-traj-512pix-7perc-500nm-5mic-02Hz-250prescan-isconnect_out_3-7-2019-01.csv'  
};


pixmaps = cell(length(raster_files)+length(cs_files), 1);
for k=1:length(raster_files)
  raster_paths = get_raster_paths(dat_root, raster_files{k});
  rast_exps{k} = RasterExp(raster_paths);
end

data_root = PATHS.cs_image_data(size_dir, exp_date);

for k=1:length(cs_files)
  cs_paths = get_cs_paths(data_root, cs_files{k});
  cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth);% 'gg', gg);
  % cs_exp1 = CsExp(cs_paths, 'feature_height', hole_depth);% 'gg', gg);
  % cs_exp1.uz = detrend(cs_exp1.uz);
  % cs_exp1.print_state_times();
end

npix = rast_exps{1}.npix;
Ts = rast_exps{1}.Ts;
width = rast_exps{1}.width;
stit = sprintf('Scan %d', k)
%%
clc
N_damage = 100;
for k=1:length(rast_exps)
  ze{k} = rast_exps{k}.ze;
  ze{k} = resample(ze{k}, 1, 2); % half the sample rate
  ze{k} = ze{k} - (-0.3); %mean(ze1);
  ze_sort{k} = sort(ze{k}(ze{k}>0), 'descend');
  %   sum(ze_sort{k}(1:N_damage))/N_damage
  sum(ze{k}(ze{k}>0)) / length(ze{k}(ze{k}>0))
end

for k=1:length(cs_exps)
  ze_cs_scan = [];
  for j=1:length(cs_exps{k}.idx_state_s.scan)
    ze_cs_scan = [ze_cs_scan; cs_exps{k}.ze(cs_exps{k}.idx_state_s.scan{j})];
  end
  ze_cs{k} = ze_cs_scan;
  ze_cs{k} = ze_cs{k} - (-0.3); %- mean(ze_cs);
  zecs_sort{k} = sort(ze_cs{k}(ze_cs{k}>0), 'descend');
  %   sum(zecs_sort{k}(1:N_damage))/N_damage
  sum(ze_cs{k}(ze_cs{k}>0))/length(ze_cs{k}(ze_cs{k}>0))
end


%%
clc

%%
ze1_pos = ze1(ze1>0);
ze2_pos = ze2(ze2>0);


sum(ze1_pos)/length(ze1_pos)
sum(ze2_pos)/length(ze2_pos)
ze_cs = cs_exp1.ze;
ze_cs = ze_cs - (-0.3); %- mean(ze_cs);
ze_cs_pos = ze_cs(ze_cs>0);
sum(ze_cs_pos)/length(ze_cs_pos)
%%
clc
x1s = [47, 52, 54, 57];
x2s = [459, 468, 471 473];
figbase = 10;
for k=1:length(rast_exps)
  rast_exp2 = copy(rast_exps{k});
%   uz = detrend(rast_exp2.uz);
  rast_exp2.bin_raster_really_slow(@detrend);
  
  pixmats_raw{k} = rast_exp2.pix_mat(1:end-1, 1:end-1);
%   pixmats{k} = pixmats_raw{k};
  pixmats{k} = pin_along_column(pixmats_raw{k}, x1s(k), x2s(k));
  pixmats{k} = pixmats{k} - mean(pixmats{k}(:));
  stit = sprintf('(raster) %.2f Hz', rast_exps{k}.meta_data.raster_freq);
  plot_raster_data(pixmats{k}, figbase*k, stit)
%   stit = sprintf('(raw) scan %d', k);
%   plot_raster_data(pixmats_raw{k}, (figbase-5)*k, stit)
  
end

%%
[~, axs] = make_cs_traj_figs(figbase, 2);
cs_exps{1}.plot_all_cycles(axs{1:2});


%%



%%
bp = true;
for k=1:length(cs_exps)
  cs_exps{k}.process_cs_data(false, []);
  fprintf('finished processing raw CS data...\n');
  fprintf('nperc=%.3f\n', sum(cs_exps{k}.pix_mask(:))/cs_exps{k}.npix^2);
  ht = cs_exps{k}.feature_height;
  if bp
    cs_exps{k}.solve_bp(true, true);
    
    fprintf('Finished solving bp problem #%d\n', k);
    stit = sprintf('(CS) %.2f Hz equiv, %% %.2f sampling\nTotal time: %.2f',...
      cs_exps{k}.meta_in.tip_velocity/(cs_exps{2}.meta_in.width * 2),...
      cs_exps{k}.meta_in.actual_sub_samble_perc,...
      length(cs_exps{k}.x)*cs_exps{k}.Ts);
    
    f10=figure(1001 + 2*k); clf
    ax4 = subplot(3,1,[1,2]);
    ax4_2 =subplot(3,1,3);
    
    ImshowDataView.setup(f10);
    cb_exp =  @(event_obj)cs_exps{k}.dataview_callback(event_obj, ax4, ax4_2);
    ImshowDataView.imshow(cs_exps{k}.Img_bp, [], ax4, ax4_2, cb_exp)
    title(ax4, stit)
  end

end
%%
Img_filts = {};
for k=1:length(cs_exps)
Img_filts{k} = cs_exps{k}.Img_bp - min(cs_exps{k}.Img_bp(:));

Img_filts{k} = SplitBregmanROF(Img_filts{k}, 30, 0.001);
figure(101+k)
mesh(Img_filts{k}), colormap('gray')
end
%%
slice = 20:512-20;

im_master = pixmats{1} - mean(pixmats{1}(:));
im_master = SplitBregmanROF(im_master, 30, 0.001);
fprintf('---------------------------------------------------\n');
for k=2:length(rast_exps)
  imk = pixmats{k} - mean(pixmats{k}(:));
  imk = SplitBregmanROF(imk, 30, 0.001);
  imk_slice = imk(slice, slice);
  im1_ontok_fit = norm_align(imk_slice, im_master);
  [psn_1k, ssm_1k] = ssim_psnr_norm(im1_ontok_fit, imk_slice);

fprintf('(raster %.1f Hz) psnr: %.4f, ssm: %.4f\n',...
  rast_exps{k}.meta_data.raster_freq, psn_1k, ssm_1k);

end
%%
fprintf('---------------------------------------------------\n');
for k=1:length(cs_exps)
  imk = cs_exps{k}.Img_bp - mean(cs_exps{k}.Img_bp(:));
  imk = SplitBregmanROF(imk, 30, 0.001);
  imk_slice = imk(slice, slice);
  im1_ontok_fit = norm_align(imk_slice, im_master);
  [psn_1k, ssm_1k] = ssim_psnr_norm(im1_ontok_fit, imk_slice);

stit = sprintf('(CS: %.2f Hz, %% %.2f sampl., time: %.2f)',...
      cs_exps{k}.meta_in.tip_velocity/(cs_exps{2}.meta_in.width * 2),...
      cs_exps{k}.meta_in.actual_sub_samble_perc,...
      length(cs_exps{k}.x)*cs_exps{k}.Ts);
fprintf('(CS %s) psnr: %.4f, ssm: %.4f\n',...
  stit, psn_1k, ssm_1k);

end
%%


%%
function plot_raster_data(pixmat2, figbase, stit)

  thresh = (20/7)*(1/1000)*20;
  pixmat2 = pixmat2 - mean(pixmat2(:));
  F10 = figure(figbase+2); clf
  ax1 = subplot(3,1,[1,2])
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
  
  figure(figbase+4)
  mesh(pixmat2)
  xlabel('x')
  ylabel('y')
  title(stit)
end
