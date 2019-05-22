clc
clear
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


exp_date = '5-22-2019'
% ----------------------- Load and Process CS-data -----------------------------
dat_root = PATHS.raster_image_data(size_dir, exp_date);

% ----------------------- Load and Process raster-data -------------------------
raster_files = {...
'raster_scan_64pix_5mic_10Hz_fouriery0_out_5-22-2019-01.csv',...
'raster_scan_64pix_5mic_10Hz_y0_out_5-22-2019-01.csv',...
'raster_scan_64pix_5mic_15Hz_fouriery0_out_5-22-2019-01.csv',...
'raster_scan_64pix_5mic_15Hz_y0_out_5-22-2019-01.csv',...
};



pixmaps = cell(length(raster_files), 1);
for k=1:length(raster_files)
  raster_paths = get_raster_paths(dat_root, raster_files{k});
  rast_exps{k} = RasterExp(raster_paths);
  
  rast_exps{k}.uz  = rast_exps{k}.uz*AFM.volts2nm_z;
  rast_exps{k}.xref = rast_exps{k}.xref * AFM.volts2mic_xy;
end

name_s = {'CR 10hz fourier', 'CR 10hz w/o fourier', 'CR 15 hz w/ fourier', 'CR 15 Hz w/o fourier'};
%%
close all
figure(1); clf
ax1 = gca()
for k=1:2 %length(rast_exps)
    plot(rast_exps{k}.x*512)
    hold on
end



% figure(2); clf
% ax2 = gca()
% for k=1:length(rast_exps)
%     plot(rast_exps{k}.uz)
%     hold on
% end
% 
% figure(3); clf
% ax3 = gca()
% for k=1:length(rast_exps)
%     plot(rast_exps{k}.ze)
%     hold on
% end

linkaxes([ax1, ax2, ax3], 'x')

%%
clc
npix_x = 512;
rast_exps{1}.bin_raster_really_slow(@detrend, [], npix_x);
rast_exps{2}.bin_raster_really_slow(@detrend, [], npix_x);
rast_exps{3}.bin_raster_really_slow(@detrend, [], npix_x);
rast_exps{4}.bin_raster_really_slow(@detrend, [], npix_x);

figure(10), imagesc(rast_exps{1}.pix_mat)
figure(11), imagesc(rast_exps{2}.pix_mat)
figure(12), imagesc(rast_exps{3}.pix_mat)
figure(13), imagesc(rast_exps{4}.pix_mat)


%%

figure(5); clf
hands = gobjects(4, 1);
for k=1:length(rast_exps)
    fig = figure(k);
    [uzk1, freqs] = psd_by_line(rast_exps{k}, @detrend, fig);
    figure(5)
    hands(k) = semilogx(freqs, 10*log10(abs(uzk1)));
    hands(k).DisplayName = name_s{k};
    hold on

end

N = length(rast_exps{1}.xref) ;

[Pxx, freqs] = periodogram(rast_exps{1}.xref, [], N, 1/AFM.Ts);

idx = find(abs(Pxx) > 1e-12);
semilogx(freqs(idx), 10*log10(abs(Pxx(idx))), 'bo')

% -----------------------------
N = length(rast_exps{2}.xref);
[Pxx, freqs] = periodogram(rast_exps{2}.xref, [], N, 1/AFM.Ts);
idx = find(abs(Pxx) > 1e-12);
semilogx(freqs(idx), 10*log10(abs(Pxx(idx))), 'kx')

% -----------------------------
N = length(rast_exps{3}.xref);
[Pxx, freqs] = periodogram(rast_exps{3}.xref, [], N, 1/AFM.Ts);

idx = find(abs(Pxx) > 1e-12);
semilogx(freqs(idx), 10*log10(abs(Pxx(idx))), 'ro')

% -----------------------------
N = length(rast_exps{4}.xref);
[Pxx, freqs] = periodogram(rast_exps{4}.xref, [], N, 1/AFM.Ts);

idx = find(abs(Pxx) > 1e-12);
semilogx(freqs(idx), 10*log10(abs(Pxx(idx))), 'rx')


legend(hands)
grid on

% [uzk2, freqs] = psd_by_line(rast_exps{2}, @detrend);

% figure(5)
% semilogx(freqs, 10*log10(abs(uzk1)))
% hold on
% semilogx(freqs, 10*log10(abs(uzk2)))
% 
% figure(3)
% mesh(rast_exps{1}.pix_mat)
% figure(4)
% mesh(rast_exps{2}.pix_mat)



function [uzk, freqs] = psd_by_line(rast_exp, line_detrender, fig)

%   xpix = rast_exp.npix; 
%   ypix = rast_exp.npix; % TODO: make this work with rectangular image.
    
  % Get the indeces corresponding to trace data only.
%   trace_inds = rast_exp.get_trace_indeces();
%   xdat_trace = rast_exp.x(trace_inds);
%   ydat_trace = rast_exp.y(trace_inds);
%   udat_trace = rast_exp.uz(trace_inds);
  
  offset = rast_exp.find_lag();
  samps_per_line = rast_exp.samps_per_line
  uzk = 0;
  for k=0:rast_exp.npix_y-1
    ind_start = ceil(k*(samps_per_line)+1);
    ind_end = floor((k+1)*(samps_per_line));
    idx = (ind_start:ind_end) + offset;
    
    uz_ = line_detrender(rast_exp.uz(idx));
    [uzk_, freqs] = power_spectrum_local(uz_);
    uzk = uzk + uzk_;
    
    figure(fig)
    subplot(2,1,1)
    plot(rast_exp.x(idx))
    hold on
    subplot(2,1,2)
    plot(uz_)
    hold on
%     keyboard
  end
  
  uzk = uzk / rast_exp.npix_y;
end

function [Pxx, freqs] = power_spectrum_local(X)

  N = length(X);

  [Pxx, freqs] = periodogram(detrend(X), [], N, 1/AFM.Ts);
end 