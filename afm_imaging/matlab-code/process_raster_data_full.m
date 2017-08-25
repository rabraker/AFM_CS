clc
close all
% dat_root = '/home/arnold/gradschool/afm-cs/afm_imaging/data';
dat_root = 'C:\Users\arnold\Documents\labview\afm_imaging\data';
% data_root = 'C:\Users\arnold\Documents\labview\afm_imaging\data'
% datmat = csvread(fullfile(data_root, 'raster_8-1-2017_v3.csv'));
dat_name = 'raster-scan-8-24-2017_7Hz-full.csv';
parent_name = 'raster_traj_singleperiod_7Hz.csv';

dat_path = fullfile(dat_root, dat_name);
parent_path = fullfile(dat_root, parent_name);
datmat = csvread(dat_path);
parent_dat = csvread(parent_path);

xyref = reshape(parent_dat', 2, [])';
xref = xyref(:,1);
%%
Ts = 40e-6;
samps_per_period = size(parent_dat,1)/2; % twice as many in here for x & y.
samps_per_line = samps_per_period/2;

nperiods = 256;
pix = nperiods;
datmat = datmat([1:nperiods*samps_per_period], :);


trace_inds = get_trace_indeces(nperiods, samps_per_period);
xdat = reshape(datmat(trace_inds,1), [], nperiods);
ydat = reshape(datmat(trace_inds,2), [], nperiods);
udat = reshape(datmat(trace_inds,4), [], nperiods);

%
udat = datmat(:,4);
pixmat = bin_raster_fast(udat, pix, samps_per_period);

lo = min(min(pixmat));
hi = max(max(pixmat));
figure(1)
imshow(pixmat, [lo, hi])

I_fit = detrend_plane(pixmat);

figure(3)
lo = min(min(I_fit));
hi = max(max(I_fit));
imshow(I_fit, [lo, hi])
%%
% visualize tracking error.
np = 3;
xx = datmat(:,1);
x_np = xx(1:np*length(xref));
x_np = x_np - min(x_np);
figure(200); clf; hold on
t = [0:1:length(xref)*np-1]'*Ts;
xref_np = repmat(xref, np, 1);
p1 = plot(t, xref_np*volts2micron);
p1.DisplayName = '$x_{ref}$';
p2 = plot(t, x_np*volts2micron);
p2.DisplayName = '$x(k)$';

ylm = ylim;
ylim([0, ylm(2)+.1*ylm(2)])
ylabel('x-dir [\mu m]')
xlabel('time [s]')
leg1 = legend([p1, p2]);
set(leg1, 'FontSize', 14, 'interpreter', 'latex', 'orientation', 'horizontal')
leg1.Position=[0.6436    0.8590    0.2611    0.0640];

%% Now try by actually using the x-y measured data;



width = 5;
volts2micron = 50/10;
micron2pix = pix/width;
volts2pix = volts2micron * micron2pix;


pixmat2 = bin_raster_slow(datmat(:,[1,2,4]), pix, samps_per_period, volts2pix);



% pixmat2 = detrend_plane(pixmat2);
figure(5)
lo = min(min(pixmat2));
hi = max(max(pixmat2));
imshow(pixmat2, [lo, hi])



