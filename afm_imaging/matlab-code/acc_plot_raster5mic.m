clc
close all
addpath('functions')
% dat_root = 'C:\Users\arnold\Documents\labview\afm_imaging\data';
dat_root = getdataroot();
% data_root = 'C:\Users\arnold\Documents\labview\afm_imaging\data'
% datmat = csvread(fullfile(data_root, 'raster_8-1-2017_v3.csv'));

% -------- Constants --------
Ts = 40e-6;
width = 5;
nperiods = 256;
npix = nperiods;
volts2micron = 50/10;
micron2pix = npix/width;
volts2pix = volts2micron * micron2pix;

% ------------------- 1Hz 
% dat_name_s{1} = 'raster_scan_5mic_1Hz_out_8-26-2017-01-full.csv';
% parent_name_s{1} = 'raster_scan_5mic_1Hz.csv';
dat_name_s{1} =  'raster_scan_5mic_1Hz_out_9-4-2017-01-full.csv';
parent_name_s{1} = 'raster_scan_5mic_1Hz.csv';
% ------------------- 5Hz 
dat_name_s{2} = 'raster_scan_5mic_5Hz_out_8-26-2017-01-full.csv';
parent_name_s{2} = 'raster_scan_5mic_5Hz.csv';
% ------------------- 10Hz 
dat_name_s{3} = 'raster_scan_5mic_10Hz_out_8-26-2017-01-full.csv';
parent_name_s{3} = 'raster_scan_5mic_10Hz.csv';

freqs = [1,5, 10];
for i=1:length(freqs)
    dat_path_s{i} = fullfile(dat_root,'data/raster', dat_name_s{i});
    parent_path_s{i} = fullfile(dat_root, 'data/raster', parent_name_s{i});
end
raster_dat_s = [];
for i=1:length(freqs)
    rdat.datmat = csvread(dat_path_s{i});
    rdat.parent_dat = csvread(parent_path_s{i});

    xyref = reshape(rdat.parent_dat, 2, [])';
    rdat.xref = xyref(:,1);
    rdat.Ts = Ts;
    rdat.npix = npix;
    rdat.volts2pix = volts2pix;
    rdat.width = width;
    
    rdat.samps_per_period = size(rdat.parent_dat,1)/2; % twice as many in here for x & y.
    rdat.samps_per_line = rdat.samps_per_period/2;
    rdat.datmat = rdat.datmat([1:rdat.npix*rdat.samps_per_period], :);

    [pixmat, pixelifsampled] = bin_raster_slow(rdat.datmat(:,[1,2,4]),...
                                rdat.npix, rdat.samps_per_period, rdat.volts2pix);
    I_fit = detrend_sampled_plane(pixmat, pixelifsampled);
    I_fit = (I_fit - min(min(I_fit))).*pixelifsampled;
    rdat.pixmat = pixmat;
    rdat.I_fit = I_fit;
    rdat.pixelifsampled = pixelifsampled;
    rdat.freq = freqs(i);
    raster_dat_s = [raster_dat_s, rdat];

end



% [pixmat2, pixelifsampled] = bin_raster_slow(datmat1hz(:,[1,2,4]), npix, samps_per_period, volts2pix);


%%
clc
F1 = figure(1); clf
% axes('FontName','Times New Roman') % Set axis font style
% box('on'); % Define box around whole figure


figwidth = 7.5;
figheight = 4.5;
set(F1, 'Units', 'Inches', 'Position', [3,3, figwidth, figheight],...
    'PaperUnits', 'Inches', 'PaperSize', [figwidth, figheight])
set(F1, 'Color', 'w');


% subplot(4,3,[1,4])
xpad = .03;
wd = .28;

lft1 = .073; % start of first col;

lft2 = lft1+wd+xpad;
lft3 = lft2+wd+xpad;


ht_im = .37; % image heights
bt_im = .58;
bt_r2 = .26;
bt_r3 = .051;

ht_r2 = .2;
ht_r3 = .15;
ax1 = axes('Position', [lft1, bt_im, wd, ht_im]);

% subplot(4,3,[2,5])
ax2 = axes('Position', [lft2, bt_im, wd, ht_im]);

% subplot(4,3, [3,6])
ax3 = axes('Position', [lft3, bt_im, wd, ht_im]);

% row 2
ax4=axes('Position', [lft1, bt_r2, wd, ht_r2]);
ax5 = axes('Position', [lft2, bt_r2, wd, ht_r2]);
set(ax5, 'YTickLabel', [])
ax6 = axes('Position', [lft3, bt_r2, wd, ht_r2]);
set(ax6, 'YTickLabel', [])

ax7 = axes('Position', [lft1, bt_r3, wd, ht_r3]);
ax8 = axes('Position', [lft2, bt_r3, wd, ht_r3]);
set(ax8, 'YTickLabel', [])
ax9 = axes('Position', [lft3, bt_r3, wd, ht_r3]);
set(ax9, 'YTickLabel', [])
%
% ---------------- Start Plotting -----------------------------
% ------------------ Col 1: 1Hz -------------------------------
axr1 = [ax1, ax2, ax3]';
for iter = 1:length(axr1)
    ax_iter = axr1(iter);
    F1.CurrentAxes = ax_iter;
    
    lo = min(min(raster_dat_s(iter).I_fit));
    hi = max(max(raster_dat_s(iter).I_fit));

    % ax1.XTick = [0 2 4]
    % ax1.YTick = [0 2 4]
    xdata = [0, width];
    ydata = [0, width];

    imshow(raster_dat_s(iter).I_fit, [lo, hi],...
        'XData', xdata, 'YData', ydata, 'Parent', ax_iter)
    % imshow(I_fit, [lo, hi], 'Parent', ax1)


    axis('on')
    % ax1.YTickLabel = flipud({ax1.YTickLabel{2:end}})

    xlabel('x-dir [$\mu$m]', 'interpreter', 'latex')
    ylabel('y-dir [$\mu$m]',  'interpreter', 'latex')

    stit = sprintf('%d~Hz', raster_dat_s(iter).freq)
    title(stit, 'interpreter', 'latex')
end
%
% ---------------------------- second row, raster tracking --------------
axr2 = [ax4, ax5, ax6]
% subplot(4,3,7); hold on
for iter=1:length(axr2)
    ax_iter = axr2(iter);
    F1.CurrentAxes = ax_iter;
    hold on
    np = 2;
    xx = raster_dat_s(iter).datmat(:,1);
    x_np = xx(1:np*length(raster_dat_s(iter).xref));
    x_np = x_np - min(x_np);
    t = [0:1:length(raster_dat_s(iter).xref)*np-1]'*Ts;
    xref_np = repmat(raster_dat_s(iter).xref, np, 1);

    p1 = plot(ax_iter, t, xref_np*volts2micron, 'k');
    p1.DisplayName = '$x_{ref}$';
    p2 = plot(ax_iter, t, x_np*volts2micron, '--r');
    p2.DisplayName = '$x(k)$';

    ylm = ylim;
    ylim([0, ylm(2)+.1*ylm(2)])
    if iter>1
%         ylabel('x-dir~[$\mu$ m]', 'interpreter', 'latex')
    end
    xlabel('time [s]', 'interpreter', 'latex')
    % leg1 = legend([p1, p2]);
    % set(leg1, 'FontSize', 14, 'interpreter', 'latex', 'orientation', 'horizontal')
    % leg1.Position=[0.6436    0.8590    0.2611    0.0640];

end

axr3 = [ax7, ax8, ax9];
nn_s = [11,11,21]
for iter = 1:length(axr3)
    ax_iter = axr3(iter);
    F1.CurrentAxes = ax_iter

    samps_per_period = raster_dat_s(iter).samps_per_period;
    samps_per_line= raster_dat_s(iter).samps_per_line;
    uz = raster_dat_s(iter).datmat(:,4);
    t2 = [nn*samps_per_period:nn*samps_per_period+samps_per_line-1]'*Ts;
    nn = nn_s(iter);
    plot(ax_iter, t2, uz(nn*samps_per_period:nn*samps_per_period+samps_per_line-1))
    xlim([t2(1), t2(end)])
    if iter > 1
%         ylabel('$u_z$', 'interpreter', 'latex')
        set(ax_iter, 'YTickLabel', [])
    end
end

%%
fig_path = fullfile(getfigroot, 'rasterscans.pdf');
export_fig(F1, fig_path, '-q101')

