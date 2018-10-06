clc
clear
close all
addpath('functions')

% dat_root = getdataroot();
data_root = '/media/labserver/acc2018-data';
% -------- Constants --------
Ts = 40e-6;
width = 20;
nperiods = 512;
npix = nperiods;
volts2micron = 50/10;
micron2pix = npix/width;
volts2pix = volts2micron * micron2pix;

% ---------------------- Raster Data ----------------------------------- %
dat_name_s{1} = 'raster_scan_512pix_20mic_1.00e-01Hz_out_9-18-2017-02.csv';
% dat_name_s{2} = 'raster_scan_512pix_20mic_2.50e-01Hz_out_9-18-2017-01.csv';
dat_name_s{2} = 'raster_scan_512pix_20mic_01Hz_out_9-18-2017-01.csv';
dat_name_s{3} = 'raster_scan_512pix_20mic_2.50e+00Hz_out_9-18-2017-01.csv';
dat_name_s{4} = 'raster_scan_512pix_20mic_05Hz_out_9-20-2017-01.csv';
sub_fold = '20microns';



raster_dat_s = [];
max_s = zeros(1,4);
min_s = zeros(1,4);
for i=1:length(dat_name_s)
    data_path= fullfile(data_root,'raster', sub_fold, dat_name_s{i});
    parent_name = get_parent_name(dat_name_s{i}, '_out_');
    parent_path = fullfile(data_root, 'raster', sub_fold, parent_name);
    parent_meta_path = strrep(parent_path, '.csv', '.mat');
    meta_path = strrep(data_path, '.csv', '-meta.mat');
    mat_path = strrep(data_path, '.csv', '.mat');
<<<<<<< HEAD

    fprintf('File: %s\n', dat_name_s{i});

=======
    
    fprintf('File: %s\n', dat_name_s{i}); 
    
>>>>>>> a7af7a1f535a2ec39beb6a80c571215473c615e9
    load(mat_path);

    % ------- Do some more processing so we get even color distribution----
    thresh = (20/7)*(1/1000)*120*2;
    I_fit = rdat.I_fit;
    [kmin, kmax] = find_raster_extents(rdat);
    I_temp = I_fit(:,kmin:kmax);
    I_fit = I_fit-mean(I_temp(:));
    I_fit = max(I_fit, -thresh);
    I_fit = min(I_fit, thresh);
    I_temp = I_temp-mean(I_temp(:));
    I_temp = I_fit(:,kmin:kmax);
    [min_s(i), ind_min] = min(min(I_temp));
     max_s(i) = max(max(I_fit));
<<<<<<< HEAD

    rdat.I_fit_normalized = I_fit;
    % ------ Pull out a couple cycles of the x-direction data -------

    xsamps = rdat.samps_per_period*3;

    x_data = csvread(data_path, 0,0, [0, 0, xsamps, 0]);
    rdat.x_data = x_data(end-rdat.samps_per_period*2+1:end);

    xy_ref = reshape(csvread(parent_path), 2, [])';
    rdat.x_ref = [xy_ref(:,1);xy_ref(:,1)];

    raster_dat_s = [raster_dat_s, rdat];
    figure(i+100)

=======
    
    rdat.I_fit_normalized = I_fit;
    % ------ Pull out a couple cycles of the x-direction data -------
    
    xsamps = rdat.samps_per_period*3;
    
    x_data = csvread(data_path, 0,0, [0, 0, xsamps, 0]);
    rdat.x_data = x_data(end-rdat.samps_per_period*2+1:end);
    
    xy_ref = reshape(csvread(parent_path), 2, [])';
    rdat.x_ref = [xy_ref(:,1);xy_ref(:,1)];
    
    raster_dat_s = [raster_dat_s, rdat];
    figure(i+100)
    
>>>>>>> a7af7a1f535a2ec39beb6a80c571215473c615e9
    imshow(raster_dat_s(i).I_fit_normalized, [min_s(i), max_s(i)])
end





%%
clc
F1 = figure(1); clf
% axes('FontName','Times New Roman') % Set axis font style
% box('on'); % Define box around whole figure


figwidth = 3.4;
figheight = 4.5;
set(F1, 'Units', 'Inches', 'Position', [3,3, figwidth, figheight],...
    'PaperUnits', 'Inches', 'PaperSize', [figwidth, figheight])
set(F1, 'Color', 'w');

xpad = .03;
ypad = 0.1;
wd = .82;

lft1 = .15; % start of first col;

ht1 = .37; % image heights
bt2 = .13;
bt1 = bt2+ypad+ht1;


ht2 = .37;

ax1 = axes('Position', [lft1, bt1, wd, ht1]);
% subplot(4,3,[2,5])
ax2 = axes('Position', [lft1, bt2, wd, ht2]);

F1.CurrentAxes = ax1;
%
hold on
offsets = ([1:.5:length(dat_name_s)]-1)*2;
for i=1:length(dat_name_s)
    xdat = raster_dat_s(i).x_data;
    x_ref = raster_dat_s(i).x_ref;
    t = [0:1:length(x_ref)-1]'*Ts;
    T = 1/raster_dat_s(1).meta.raster_freq;
%     T = 1
    plot(t/t(end), (x_ref+offsets(i))*volts2microns, 'k')
    plot(t/t(end), (xdat + offsets(i))*volts2microns, '--r')
end
xlabel('normalized time', 'interpreter', 'latex')
ylabel('x(t) [$\mu$m]',  'interpreter', 'latex')

F1.CurrentAxes = ax2;
hold on

kk = 100;
for i=1:length(dat_name_s)
    udat = raster_dat_s(i).I_fit_normalized(kk,:);
    xs = [0:1:raster_dat_s(i).npix-1]';
<<<<<<< HEAD

=======
    
>>>>>>> a7af7a1f535a2ec39beb6a80c571215473c615e9
    plot(xs, udat)

end
xlim([0,raster_dat_s(1).npix])
ylabel('$u_z$ [v]', 'interpreter', 'latex')
xlabel('$x$-direction pixel', 'interpreter', 'latex')

fig1_path = fullfile(getfigroot, '20micron_x_u_data.pdf');
export_fig(F1, fig1_path, '-q101')
%%
clc
F2 = figure(1); clf
% axes('FontName','Times New Roman') % Set axis font style
% box('on'); % Define box around whole figure


figwidth = 7.0;
figheight = 2.05;
set(F2, 'Units', 'Inches', 'Position', [-10,3, figwidth, figheight],...
    'PaperUnits', 'Inches', 'PaperSize', [figwidth, figheight])
set(F2, 'Color', 'w');


% subplot(4,3,[1,4])
xpad = .02;
wd = .2071;

lft1 = .035; % start of first col;

lft2 = lft1+wd+xpad;
lft3 = lft2+wd+xpad;
lft4 = lft3+wd+xpad;

ht_im = .7879; % image heights
bt_im = 0.09;


% subplot(1,3,1)
% ax1 = gca();
ax1 = axes('Position', [lft1, bt_im, wd, ht_im]);

% subplot(1,3,2)
% ax2 = gca();
ax2 = axes('Position', [lft2, bt_im, wd, ht_im]);

% subplot(1,3,3)
% ax3=gca();
ax3 = axes('Position', [lft3, bt_im, wd, ht_im]);

ax4 = axes('Position', [lft4, bt_im, wd, ht_im]);
% row 2

%
% ---------------- Start Plotting -----------------------------
% ------------------ Col 1: 1Hz -------------------------------
% pix_starts = [2,10, 22]
% pix_ends = [4, 16, 40]
pix_starts = [1,1, 1, 1];
pix_ends = [0, 0, 0, 0];

axr1 = [ax1, ax2, ax3, ax4]';
% minn = min(min_s);
% maxx = max(max_s);
minn = -thresh
maxx = thresh

% raster_dat_s = [raster_dat_s, raster_dat_s(end)];
for iter = 1:length(axr1)
    ax_iter = axr1(iter);
    F2.CurrentAxes = ax_iter;
<<<<<<< HEAD

=======
    
>>>>>>> a7af7a1f535a2ec39beb6a80c571215473c615e9
    I_plot = raster_dat_s(iter).I_fit_normalized;
    I_plot = I_plot(:, pix_starts(iter):end-pix_ends(iter));
    lo = min(min(I_plot));
    hi = max(max(I_plot));
%
    % ax1.XTick = [0 2 4]
    % ax1.YTick = [0 2 4]
    xdata = [0, raster_dat_s(iter).meta.width];
    ydata = [0, raster_dat_s(iter).meta.width];
<<<<<<< HEAD


=======
    
    
>>>>>>> a7af7a1f535a2ec39beb6a80c571215473c615e9
    imshow(I_plot, [minn, maxx],...
        'XData', xdata, 'YData', ydata, 'Parent', ax_iter)
    % imshow(I_fit, [lo, hi], 'Parent', ax1)


    axis('on')
    % ax1.YTickLabel = flipud({ax1.YTickLabel{2:end}})

%     xlabel('x-dir [$\mu$m]', 'interpreter', 'latex');
    if iter ==1
%         ylabel('y-dir [$\mu$m]',  'interpreter', 'latex');
    else
<<<<<<< HEAD
       set(ax_iter, 'YTickLabel', []);
=======
       set(ax_iter, 'YTickLabel', []); 
>>>>>>> a7af7a1f535a2ec39beb6a80c571215473c615e9
    end

    stit = sprintf('%.2f~Hz', raster_dat_s(iter).freq);
    title(stit, 'interpreter', 'latex');

end
lft5 = lft4-.02 +wd;
wd2 = wd/4;
%%
ax5 = axes('Position', [lft5, bt_im+0.02, wd2, ht_im - 0.05]);

axis('off')
h = colorbar;
colormap gray
% convert the threshold to nanometers
% (20/7)*(1/1000)*120*2;
thresh_color = 0.5*thresh*(7/20)*(1000/1);
caxis([-thresh_color, thresh_color])
% set(h, 'ylim', [-thresh_color, thresh_color])


fig_path = fullfile(getfigroot, '20micron_rasterscans.pdf');
export_fig(F2, fig_path, '-q101')
<<<<<<< HEAD
=======

>>>>>>> a7af7a1f535a2ec39beb6a80c571215473c615e9
