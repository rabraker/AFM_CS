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
dat_name_s{1} = 'raster_scan_512pix_5mic_2.50e-01Hz_out_9-23-2017-02.csv';
dat_name_s{2} = 'raster_scan_512pix_5mic_01Hz_out_9-23-2017-04.csv';
dat_name_s{3} = 'raster_scan_512pix_5mic_05Hz_out_9-23-2017-03.csv';
dat_name_s{4} = 'raster_scan_512pix_5mic_10Hz_out_9-23-2017-01.csv';
sub_fold = '5microns/9-22-2017';



raster_dat_s = [];
max_s = zeros(1,4);
min_s = zeros(1,4);
for i=1:length(dat_name_s)
    data_path= fullfile(data_root,'raster', sub_fold, dat_name_s{i});
    parent_name = get_parent_name(dat_name_s{i}, '_out_');
    parent_path = fullfile(data_root, 'raster', sub_fold, parent_name);
    parent_meta_path = strrep(parent_path, '.csv', '.mat');
    meta_path = strrep(data_path, '.csv', '-meta.mat'); % Labview output
    mat_path = strrep(data_path, '.csv', '.mat'); % processed matlab output
    
    fprintf('File: %s\n', dat_name_s{i}); 
    
    load(mat_path);

    % ------- Do some more processing so we get even color distribution----
    thresh = (20/7)*(1/1000)*22;
    I_fit = rdat.I_fit;
    
%     figure; imagesc(I_fit); colormap gray
    
    [kmin, kmax] = find_raster_extents(rdat);
    I_temp = I_fit(:,kmin:kmax);
    I_mu = mean(I_fit(:));
    I_fit = I_fit - I_mu;
    I_fit = max(I_fit, -thresh);
    I_fit = min(I_fit, thresh);
    I_fit = I_fit - mean(I_fit(:));
%     I_temp = I_temp-mean(I_temp(:));
%     I_temp = I_fit(:,kmin:kmax);
%     [min_s(i), ind_min] = min(min(I_temp));
    [min_s(i), ind_min] = min(min(I_fit));
     max_s(i) = max(max(I_fit));
    
     
    rdat.I_fit_normalized = I_fit;
    % ------ Pull out a couple cycles of the x-direction data -------
    
    xsamps = rdat.samps_per_period*3;
    
    x_data = csvread(data_path, 0,0, [0, 0, xsamps, 0]);
    rdat.x_data = x_data(end-rdat.samps_per_period*2+1:end);
    
    xy_ref = reshape(csvread(parent_path), 2, [])';
    rdat.x_ref = [xy_ref(:,1);xy_ref(:,1)];
    
    raster_dat_s = [raster_dat_s, rdat];
    figure(i+100)
    
    imshow(raster_dat_s(i).I_fit_normalized, [min_s(i), max_s(i)])
end




%%
close all
clc
F1 = figure(1); clf
% axes('FontName','Times New Roman') % Set axis font style
% box('on'); % Define box around whole figure


figwidth = 3.4;
figheight = 4.5;
set(F1, 'Units', 'Inches', 'Position', [-13,3, figwidth, figheight],...
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
ylabel('x(t)',  'interpreter', 'latex')

F1.CurrentAxes = ax2;
hold on

kk = 23;
hands = []
for i=1:length(dat_name_s)
    udat = raster_dat_s(i).I_fit_normalized(kk,:);
    xs = [0:1:raster_dat_s(i).npix-1]';
    
    hi = plot(xs, udat);
    hi.DisplayName = sprintf('%.2fHz', raster_dat_s(i).freq);
    hands = [hands;hi];
end
%     plot([xs(1), xs(end)], [-thresh, -thresh], '--k')
%     plot([xs(1), xs(end)], [thresh, thresh], '--k')

ylim([-0.06, 0.08])
xlim([0,raster_dat_s(1).npix])
hylab = ylabel('$u_z$ [v]', 'interpreter', 'latex');
set(hylab, 'Units', 'normalized', 'Position', [-0.1223 0.5000 0])
xlabel('$x$-direction pixel', 'interpreter', 'latex')


leg1 = legend(ax2, hands(1:2));
leg1.Box = 'off';
leg1.Position = [0.3978 0.3893 0.2791 0.0799];
ax3 = axes('position', get(ax2, 'position'), 'visible', 'off');
leg2 = legend(ax3, hands(3:4));
leg2.Position = [0.6708 0.3893 0.2577 0.0799];

leg2.Box = 'off';


%%
fig1_path = fullfile(getfigroot, '5micron_x_u_data.pdf');
export_fig(F1, fig1_path, '-q101')
% saveEps(F1, fig1_path)
%%
clc
F2 = figure(1); clf
% axes('FontName','Times New Roman') % Set axis font style
% box('on'); % Define box around whole figure


figwidth = 7.0;
figheight = 1.9;
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
bt_im = 0.105;


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
% minn = -thresh
% maxx = thresh

% raster_dat_s = [raster_dat_s, raster_dat_s(end)];
for iter = 1:length(axr1)
    ax_iter = axr1(iter);
    F2.CurrentAxes = ax_iter;
    
    I_plot = raster_dat_s(iter).I_fit_normalized;
    mean(I_plot(:))
    I_plot = I_plot(:, pix_starts(iter):end-pix_ends(iter));
    lo = min(min(I_plot));
    hi = max(max(I_plot));
%
    % ax1.XTick = [0 2 4]
    % ax1.YTick = [0 2 4]
    xdata = [0, raster_dat_s(iter).meta.width];
    ydata = [0, raster_dat_s(iter).meta.width];
    
    
    imshow(I_plot, [-thresh, thresh],...
        'XData', xdata, 'YData', ydata, 'Parent', ax_iter)
    % imshow(I_fit, [lo, hi], 'Parent', ax1)


    axis('on')
%     xlabel('x-dir [$\mu$m]', 'interpreter', 'latex');
    if iter ==1
%         ylabel('y-dir [$\mu$m]',  'interpreter', 'latex');
    else
       set(ax_iter, 'YTickLabel', []); 
    end

    stit = sprintf('%.2f~Hz', raster_dat_s(iter).freq);
    title(stit, 'interpreter', 'latex');

end
lft5 = lft4-.02 +wd;
wd2 = wd/4;

ax5 = axes('Position', [lft5, bt_im+0.02, wd2, ht_im - 0.05]);

axis('off')
h = colorbar;
colormap gray
% convert the threshold to nanometers
% (20/7)*(1/1000)*120*2;
thresh_color = thresh*(7/20)*(1000/1); % should give 22
caxis([-thresh_color, thresh_color])
% set(h, 'ylim', [-thresh_color, thresh_color])


fig_path = fullfile(getfigroot, '5micron_rasterscans_v2.eps');
% export_fig(F2, fig_path, '-q101')
saveEps(F2, fig_path)
%%


% -----------------------------------------------------------------------%
% Calculate metrics of ground truth against the faster raster scans.
clc
REF_ind = 1;
for iter=1:length(raster_dat_s)
    ssim_s = [];
    ind_s = {};
    j = 1;
    for k=1:25
        for L=1:25

        REF = raster_dat_s(ref_ind).I_fit_normalized(1:end-L+1, k:end);
        Y = raster_dat_s(iter).I_fit_normalized(L:end, 1:end-k+1);
        
            DRng = max(REF(:)) - min(REF(:));

            ssim_s = [ssim_s; ssim(Y, REF, 'DynamicRange', DRng)];

            ind_s{j} = [L, k];
            j = j+1;
        end
    end


    [ssim_max, ind_max] = max(ssim_s);
    L = ind_s{ind_max}(1);
    k = ind_s{ind_max}(2);
%     REF = raster_dat_s(REF_ind).I_fit_normalized(L:end, 1:end-k+1);
%     Y = raster_dat_s(iter).I_fit_normalized(1:end-L+1, k:end);
    REF = raster_dat_s(ref_ind).I_fit_normalized(1:end-L+1, k:end);
    Y = raster_dat_s(iter).I_fit_normalized(L:end, 1:end-k+1);
    psnr_max = psnr(Y, REF, DRng);
    freq = raster_dat_s(iter).freq;
    fprintf('Freq: %.2f,  ssim: %.3f,  psnr: %.2f\n', freq, ssim_max, psnr_max)
    fprintf('opt L, k: %d, %d\n', L,k)

end


% % figure(5);clf
% % plot(REF(41, :))
% % hold on
% % plot(Y(41,:))

%%
% ------------------------------------------------------------------------%
% Calculate Metrics of ground truth raster image against the CS images. 

cs_root =  fullfile(data_root,'cs-data');
load(fullfile(cs_root, sub_fold, 'cs_img_data_s.mat'))

figure(200)
imshow(cs_img_data_s{1}.bp_im_normalized,[minn, maxx] )
%%
clc
ref_ind = 1;

for i=1:3
   I = cs_img_data_s{i}.bp_im_normalized; 
   for krow = 1:size(I,1)
       the_row = detrend(I(krow,:)')'; 
%        keyboard
      I(krow, :) = the_row;
   end
   I = I - mean(I(:));
    cs_img_data_s{i}.I_metric = I;
    
end

%%

for ind = 1:3
ssim_s = [];
ind_s = {};
j = 1;
for k=1:30
    for L=1:25

        REF = raster_dat_s(ref_ind).I_fit_normalized(1:end-L+1, k:end);
%         Y = cs_img_data_s{ind}.bp_im_normalized(L:end, 1:end-k+1);
        Y = cs_img_data_s{ind}.I_metric(L:end, 1:end-k+1);

        DRng = max(REF(:)) - min(REF(:));


        ssim_s = [ssim_s; ssim(Y, REF, 'DynamicRange', DRng)];

        ind_s{j} = [L, k];
        j = j+1;
    end
end

[ssim_max, ind_max] = max(ssim_s);

L = ind_s{ind_max}(1);
k = ind_s{ind_max}(2);
REF = raster_dat_s(ref_ind).I_fit_normalized(1:end-L+1, k:end);
% Y = cs_img_data_s{ind}.bp_im_normalized(L:end, 1:end-k+1);
Y = cs_img_data_s{ind}.I_metric(L:end, 1:end-k+1);

% REF = raster_dat_s(ref_ind).I_fit_normalized(L:end, 1:end-k+1);
% Y = cs_img_data_s{ind}.bp_im_normalized(1:end-L+1, k:end);

psnr_max = psnr(Y, REF, DRng);

perc = cs_img_data_s{ind}.img_data.perc;


fprintf('perc: %.2f,  ssim: %.3f,  psnr: %.2f\n', perc, ssim_max, psnr_max)
end
%%
figure(5);clf
plot(REF(25, :))
hold on
plot(Y(25,:))