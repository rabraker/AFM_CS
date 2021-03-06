% ------------------- DRAW FIGURES FOR 20 MICRON CS IMAGES ----------------

addpath('functions')

data_root = '/media/labserver/acc2018-data';
% dat_root = '/media/labserver/acc2018-data/data/';


% -------- Constants --------
Ts = 40e-6;
cs_root =  fullfile(data_root,'cs-data');
cs_name_s{1} = 'cs-traj-512pix-5perc-500nm-5mic-01Hz_out_9-23-2017-04.csv';
cs_name_s{2} = 'cs-traj-512pix-10perc-500nm-5mic-01Hz_out_9-23-2017-04.csv';
cs_name_s{3} = 'cs-traj-512pix-15perc-500nm-5mic-01Hz_out_9-23-2017-03.csv';

sub_fold = '5microns/9-22-2017';
width = 20;

xdata = [0, width];
ydata = [0, width];

max_cs = zeros(1,length(cs_name_s));
min_cs = zeros(1,length(cs_name_s));
for i=1:length(cs_name_s)
    img_data_name = strrep(cs_name_s{i}, '.csv', '_img-data.mat');
    parent_name = get_parent_name(cs_name_s{i}, '_out_');
    meta_name = strrep(parent_name, '.csv', '.mat');
    load(fullfile(cs_root,sub_fold, meta_name));
    
    load(fullfile(cs_root, sub_fold, img_data_name))
    img_data.CsExpMetaIn = CsExpMetaIn;
    pixelifsampled = img_data.CsExpMetaIn.pixelifsampled;
    pix = img_data.CsExpMetaIn.npix;
    actual_sub_sample_frac = length(find(pixelifsampled == 1))/pix^2;
    
    cs_img_data_s{i} = img_data;
    cs_img_data_s{i}.img_data.perc = actual_sub_sample_frac*100;
    bp_im = img_data.bp_im;
    
    bp_im = bp_im - mean(bp_im(:));
%     bp_im = bp_im - min(bp_im(:));
    min_cs(i) = min(min(bp_im));
    max_cs(i) = max(max(bp_im));
    
    cs_img_data_s{i}.bp_im = bp_im;
end

min_s = zeros(length(cs_name_s), 1);
max_s = zeros(length(cs_name_s), 1);

for i=1:length(cs_name_s)
    cs_img_data_s{i}.bp_im_normalized = cs_img_data_s{i}.bp_im - mean(cs_img_data_s{i}.bp_im(:));
    figure(i+length(cs_name_s)+100);
    imshow(cs_img_data_s{i}.bp_im_normalized, [min_cs(i), max_cs(i)])
%     [0, DRng])
%     ,  'XData', xdata, 'YData', ydata)
    axis('on') 
    drawnow
    title(sprintf('bp: %d%%', floor(cs_img_data_s{i}.img_data.perc)))
    ax = gca;

    grid on
      set(ax,'xminorgrid','on','yminorgrid','on')
      colorbar
      min_s(i) = min(cs_img_data_s{i}.bp_im_normalized(:));
      max_s(i) = max(cs_img_data_s{i}.bp_im_normalized(:));
   
end
save(fullfile(cs_root, sub_fold, 'cs_img_data_s.mat'), 'cs_img_data_s')
minn = min(min_s);
maxx = max(max_s);


%%

F3 = figure(3+100); clf
figwidth = 7.5;
figheight = 2.1;
set(F3, 'Units', 'Inches', 'Position', [-10,3, figwidth, figheight],...
    'PaperUnits', 'Inches', 'PaperSize', [figwidth, figheight])
set(F3, 'Color', 'w');


% subplot(4,3,[1,4])
xpad = -.04;
wd = .27071;
wd2 = wd/4;
lft1 = .15; % start of first col;

lft2 = lft1+wd+xpad;
lft3 = lft2+wd+xpad;
lft4 = lft3+wd+-0.04;

ht_im = .75; % image heights
bt_im = 0.13;

ax0 = axes('Position', [0, bt_im, lft1, ht_im]);

annotation('textarrow', [0.05, 0.05], [.5, .75]);
annotation('textarrow', [0.05, 0.125], [.5, .5]);
text(0.5, 0.4, '$x$', 'interpreter', 'latex', 'FontSize', 14)
text(0.15, 0.6, '$y$', 'interpreter', 'latex', 'FontSize', 14)
set(ax0, 'Visible', 'off')

ax1 = axes('Position', [lft1, bt_im, wd, ht_im]);
ax2 = axes('Position', [lft2, bt_im, wd, ht_im]);
ax3 = axes('Position', [lft3, bt_im, wd, ht_im]);

pix_starts = [1,1, 1, 1];
pix_ends = [0, 0, 0, 0];

axr1 = [ax3, ax2, ax1]';

% plot a sampling pattern
F3.CurrentAxes = ax1;
xdata = [0, cs_img_data_s{1}.width];
ydata = [0, cs_img_data_s{1}.width];


% axis('on')
thresh = (20/7)*(1/1000)*22;
%
% cs_img_data_s{3} = cs_img_data_s{2};
for iter = 1:length(axr1)
    ax_iter = axr1(iter);
    F3.CurrentAxes = ax_iter;
    
    I_plot = cs_img_data_s{iter}.bp_im_normalized;
    I_plot = I_plot(:, pix_starts(iter):end-pix_ends(iter));
    lo = min(min(I_plot));
    hi = max(max(I_plot));

    xdata = [0, cs_img_data_s{i}.width];
    ydata = [0, cs_img_data_s{i}.width];
    
    imshow(I_plot, [-thresh, thresh],...
        'XData', xdata, 'YData', ydata, 'Parent', ax_iter)

    axis('on')
   if iter ~=3
       set(ax_iter, 'YTickLabel', []); 
   end
    stit = sprintf('$\\mu$-path: %.2f\\%%', cs_img_data_s{iter}.img_data.perc);
    title(stit, 'interpreter', 'latex');

end
F3.CurrentAxes = ax1;
set(ax1, 'YTick', [0, 2, 4])

ax4 = axes('Position', [lft4, bt_im+0.02, wd2, ht_im - 0.05]);

axis('off')
h = colorbar;
colormap gray
% convert the threshold to nanometers
% (20/7)*(1/1000)*120*2;
thresh_color = thresh*(7/20)*(1000/1); %should give 22
caxis([-thresh_color, thresh_color])
text(10.5, .3, '$z$ [nm]', 'rot', 90, 'interpreter', 'latex',...
    'Units', 'normalized', 'FontSize', 14)

fig_path = fullfile(getfigroot, '5micron_csscans_v2.eps');
% export_fig(F1, fig_path, '-q101')
saveEps(F3, fig_path)



%%
% % F3 = figure(3+100); clf
% % figwidth = 6;
% % figheight = 2.1;
% % set(F3, 'Units', 'Inches', 'Position', [-10,3, figwidth, figheight],...
% %     'PaperUnits', 'Inches', 'PaperSize', [figwidth, figheight])
% % set(F3, 'Color', 'w');
% % 
% % 
% % % subplot(4,3,[1,4])
% % xpad = .03;
% % wd = .27071;
% % wd2 = wd/4;
% % lft1 = .035; % start of first col;
% % 
% % lft2 = lft1+wd+xpad;
% % lft3 = lft2+wd+xpad;
% % lft4 = lft3+wd+-0.02;
% % 
% % ht_im = .7879; % image heights
% % bt_im = 0.13;
% % 
% % ax1 = axes('Position', [lft1, bt_im, wd, ht_im]);
% % ax2 = axes('Position', [lft2, bt_im, wd, ht_im]);
% % ax3 = axes('Position', [lft3, bt_im, wd, ht_im]);
% % 
% % pix_starts = [1,1, 1, 1];
% % pix_ends = [0, 0, 0, 0];
% % 
% % axr1 = [ax3, ax2, ax1]';
% % 
% % % plot a sampling pattern
% % F3.CurrentAxes = ax1;
% % xdata = [0, cs_img_data_s{1}.width];
% % ydata = [0, cs_img_data_s{1}.width];
% % 
% % 
% % % axis('on')
% % thresh = (20/7)*(1/1000)*22;
% % %
% % % cs_img_data_s{3} = cs_img_data_s{2};
% % for iter = 1:length(axr1)
% %     ax_iter = axr1(iter);
% %     F3.CurrentAxes = ax_iter;
% %     
% %     I_plot = cs_img_data_s{iter}.bp_im_normalized;
% %     I_plot = I_plot(:, pix_starts(iter):end-pix_ends(iter));
% %     lo = min(min(I_plot));
% %     hi = max(max(I_plot));
% % 
% %     xdata = [0, cs_img_data_s{i}.width];
% %     ydata = [0, cs_img_data_s{i}.width];
% %     
% %     imshow(I_plot, [-thresh, thresh],...
% %         'XData', xdata, 'YData', ydata, 'Parent', ax_iter)
% % 
% %     axis('on')
% %    if iter ~=3
% %        set(ax_iter, 'YTickLabel', []); 
% %    end
% %     stit = sprintf('$\\mu$-path: %.2f\\%%', cs_img_data_s{iter}.img_data.perc);
% %     title(stit, 'interpreter', 'latex');
% % 
% % end
% % F3.CurrentAxes = ax1;
% % set(ax1, 'YTick', [0, 2, 4])
% % 
% % ax4 = axes('Position', [lft4, bt_im+0.02, wd2, ht_im - 0.05]);
% % 
% % axis('off')
% % h = colorbar;
% % colormap gray
% % % convert the threshold to nanometers
% % % (20/7)*(1/1000)*120*2;
% % thresh_color = thresh*(7/20)*(1000/1); %should give 22
% % caxis([-thresh_color, thresh_color])

%%
fig_path = fullfile(getfigroot, '5micron_csscans_v2.eps');
% export_fig(F1, fig_path, '-q101')
saveEps(F3, fig_path)
%%
clc
for i =1:length(cs_img_data_s)
    
time_s = cs_img_data_s{i}.meta.state_counts*Ts;
time_s(2) = time_s(2) + time_s(3);
time_s(3) = [];
time_tot = sum(time_s);
s = '';
for k = 1:length(time_s)
   s = sprintf('%s%.2f & ', s,(time_s(k)/time_tot)*100); 
end
s = sprintf('%s&%.2f\\\\\n', s, time_tot);
fprintf('State times as percentage (%.2f sampling)\n%s',cs_img_data_s{i}.img_data.perc, s)
end


