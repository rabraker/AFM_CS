clc
clear
addpath('functions')
data_root = fullfile(getdataroot, 'cs-data');

cs_exp_data_name1 = 'cs-traj-8perc-500nm-5mic-1Hz_out_9-4-2017-08.csv';
cs_exp_data_name2 = 'cs-traj-8perc-500nm-5mic-1Hz_out_9-4-2017-09.csv';
% cs_exp_data_name = 'cs-traj10-500_8-22-2017_03.csv';

%----------------- Data set 1 -------------
img_data_file_name1 = strrep(cs_exp_data_name1, '.csv', '_img-data.mat');
img_data_path1 = fullfile(data_root, img_data_file_name1);

cs_data_path1 = fullfile(data_root, cs_exp_data_name1);
% latex tables
table_name1 = strrep(cs_exp_data_name1, '.csv', '_tabledata.tex');
table_path1 = fullfile(gettablesroot, table_name1);

doc_fig_name1 = strrep(cs_exp_data_name1, '.csv', '-compare_fig.pdf');
load(img_data_path1);
img_data1 = img_data;

% -------------- Data set 2 ----------------
img_data_file_name2 = strrep(cs_exp_data_name2, '.csv', '_img-data.mat');
img_data_path2 = fullfile(data_root, img_data_file_name2);
cs_data_path2 = fullfile(data_root, cs_exp_data_name2);
% latex tables
table_name2 = strrep(cs_exp_data_name2, '.csv', '_tabledata.tex');
table_path2 = fullfile(gettablesroot, table_name2);
doc_fig_name2 = strrep(cs_exp_data_name2, '.csv', '-fig.pdf');

load(img_data_path2);
img_data2 = img_data;


saveall = 0;

% -------------- Get line data1-----------------------
dat1 = csvread(cs_data_path1); 
% Drop all data corresponding to movement between points.
ind_meas = find(dat1(:, 6) > 0);  % Index of the movements are all 0.
% meta_ind = dat1(:,6);
dat_meas = dat1(ind_meas, :);
dat_meas(:,5) = detrend(dat_meas(:,5)); % detrend u_z as a long vector

line_data = {};
for k = 1:max(dat_meas(:,end))
    inds = find(dat_meas(:,end) == k); 
    line_data{k} = dat_meas(inds, 5);        
end

% -------------- Get line data2-----------------------
dat2 = csvread(cs_data_path2); 
% Drop all data corresponding to movement between points.
ind_meas2 = find(dat2(:, 6) > 0);  % Index of the movements are all 0.
% meta_ind = dat1(:,6);
dat_meas2 = dat2(ind_meas2, :);
dat_meas2(:,5) = detrend(dat_meas2(:,5)); % detrend u_z as a long vector

line_data2 = {};
for k = 1:max(dat_meas2(:,end))
    inds = find(dat_meas2(:,end) == k); 
    line_data2{k} = dat_meas2(inds, 5);        
end

%%
% ---------------------- Plotting ----------------------------------
% --------------------------- Build Figure ------------------------------
width = 5; % This really needs to get saved into the meta data

close all
figwidth = 3.4;
figheight = 3.4;
F1 = figure(1); clf;
set(F1, 'Units', 'Inches', 'Position', [12,0, figwidth, figheight],...
    'PaperUnits', 'Inches', 'PaperSize', [figwidth, figheight])
set(F1, 'Color', 'w');

xpad = .04;
lft1 = .12;
wd = .4;
lft2 = lft1 + xpad + wd;
ht1 = .85;
bt2 = .2;
% bt1 = bt2 + ypad + ht1;
bt1 = .31;

bt0 = .1;
ht2 = .25;


ax1 = axes('Position', [lft1, bt1, wd, ht1]);
ax2 = axes('Position', [lft2, bt1, wd, ht1]);
ax3 = axes('Position', [lft1, bt0, wd, ht2]);
ax4 = axes('Position', [lft2, bt0, wd, ht2]);

set(ax2, 'YTickLabel', [])

height = width;
xdata = [0, width];
ydata = [0, height];

F1.CurrentAxes = ax1;
lo = min(min(img_data1.bp_im));
hi = max(max(img_data1.bp_im));
imshow(img_data1.bp_im, [lo, hi], 'XData', xdata, 'YData', ydata, 'Parent', ax1);
stit1 = sprintf('$\\delta_z = %f$', img_data1.meta.ramp_rate);
title(stit1, 'interpreter', 'latex');
axis('on')
xlabel('$x$-direction [$\mu$m]', 'interpreter', 'latex')
ylabel('$y$-direction [$\mu$m]', 'interpreter', 'latex')

clc

F1.CurrentAxes = ax2;
lo = min(min(img_data2.bp_im));
hi = max(max(img_data2.bp_im));

imshow(img_data2.bp_im,[lo, hi], 'XData', xdata, 'YData', ydata, 'Parent', ax2);
stit2 = sprintf('$\\delta_z = %f$', img_data2.meta.ramp_rate);
title(stit2, 'interpreter', 'latex');

axis('on')
xlabel('$x$-direction [$\mu$m]', 'interpreter', 'latex')
set(ax2, 'YTickLabel', []) 
% ylabel('interaction force F(r) [N]', 'interpreter', 'latex')

F1.CurrentAxes = ax3;
k = 10;
t = [0:1:length(line_data{k})-1]'*img_data1.Ts;
plot(t, line_data{k}-max(line_data{k}))
xlim([t(1), t(end)])
xl1 = xlabel('time [s]', 'interpreter', 'latex');
set(xl1, 'Position', [0.0253 -0.11 -1.0000])
F1.CurrentAxes = ax4;
k = 12;
t = [0:1:length(line_data2{k})-1]'*img_data1.Ts;
plot(t, line_data2{k}-max(line_data2{k}))
xlim([t(1), t(end)])
linkaxes([ax3, ax4], 'y')
set(ax4, 'YTickLabel', [])
xl2 = xlabel('time [s]', 'interpreter', 'latex');
set(xl2, 'Position', [0.0253 -0.11 -1.0000])

if saveall
    s_latex1 = metadata2latex(img_data1.meta, img_data1.Ts);
    s_latex2 = metadata2latex(img_data2.meta, img_data2.Ts);
    disp(s_latex1)

    fid = fopen(table_path1, 'w');
    fprintf(fid, '%s', s_latex1);
    fclose(fid);
    fid = fopen(table_path2, 'w');
    fprintf(fid, '%s', s_latex2);
    fclose(fid);

    export_fig(F1, fullfile(getfigroot(), doc_fig_name1), '-q101')

end