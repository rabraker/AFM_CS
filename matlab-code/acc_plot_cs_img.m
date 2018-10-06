clc
clear
addpath('functions')
data_root = fullfile(getdataroot, 'cs-data');
cs_exp_data_name = 'cs-traj10-500_out_8-25-2017-14.csv';
cs_exp_data_name = 'cs-traj10-500_8-22-2017_08.csv';

% cs_exp_data_name = 'cs-traj10-500_8-22-2017_03.csv';

img_data_file_name = strrep(cs_exp_data_name, '.csv', '_img-data.mat');
img_data_path = fullfile(data_root, img_data_file_name);

% latex tables
table_name = strrep(cs_exp_data_name, '.csv', '_tabledata.tex');
table_path = fullfile(gettablesroot, table_name);

doc_fig_name = strrep(cs_exp_data_name, '.csv', '-fig.pdf');

load(img_data_path);

saveall = 1;

%%
% ---------------------- Plotting ----------------------------------
% --------------------------- Build Figure ------------------------------
width = 5; % This really needs to get saved into the meta data

close all
figwidth = 3.4;
figheight = 2;
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
bt1 = .13




ax1 = axes('Position', [lft1, bt1, wd, ht1]);
ax2 = axes('Position', [lft2, bt1, wd, ht1]);
set(ax2, 'YTickLabel', [])

height = width;
xdata = [0, width];
ydata = [0, height];

F1.CurrentAxes = ax1;
lo = min(min(img_data.cs_im));
hi = max(max(img_data.cs_im));
imshow(img_data.cs_im, [lo, hi], 'XData', xdata, 'YData', ydata, 'Parent', ax1);
title('sample');
axis('on')
xlabel('$x$-direction [$\mu$m]', 'interpreter', 'latex')
ylabel('$y$-direction [$\mu$m]', 'interpreter', 'latex')


F1.CurrentAxes = ax2;
lo = min(min(img_data.bp_im));
hi = max(max(img_data.bp_im));

imshow(img_data.bp_im,[lo, hi], 'XData', xdata, 'YData', ydata, 'Parent', ax2);
title('BP reconst.');
axis('on')
xlabel('$x$-direction [$\mu$m]', 'interpreter', 'latex')
set(ax2, 'YTickLabel', []) 
% ylabel('interaction force F(r) [N]', 'interpreter', 'latex')




if saveall
    s_latex = metadata2latex(img_data.meta, img_data.Ts);
    disp(s_latex)

    fid = fopen(table_path, 'w');

    fprintf(fid, '%s', s_latex);
    fclose(fid);

    export_fig(F1, fullfile(getfigroot(), doc_fig_name), '-q101')

end