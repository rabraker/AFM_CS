clc
clear

data_root = fullfile(getdataroot, 'data', 'cs-data');
cs_exp_data_name = 'cs-traj10-500_out_8-25-2017-14.csv';
cs_exp_data_name = 'cs-traj10-500_8-22-2017_08.csv';

img_data_file_name = strrep(cs_exp_data_name, '.csv', '_img-data.mat');
img_data_path = fullfile(data_root, img_data_file_name);

load(img_data_path)
%%

f5 = figure(6); clf
width = img_data.width;
subplot(2,3,1)
ax3 = gca();
imshow_sane(img_data.cs_im, ax3, width, width);
title('sample');

subplot(2,3,2)
ax4 = gca();
% imshow_sane(PixelVectorToMatrix(Ir_bp,[n m]), ax4, width, width);
imshow_sane(img_data.bp_im, ax4, width, width);
title('BP reconstruction');

subplot(2,3,3)
ax5 = gca();
imshow_sane(img_data.smp_im, ax5, width, width)
title('SMP reconstruction');

s = metadata2text(img_data.meta, img_data.Ts);
disp(s)
