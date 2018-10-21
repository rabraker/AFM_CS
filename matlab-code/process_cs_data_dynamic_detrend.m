clc
clear

% This path should be set to the root of your github working directory. 


addpath(fullfile(getMatPath(), 'afm-cs', 'matlab-code', 'functions'))
addpath(genpath(fullfile(getMatPath(), 'afm-cs', 'reconstruction', 'BP')))
addpath(fullfile(getMatPath(), 'afm-cs', 'reconstruction', 'SMP_1D'))
Ts = 40e-6;


cs_exp_data_name_s{1} = 'cs-traj-z-bounce_out_10-17-2018-03.csv';
cs_exp_data_name_s{1} = 'cs-traj-512pix-10perc-500nm-5mic-01Hz_out_10-15-2018-01.csv';

sub_dir = '5microns/10-14-2018';
data_root = fullfile(PATHS.exp(), 'imaging', 'cs-imaging', sub_dir);

gk = load('g_k.mat', 'g_k');
gk = gk.g_k


cs_paths = cs_exp_paths(data_root, cs_exp_data_name_s{1});

clear CsExp
hole_depth = (20);
chan_map = ChannelMap(1:5);

cs_exp = CsExp(cs_paths, chan_map, Ts, hole_depth);
cs_exp.idx_state_s

fprintf('Total Imaging time: %.2f\n', cs_exp.time_total)

%%
% f2 = figure(2+10*i); clf
% ylabel('u_z')
% xlabel('pixel')
% ax2 = gca();
% hold on;

% Since we want to convert xy-coords to pixels, move all data to the
% positive orthant:

clc
clear CsExp;
cs_exp.Gz = gk;
% cs_exp.Gz = zpk([], [], 1, Ts);
figs{1}= figure(1000); clf; hold on, grid on;
figs{2}= figure(2000); clf; hold on, grid on;
figs{3}= figure(3000); clf; hold on, grid on;

cs_exp.process_cs_data(false, figs);

%%
pixelifsampled = cs_exp.pix_mask;
I = cs_exp.Img_raw;

figure(100 + 10*i);
ax = gca();
imshow(cs_exp.Img_raw, [min(cs_exp.Img_raw(:)), max(cs_exp.Img_raw(:))])


bp = true;

% ********* SMP *************
%%
clear CsExp
cs_exp.solve_smp1d();
%%


f5 = figure(9 + 10*i); clf
subplot(2,3,1)
ax3 = gca();
imshow_sane(cs_exp.Img_raw, ax3, cs_exp.width, cs_exp.width);
title('sample');

subplot(2,3,3)
ax5 = gca();
imshow_sane(cs_exp.Img_smp1d, ax5, cs_exp.width, cs_exp.width)
title('SMP reconstruction');
drawnow

if bp
    cs_exp.solve_basis_pursuit();
    subplot(2,3,2)
    ax4 = gca();
    % imshow_sane(PixelVectorToMatrix(Ir_bp,[n m]), ax4, cs_exp.width, cs_exp.width);
    
    imshow_sane(cs_exp.Img_bp, ax4, cs_exp.width, cs_exp.width);

    title('BP reconstruction');
end

clc
s = metadata2text(cs_exp.meta_exp, Ts)
s = sprintf('%s\nperc=%.3f', s, sum(sum(pixelifsampled))/cs_exp.npix/cs_exp.npix);
s2 = sprintf('%s', cs_paths.data_path)
disp(s)
subplot(2,3,[4,5,6])
ax4 = gca();
ax4.Visible = 'off';

t1 = text(0,.5, s, 'Units', 'normalized');
t2 = text(0, t1.Extent(2)-.1, s2, 'Units', 'normalized', 'interpreter', 'none');

% text(0,-1.2, s, 'Units', 'normalized')
%%
savedata = 1;
if savedata
   img_data.cs_im = I;
   img_data.bp_im = bp_im;

   img_data.smp_im = Ir_smp;
   img_data.pixelifsampled = pixelifsampled;
   
   img_data.width = width;
   img_data.meta = ExpMetaData;
   img_data.Ts = Ts;
   img_data_file_name = strrep(cs_exp_data_name, '.csv', '_img-data.mat');
   img_data_path = fullfile(data_root, img_data_file_name);
   
   save(img_data_path, 'img_data')
    
end
%
% fig_root = 'C:\Users\arnold\Documents\labview\afm_imaging\matlab-code\figures';
cs_exp_fig_name = strrep(cs_data_path, '.csv', '-fig.fig')

fig_path = fullfile(cs_exp_fig_name);
saveas(f5, fig_path)

