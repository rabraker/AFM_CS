clc
clear

% This path should be set to the root of your github working directory. 
if ~ispc
  addpath /home/arnold/gradschool/sysID/matlab/functions
end

addpath(fullfile(getMatPath(), 'afm-cs', 'matlab-code', 'functions'))
addpath(genpath(fullfile(getMatPath(), 'afm-cs', 'reconstruction', 'BP')))
addpath(fullfile(getMatPath(), 'afm-cs', 'reconstruction', 'SMP_1D'))
Ts = 40e-6;


% cs_exp_data_name_s{1} = 'cs-traj-z-bounce_out_10-17-2018-03.csv';
cs_exp_data_name_s{1} = 'cs-traj-512pix-10perc-500nm-5mic-01Hz_out_10-22-2018-06.csv';
sub_dir = '5microns/10-21-2018';
data_root = fullfile(PATHS.exp(), 'imaging', 'cs-imaging', sub_dir);
chan_map = ChannelMap([1:5]);

% cs_exp_data_name_s{1} = 'cs-traj-512pix-10perc-500nm-5mic-01Hz_out_9-23-2017-04.csv';
% data_root = '/media/labserver/acc2018-data/cs-data/5microns/9-22-2017';
% chan_map = ChannelMap([1:3, 5,6])

cs_paths = cs_exp_paths(data_root, cs_exp_data_name_s{1});

hole_depth = (20);

cs_exp = CsExp(cs_paths, chan_map, Ts, hole_depth);

fprintf('Total Imaging time: %.2f\n', cs_exp.time_total)
figbase = 100;


% Dz = zpk([0], [1], cs_exp.meta_exp.Ki, Ts);
% 
% models = load('G_zdir.mat');
% G2 = models.G2;
% LPF = G2*G2*G2*G2;
% 
% G1 = zpk(models.G1*models.G2);
% models = load(fullfile(PATHS.sysid, 'x-axis_sines_infoFourierCoef_10-21-2018-03.mat'));
% G1 = -models.modelFit.G_zdir;
% H2 = 1 + G1*Dz;
% H1 = minreal(H2/Dz);
% cs_exp.Gz; %= -H1;
% 
% cs_exp.Gz = zpk([], [], 1, Ts);
verbose = false;
if verbose
  figs{1}= figure(1000); clf; hold on, grid on;
  figs{2}= figure(2000); clf; hold on, grid on;
  figs{3}= figure(3000); clf; hold on, grid on;
else
  figs = [];
end
cs_exp.process_cs_data(verbose, figs);

%
pixelifsampled = cs_exp.pix_mask;
I = cs_exp.Img_raw;
fprintf('finished processing raw CS data...\n');
%%

bp = true;
% ********* SMP *************
clear CsExp
cs_exp.solve_smp1d();

ht = cs_exp.feature_height;
figure(10+figbase)
ax = gca();
figure(11+figbase)
axx = gca();
imshow_dataview(cs_exp.Img_raw, [-ht, ht], ax, axx)

figure(12+figbase)
ax = gca();
figure(13+figbase)
axx = gca();
imshow_dataview(cs_exp.Img_smp1d - mean(cs_exp.Img_smp1d(:)), [-ht, ht], ax, axx)


ht = cs_exp.feature_height;
f5 = figure(9 + 10*i); clf
subplot(2,3,1)
ax3 = gca();

subplot(2,3,2)
ax4 = gca();

subplot(2,3,3)
ax5 = gca();
    
subplot(2,3,[4,5,6])
ax6 = gca();
ax6.Visible = 'off';

cs_exp.Img_raw = cs_exp.Img_raw - mean(cs_exp.Img_raw(:));
cs_exp.Img_smp1d = cs_exp.Img_smp1d - mean(cs_exp.Img_smp1d(:));
imshow_sane(cs_exp.Img_raw, ax3, cs_exp.width, cs_exp.width, [-ht, ht]);
title(ax3, 'sample');

imshow_sane(cs_exp.Img_smp1d, ax5, cs_exp.width, cs_exp.width, [-ht, ht])
title(ax5, 'SMP reconstruction');
drawnow
%%
if bp
    cs_exp.solve_basis_pursuit();

    % imshow_sane(PixelVectorToMatrix(Ir_bp,[n m]), ax4, cs_exp.width, cs_exp.width);
    
    
    imshow_sane(cs_exp.Img_bp, ax4, cs_exp.width, cs_exp.width, [-ht, ht]);

    title(ax4, 'BP reconstruction');
end

figure(14+figbase)
ax = gca();
figure(15+figbase)
axx = gca();
imshow_dataview(cs_exp.Img_smp1d - mean(cs_exp.Img_smp1d(:)), [-ht, ht], ax, axx)

%%
axes(ax6)
s = metadata2text(cs_exp.meta_exp, Ts)
s = sprintf('%s\nperc=%.3f', s, sum(sum(pixelifsampled))/cs_exp.npix/cs_exp.npix);
s2 = sprintf('%s', cs_paths.data_path)
disp(s)
ax6.Visible = 'off';

t1 = text(0,.5, s, 'Units', 'normalized');
t2 = text(0, t1.Extent(2)-.1, s2, 'Units', 'normalized', 'interpreter', 'none');


figure(12)
ax = gca();
figure(13)
axx = gca();
imshow_dataview(cs_exp.Img_raw, [-ht, ht], ax, axx)

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

