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

cs_exp.process_cs_data(true, figs);

%%
pixelifsampled = cs_exp.pix_mask;
I = cs_exp.Img;

figure(100 + 10*i);
ax = gca();
imshow(I, [min(min(I)), max(max(I))])


bp = true;

% ********* SMP *************

[n,m] = size(I);

% Ir_smp = SMP(I.*pixelifsampled,pixelifsampled,round(0.02*n*m));

tic
samp_frac = sum(sum(pixelifsampled))/(n*m)
reduce_frac = 0.1;  %It's my understanding Yufan says 1/10 of sub-sample fraction.
maxiter = round(samp_frac*reduce_frac*n*m)
Ir_smp = SMP_1D(I, pixelifsampled, maxiter);
time_smp = toc;

f5 = figure(9 + 10*i); clf
subplot(2,3,1)
ax3 = gca();
imshow_sane(I, ax3, cs_exp.width, cs_exp.width);
title('sample');

subplot(2,3,3)
ax5 = gca();
imshow_sane(Ir_smp, ax5, cs_exp.width, cs_exp.width)
title('SMP reconstruction');
drawnow

if bp
    % ********* BP *************
    [n m] = size(I);
    tic
    I_vector = PixelMatrixToVector(I);

    pixelifsampled_vec = PixelMatrixToVector(pixelifsampled);
    I_vector = I_vector(find(pixelifsampled_vec>0.5));

    A = @(x) IDCTfun(x,pixelifsampled_vec);
    At = @(x) DCTfun(x,pixelifsampled_vec);

    Ir_bp = idct(l1qc_logbarrier(At(I_vector), A, At, I_vector, 0.1));
    Ir_bp = real(Ir_bp);
    bp_im = PixelVectorToMatrix(Ir_bp,[n m]);
    time_bp = toc;
    fprintf('SMP Time: %f \nBP Time: %f\n', time_smp, time_bp);

    bp_im = detrend_plane(bp_im);
    
    subplot(2,3,2)
    ax4 = gca();
    % imshow_sane(PixelVectorToMatrix(Ir_bp,[n m]), ax4, cs_exp.width, cs_exp.width);
    
    imshow_sane(bp_im, ax4, cs_exp.width, cs_exp.width);

    title('BP reconstruction');
end
%
clc
s = metadata2text(ExpMetaData, Ts)
s = sprintf('%s\nperc=%.3f', s, sum(sum(pixelifsampled))/cs_exp.npix/cs_exp.npix);
s2 = sprintf('%s', cs_exp_data_name)
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

