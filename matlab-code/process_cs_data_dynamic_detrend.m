clc
clear

% This path should be set to the root of your github working directory. 


addpath(fullfile(getMatPath(), 'afm-cs', 'matlab-code', 'functions'))
addpath(genpath(fullfile(getMatPath(), 'afm-cs', 'reconstruction', 'BP')))
addpath(fullfile(getMatPath(), 'afm-cs', 'reconstruction', 'SMP_1D'))
Ts = 40e-6;


cs_exp_data_name_s{1} = 'cs-traj-z-bounce_out_10-17-2018-03.csv';
cs_exp_data_name_s{1} = 'cs-traj-512pix-10perc-500nm-5mic-01Hz_out_10-15-2018-02.csv';

channel_map = containers.Map({'x', 'y', 'ze', 'uz', 'meta'}, {1,2,3,4,5})
x_idx = 1;
y_idx = 2;
ze_idx = 3;
uz_idx = 4;
met_idx = 5;

sub_dir = '5microns/10-14-2018';
gk = load('g_k.mat', 'g_k');
gk = gk.g_k


cs_exp_data_name = cs_exp_data_name_s{1};
data_root = fullfile(PATHS.exp(), 'imaging', 'cs-imaging', sub_dir);

cs_exp_meta_name = strrep(cs_exp_data_name, '.csv', '-meta.mat');

cs_data_path = fullfile(data_root, cs_exp_data_name);
cs_meta_path = fullfile(data_root, cs_exp_meta_name);

dat = csvread(cs_data_path); 
load(cs_meta_path);  % Provides ExpMetaData

k = regexp(cs_data_path, 'Hz_out');
meta_in_path = sprintf('%s.mat', cs_data_path(1:k+1));
load(meta_in_path)   % Provides CsExpMetaIn

state_ticks = ExpMetaData.state_counts;
state_times = state_ticks*Ts;
time_total = sum(state_times);
fprintf('Total Imaging time: %.2f\n', time_total)

%%
clc
clear CsExp
width = CsExpMetaIn.width;
npix =  CsExpMetaIn.npix;


hole_depth = (20/7)*(1/1000)*(2*20);
chan_map = ChannelMap(1:5);

cs_exp = CsExp(dat, chan_map, npix, width, Ts);
cs_exp.idx_state_s


%%
f2 = figure(2+10*i); clf
ylabel('u_z')
xlabel('pixel')
ax2 = gca();
hold on;

% Since we want to convert xy-coords to pixels, move all data to the
% positive orthant:

%%
clc
clear CsExp;
% cs_exp.Gz = gk;
cs_exp.Gz = zpk([], [], 1, Ts);
cs_exp.process_cs_data(false, hole_depth);

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

f5 = figure(8 + 10*i); clf
subplot(2,3,1)
ax3 = gca();
imshow_sane(I, ax3, width, width);
title('sample');

subplot(2,3,3)
ax5 = gca();
imshow_sane(Ir_smp, ax5, width, width)
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
    % imshow_sane(PixelVectorToMatrix(Ir_bp,[n m]), ax4, width, width);
    
    imshow_sane(bp_im, ax4, width, width);

    title('BP reconstruction');
end
%
clc
s = metadata2text(ExpMetaData, Ts)
s = sprintf('%s\nperc=%.3f', s, sum(sum(pixelifsampled))/npix/npix);
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

