clear
clc

workingRoot = 'C:\Users\arnold\Documents\labview\afm_imaging\data\sys_id\';
% fullfile(getMatPath(), 'AFM_SS', 'System_Identification');
addpath(fullfile(workingRoot, 'scripts'));
axis = 0;

root = workingRoot;

% xdkFileName = 'z-axis_sines_in_329.csv';
xdkFileName = 'z-axis_sines_in_long.csv';
metaFileName = strrep(xdkFileName, '.csv', '_metaData.csv');

xk_FileName = 'z-axis_sines_in_long_out_9-10-2017-01.csv';
% For saving results:
modFitName   = strrep(xk_FileName, '.csv', '.mat');
modFitPath   = fullfile(root, modFitName);
% xk_xdir_FileName  = [expName, '_x_dir_dataOut.csv'];


metaFile     = fullfile(root, metaFileName); 
% xdkFileName  = [expName, '_dataIn.csv'];
xdk_path     = fullfile(root, xdkFileName);
xk_xdir_path = fullfile(root, xk_FileName);

saveOn = 0;
numChan = 2;
Ni = 1; % Number of inputs.
No = 1; % Number of outputs. 


ssOpts = sweptSinesMeta('read', metaFile);
Ts = ssOpts.Ts;

% x direction input
[xdk_xdir, xk_xdir] = sys_id_raw_Data(xdk_path, xk_xdir_path, ssOpts, numChan);
xkxdk_xdir_FC1_mat_page = rawData2FourierCoeffs(xdk_xdir, xk_xdir, ssOpts, numChan);



%%
w_s = ssOpts.freq_s*2*pi;
% Do averaging
% Cross data is only accurate if we calculate it for each measurement and
% then average the results (not average first and then calc).

G_xdir = zeros(numChan+1, numChan+1, length(ssOpts.freq_s));
for kk=1:numChan+1
   for jj = 1:numChan+1
        G_xdir(kk,jj, :) = avgFreqData(xkxdk_xdir_FC1_mat_page(:,:,kk),xkxdk_xdir_FC1_mat_page(:,:,jj), ssOpts);
       
   end
end

% G is a matrix (of pages extending into the screen) that looks like this:
% where each G(i,j) is the power spectrum of G_i_j

% Current scheme (no y-data)
% [xx, x_ux,  x_xd]
% [ux_x, ux_ux, ux_xd]
% [xd_x, xd_ux, xd_xd]


% Pre-allocate:
% 4 records for frf data.
P = zeros(No, Ni, length(ssOpts.freq_s));
% 4 records for coherence data.
S = zeros(No, Ni, length(ssOpts.freq_s));
% For now, I assign first to inuitively named variables, ala Px_ux_OL_frf,
% so I can keep track of whats going on.

%[x_ux, x_uy;
% y_ux, y_uy];
% Px_ux_OL_frf    = G_xdir(1,3,:)./G_xdir(3,3,:); % Method 1??
Px_ux_OL_frf    = G_xdir(1,2,:)./G_xdir(2,2,:); % Method 1??
Px_xd_CL_ff     = G_xdir(1,3,:)./G_xdir(3,3,:); % Method 1??
P(1,1,:)        = Px_ux_OL_frf;
Sx_ux_Coherence = frfCoherence(G_xdir, 1, 2);
S(1,1,:)        = Sx_ux_Coherence;


D = tf(-.005*[1 0], [1 -1], Ts)
dfrf = squeeze(freqresp(D, w_s));
% ==========
Gx_ux_OL_frd    = frdExt(Px_ux_OL_frf, w_s, Sx_ux_Coherence, Ts); % Method 2?

G_mimo_frd      = frdExt(P, w_s, S, Ts);
% Plotting
F1 = figure(1);
plotHandles1 = frfBode(Gx_ux_OL_frd, F1, '--k');
title(plotHandles1.AX(1), 'G_{u_x,x}')


F2 = figure(2);
plotHandles1 = frfBode(squeeze(Px_xd_CL_ff), w_s/2/pi, F2, 'r');
hold on
% title(plotHandles1.AX(1), 'G_{u_x,x}')
Px_ux_OL_frf = squeeze(Px_ux_OL_frf);
frfBode(dfrf.*Px_ux_OL_frf./(1 + dfrf.*Px_ux_OL_frf), w_s/2/pi, F2, '--k');

F3 = figure(3);
frfBode(dfrf.*Px_ux_OL_frf, w_s/2/pi, F3, 'r');

%%
figure(6)

DG = dfrf.*Px_ux_OL_frf;
k = 250
plot(real(DG(1:k)), imag(DG(1:k)), 'k', real(DG(1:k)), -imag(DG(1:k)), '--k')
grid on


%%
% G_xdir(1,3,:)./G_xdir(3,3,:)
modelFit.frf.Gx_ux_OL_frf = Gx_ux_OL_frd;


modelFit.frf.G_frf        = P;
modelFit.frf.S_coherence  = S;
modelFit.frf.G_mimo_frf   = G_mimo_frd;
modelFit.frf.w_s       = w_s;
modelFit.frf.Dz        = ssOpts.Dz;
modelFit.frf.Ts        = Ts;
modelFit.frf.freq_s    = ssOpts.freq_s;


modelFit.Sux = G_xdir(1,3,:);
modelFit.Suu = G_xdir(3,3,:);

if 1
        save(modFitPath, 'modelFit');
end










