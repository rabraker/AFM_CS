close all
clear
clc
% ----------------- Scenario 1 ---------------------------------
% In this scheme, we hold constant the following:
% -- tip velocity
% -- rho == number of pixels per micron, e.g., 256/5
% 
figure(1); clf;


rho = 256/5;
% It is my experience so far that we cant scan as fast doing the mu-path.
% At 5 microns, 8% sampling etc, 5hz raster seems to produce about an
% equivalent image as 1hz mu-path, and ~204 scans, length=500 nm. Use those
% params to produce plots.

vel_mu = 5/0.5; % 1 hz, 5 microns in .5 sec
vel_raster = 5/(.1);

% raster time. Scales like width^2.
width_s = [1:.5:20]';
Traster = 2*rho*width_s.^2/vel_raster;

subplot(1,2,1)
hold on; grid on;
h0 = plot(width_s, Traster, '-k', 'LineWidth', 2);
h0.DisplayName = 'Raster time';
xlabel('width [microns]')
ylabel('scan time [s]')

% -------- mu-path ----------------------
% Assume the time to move, tip down, tip up is constant. Approximate as:



tmu = 0.18; % From experiement, seems like essentially current lower bound.
m = 0.8;  % = 204/256;

tscan = @(ell) ell./vel_mu; % time for one scan.
Tscan = @(ell, M) M.*tscan(ell); % time for all the scans.

alphas = [0.05, 0.08, 0.1, 0.15, 0.2]';
hands = [h0];
for k=1:length(alphas)
    alpha = alphas(k);
    
    
    M_s=rho*width_s*m;
    ell_s = width_s*alpha/m;
    T_mu = Tscan(ell_s, M_s) + M_s.*tmu;


  subplot(1,2,1)
  h_i = plot(width_s, T_mu);
  h_i.DisplayName = sprintf('sample %%%.2g', alpha);
  hands = [hands; h_i];
  
  subplot(1,2,2); hold on;
  hk = plot(width_s, ell_s);
  hk.Color = h_i.Color;
  
end


leg1 = legend(hands);
leg1.Location = 'northwest';

subplot(1,2,2)
grid on
xlabel('width [microns]')
ylabel('\mu-path length')

%%
% ---------- scheme 2 ------------------------
% hold mu-length/width constant.


figure(2); clf;


rho = 256/5;
% It is my experience so far that we cant scan as fast doing the mu-path.
% At 5 microns, 8% sampling etc, 5hz raster seems to produce about an
% equivalent image as 1hz mu-path, and ~204 scans, length=500 nm. Use those
% params to produce plots.

vel_mu = 5/0.5; % 1 hz, 5 microns in .5 sec
vel_raster = 5/(.1);

% raster time. Scales like width^2.
width_s = [1:.5:20]';
Traster = 2*rho*width_s.^2/vel_raster;

subplot(1,2,1)
hold on; grid on;
h0 = plot(width_s, Traster, '-k', 'LineWidth', 2);
h0.DisplayName = 'Raster time';
xlabel('width [microns]')
ylabel('scan time [s]')

% -------- mu-path ----------------------
% Assume the time to move, tip down, tip up is constant. Approximate as:



tmu = 0.18; % From experiement, seems like essentially current lower bound.

mu_len_frac = 0.5/5;

tscan = @(ell) ell./vel_mu; % time for one scan.
Tscan = @(ell, M) M.*tscan(ell); % time for all the scans.

alphas = [0.05, 0.08, 0.1, 0.15, 0.2]';
hands = [h0];
for k=1:length(alphas)
    alpha = alphas(k);
    
    
    M_s=rho*alpha*width_s/mu_len_frac;
    ell_s = width_s*mu_len_frac;
    
    T_mu = Tscan(ell_s, M_s) + M_s.*tmu;


  subplot(1,2,1)
  h_i = plot(width_s, T_mu);
  h_i.DisplayName = sprintf('sample %%%.2g', alpha);
  hands = [hands; h_i];
  
  subplot(1,2,2); hold on;
  hk = plot(width_s, ell_s);
  hk.Color = h_i.Color;
  
end


leg1 = legend(hands);
leg1.Location = 'northwest';

subplot(1,2,2)
grid on
xlabel('width [microns]')
ylabel('\mu-path length')


%%
% ------------------------- Scheme 3 ---------------------------------
% hold mu-length/width constant.
% Change from scheme 2: Instead of hold rho = pix/length constant, hold pix
% alone constant. 


figure(3); clf;


npix = 256;
% It is my experience so far that we cant scan as fast doing the mu-path.
% At 5 microns, 8% sampling etc, 5hz raster seems to produce about an
% equivalent image as 1hz mu-path, and ~204 scans, length=500 nm. Use those
% params to produce plots.

vel_mu = 5/0.5; % 1 hz, 5 microns in .5 sec
vel_raster = 5/(.1);

% raster time. Scales like width^2.
width_s = [1:.5:20]';
Traster = 2*npix*width_s/vel_raster;

subplot(1,2,1)
hold on; grid on;
h0 = plot(width_s, Traster, '-k', 'LineWidth', 2);
h0.DisplayName = 'Raster time';
xlabel('width [microns]')
ylabel('scan time [s]')

% -------- mu-path ----------------------
% Assume the time to move, tip down, tip up is constant. Approximate as:
%

gamma = (256/5)*(0.5)
npix = 256;

tmu = 0.18; % From experiement, seems like essentially current lower bound.


tscan = @(ell) ell./vel_mu; % time for one scan.
Tscan = @(ell, M) M.*tscan(ell); % time for all the scans.

alphas = [0.05, 0.08, 0.1, 0.15, 0.2]';
hands = [h0];
subplot(1,2,2)
hold on
ax2 = gca;


for k=1:length(alphas)
    alpha = alphas(k);
    
    rho_s = npix./width_s;
    ell_s = gamma./rho_s;
    
    
    M_s = alpha*(npix*width_s)./ell_s;
    
    ell_s = width_s*mu_len_frac;
    
    T_mu = Tscan(ell_s, M_s) + M_s.*tmu;


  subplot(1,2,1)
  hold on
  h_i = plot(width_s, T_mu);
  h_i.DisplayName = sprintf('sample %%%.2g', alpha);
  hands = [hands; h_i];
  
  subplot(1,2,2);
%   yyaxis left
    yyaxis(ax2, 'left')
  hk = plot(width_s, ell_s);
   hold on;
  hk.Color = h_i.Color;
  
  yyaxis(ax2, 'right')
  plot(width_s, M_s)
   hold on;
end


leg1 = legend(hands);
leg1.Location = 'northwest';

subplot(1,2,2)
grid on
yyaxis(ax2, 'left')
xlabel('width [microns]')
ylabel('\mu-path length')

yyaxis(ax2, 'right')
ylabel('number of \mu-paths')




