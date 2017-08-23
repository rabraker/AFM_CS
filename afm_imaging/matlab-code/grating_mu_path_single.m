clc
pwd
addpath('classes');  % This will fail depending on where you run the script from.

% In my current LabView implementation, the movement only needs to be
% defined for one sample. I re-use that data until the stage has settled. 
N_mve = 1;
mu = 1
% A X-Y points at which to take measurements. 
npts = 2;
XR = rand(npts)*1;
YR = rand(npts)*2;

% XR = rand(npts)*0;
% YR = rand(npts)*0;

% XR = zeros(1, npts);
% YR = zeros(1, npts);
% How long should we sit at eache point? 
N_me = 5000;  % Probably absurdly long. 
meta_cell = repmat({N_me}, 1, npts);

% Instantiate. 
if mu
    ME = MeasEntityMu.factory([1/8000, 0]);
%     ME = MeasEntityMu.factory([0, 0]);
    MT = MasterTrajster(XR, YR, meta_cell, MoveEntityStatic.factory(N_mve),...
           ME);
else
    MT = MasterTrajster(XR, YR, meta_cell, MoveEntityStatic.factory(N_mve),...
          @MeasEntityStatic);
end
% This is the data we will write to a csv file. We don't really need this,
% but you can inspect the data this way. 
master_data_vec = MT.as_vector(); 

% Plot. The movement is only a single point, so it doesn't show up. 
MT.visualize;
%%
% File name for output of csv.
cs_data_fname = 'C:\Users\arnold\Documents\labview\afm_imaging\data\data-in_single.csv';
MT.write_csv(cs_data_fname);

% w = 1000*2*pi;
% Gs = tf(w^2, [1 2*.7*w, w^2]);
% G = c2d(Gs, Ts)

%%
% dat = csvread('C:\Users\arnold\Documents\labview\afm_imaging\data\data-out.csv');
dat = csvread('C:\Users\arnold\Documents\labview\afm_imaging\data\cs-traj10-500_8-22-2017_07.csv');
% data_z_engage\z_engage_8e4.csv');

size(dat)
% [ x(k), y(k), e_z(k),  z(k), u_z(k), FIFO_INDEX(k),  ]

k = 1;
kend = length(dat)
% kend = find(dat(:,6)==1, 1, 'last')
t = [k:1:kend]'*40e-6;
% kend = k + 6493
xdat = dat(k:kend,1);
ydat = dat(k:kend,2);
err_dat = dat(k:kend,3);
uz_dat = dat(k:kend,5);
index_dat = dat(k:kend,6);

figure(5)
subplot(4,2,1); hold on
ax1 = gca;
plot(t, ydat)
title('y')

subplot(4,2,3); hold on
ax2 = gca;
plot(t, xdat)
title('x')


subplot(4,2,[5,7]); hold on
ax3 = gca;
plot(t, uz_dat)
ylim([-.5, 3])
title('uz')

% subplot(4,2,7); hold on
% ax4 = gca
% plot(t, err_dat)
% ylim([-1, .5])
% % ylim([-.07, .07])


subplot(4,2,[2,4,6,8]); hold on
ax5 = gca();
plot(t, err_dat)
ylim([-1, .5])
title('err')

linkaxes([ax1, ax2, ax3, ax4, ax5], 'x')
%%
dat = csvread('C:\Users\arnold\Documents\labview\afm_imaging\data\data-out2.csv');

size(dat)
% [ x(k), y(k), e_z(k),  z(k), u_z(k), FIFO_INDEX(k),  ]
xdat = dat(:,1);
ydat = dat(:,2);
err_dat = dat(:,3);
uz_dat = dat(:,5);
index_dat = dat(:,6);

plot(uz_dat)
