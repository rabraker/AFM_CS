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
% XR = zeros(1, npts);
% YR = zeros(1, npts);
% How long should we sit at eache point? 
N_me = 5000;  % Probably absurdly long. 
meta_cell = repmat({N_me}, 1, npts);

% Instantiate. 
if mu
    ME = MeasEntityMu.factory([1/8000, 0]);
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
cs_data_fname = 'C:\Users\arnold\Documents\labview\afm_imaging\data\data-in2.csv';
MT.write_csv(cs_data_fname);

w = 1000*2*pi;
Gs = tf(w^2, [1 2*.7*w, w^2]);
G = c2d(Gs, Ts)

%%
% dat = csvread('C:\Users\arnold\Documents\labview\afm_imaging\data\data-out.csv');
dat = csvread('C:\Users\arnold\Documents\labview\afm_imaging\data\data-out_ontable.csv');
%%
size(dat)
% [ x(k), y(k), e_z(k),  z(k), u_z(k), FIFO_INDEX(k),  ]
% k = 636866 + 36187;
% kend = k + 5934;
% k = 100484;
k = length(dat)- 17000;
kend = length(dat)
xdat = dat(k:kend,1);
ydat = dat(k:kend,2);
err_dat = dat(k:kend,3);
uz_dat = dat(k:kend,5);
index_dat = dat(k:kend,6);

figure(2)
subplot(4,1,1)
plot(ydat)
title('y')

subplot(4,1,2)
plot(xdat)
title('x')

subplot(4,1,3)
plot(uz_dat)
ylim([-.5, 2])
title('uz')

subplot(4,1,4)
plot(err_dat)
ylim([-1, .5])
title('err')
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
