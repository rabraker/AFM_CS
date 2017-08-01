clc
pwd
addpath('classes');  % This will fail depending on where you run the script from.

% In my current LabView implementation, the movement only needs to be
% defined for one sample. I re-use that data until the stage has settled. 
N_mve = 1;
mu = 1
% A X-Y points at which to take measurements. 
npts = 20;
XR = rand(npts)*2;
YR = rand(npts)*2;
% XR = zeros(1, npts);
% YR = zeros(1, npts);
% How long should we sit at eache point? 
N_me = 8000;  % Probably absurdly long. 
meta_cell = repmat({N_me}, 1, npts);

% Instantiate. 
if mu
    ME = MeasEntityMu.factory([1/5000, 0]);
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

% File name for output of csv.
cs_data_fname = 'C:\Users\arnold\Documents\labview\afm_imaging\data\data-in2.csv';
MT.write_csv(cs_data_fname);

%%
dat = csvread('C:\Users\arnold\Documents\labview\afm_imaging\data\data-out5.csv');

size(dat)
% [ x(k), y(k), e_z(k),  z(k), u_z(k), FIFO_INDEX(k),  ]
xdat = dat(:,1);
ydat = dat(:,2);
err_dat = dat(:,3);
uz_dat = dat(:,5);
index_dat = dat(:,6);

figure(2)
subplot(2,1,1)
plot(uz_dat)
% ylim([-.7, -.2])

subplot(2,1,2)
plot(err_dat)
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
