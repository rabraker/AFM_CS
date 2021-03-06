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


