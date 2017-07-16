% This is an example of using the MasterTrajster and associated classes to
% generate CS trajectory data that is suitible for use in LabView. THe
% first example is for CS with single points. The second is for mu-paths. 


%% 
% Example 1. Discrete (static) measurement points. 
clc
pwd
addpath('classes');  % This will fail depending on where you run the script from.

% In my current LabView implementation, the movement only needs to be
% defined for one sample. I re-use that data until the stage has settled. 
N_mve = 1;

% A X-Y points at which to take measurements. 
XR = [-.5, .5, 0.25, 0];
YR = [-.4, .4, 0.15, 0];
% How long should we sit at eache point? 
N_me = 5000;  % Probably absurdly long. 
meta_cell = {N_me, N_me, N_me, N_me};

% Instantiate. 
MT = MasterTrajster(XR, YR, meta_cell, MoveEntityStatic.factory(N_mve),...
        @MeasEntityStatic);

% This is the data we will write to a csv file. We don't really need this,
% but you can inspect the data this way. 
master_data_vec = MT.as_vector(); 

% Plot. The movement is only a single point, so it doesn't show up. 
MT.visualize;

% File name for output of csv.
cs_data_fname = 'C:\Users\arnold\Documents\labview\afm_imaging\data\data-in2.csv';
MT.write_csv(cs_data_fname);

%%
% Example 2. mu-paths. 
clc

N_mve = 1;
N_me = 2000; 
XR = [-.5, .5, 0.25]*.25;
YR = [-.4, .4, 0.15]*.25;
% Since this is an example, I make all the mu-paths the same length.
% Obvisouly this is not the case in practice, which is why we have a list.
meta_cell = {N_me, N_me, N_me};

% When we generate a paramaterized MeasEntityMu, we provide a vector of
% scan-rates. This is pretty slow.
ME = MeasEntityMu.factory([1/4000, 1/4000]);
MT = MasterTrajster(XR, YR, meta_cell, MoveEntityStatic.factory(N_mve),...
        ME);

MT.visualize

cs_data_fname_mu = 'C:\Users\arnold\Documents\labview\afm_imaging\data\data-in-mu.csv';
MT.write_csv(cs_data_fname_mu)


%%
