clear
load('../publications/afmMPC/mimo_modelData.mat')

G = PLANT_x;
Ts = 40e-6;
Ki_x = 0.005;
D = tf([1 0]*Ki_x, [1 -1], Ts)

H = feedback(D*G, 1)
step(H)
%%
clc
addpath('classes')
N_mve = 1;

XR = [-.5, .5, 0.25, 0];
YR = [-.4, .4, 0.15, 0];
meta_cell = {5000, 5000, 5000, 5000};
MT = MasterTrajster(XR, YR, meta_cell, MoveEntityStatic.factory(N_mve),...
        @MeasEntityStatic);

master_data_vec = MT.as_vector(); 

MT.visualize;

% masterTrajster
cs_data_fname = 'C:\Users\arnold\Documents\labview\afm_imaging\data\data-in.csv';
% dat = csvread('C:\Users\arnold\Documents\labview\afm_imaging\data\data-out.csv');
% csvwrite(cs_data_fname, masterTrajster);
csvwrite(cs_data_fname, master_data_vec);

%%
% Test out the mu-path generator.
clc
N_mve = 50;
N_me = 50;
XR = [-.5, .5, 0.25];
YR = [-.4, .4, 0.15];
meta_cell = {N_mve, N_mve, N_mve};

ME = MeasEntityMu.factory([.1, .1]);
MT = MasterTrajster(XR, YR, meta_cell, MoveEntityStatic.factory(N_mve),...
        ME);

master_data_vec = MT.asVector(); 

MT.visualize

