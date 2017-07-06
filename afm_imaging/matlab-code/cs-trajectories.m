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
N = 2500;

% State 2
[xr, yr, ind] = CSTraj.meas_entity_static(0, 0, N, 1);   
MT = CSTraj.entity2vector(xr, yr, ind);

XR = [-.5, .5, 0.25]*0;
YR = [-.4, .4, 0.15]*0;


masterTrajster = MT;
for i = 1:length(XR)
    % state 4
    [xr, yr, ind] = CSTraj.move_entity_static(XR(i), YR(i), N);   
    MT = CSTraj.entity2vector(xr, yr, ind);
    masterTrajster = [masterTrajster; MT];
    % State 2: Define this separatly, instead of just holding at the same 
    % XR, YR so that in the future this same stuff can define, eg. a
    % mu-path.
    [xr, yr, ind] = CSTraj.meas_entity_static(XR(i), YR(i), N, i+1);   
    MT = CSTraj.entity2vector(xr, yr, ind);
    
    masterTrajster = [masterTrajster; MT];
end


% masterTrajster
cs_data_fname = 'C:\Users\arnold\Documents\labview\afm_imaging\data\data-in.csv';
csvwrite(cs_data_fname, masterTrajster);

%%
dat = csvread('C:\Users\arnold\Documents\labview\afm_imaging\data\data-out.csv');
