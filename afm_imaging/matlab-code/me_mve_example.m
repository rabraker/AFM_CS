clear
clc

% The number of samples that the Moves are Defined For. 
N_mve = 5;

% The number of samples that the measurements are defined for:
N_me = 5

XR = [-.5, .5, 0.25];
YR = [-.4, .4, 0.15];

MT = MasterTrajster(XR, YR, MoveEntityStatic.factory(N_mve),...
        MeasEntityStatic.factory(N_me));

vec = MT.asVector() 
