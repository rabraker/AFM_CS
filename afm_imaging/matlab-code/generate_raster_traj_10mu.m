% Generate x-y input waveforms.
%
% Use 0.2 hz x-dir triangle wave.
clear
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Where to save raster data
raster_root = 'C:\Users\arnold\Documents\labview\afm_imaging\data\raster';

%-------------- Location of System Model -------------------------
dataroot      = fullfile(getMatPath(), 'AFM_SS',...
                'System_Identification', 'data', 'data_xaxis'); 
expName           = ['22-Jun-2016_exp01'];
modFitName    = [expName, '.mat'];
modFitPath    = fullfile(dataroot, modFitName);
load(modFitPath, 'modelFit')


% FitNum    = 'FitNum001';
PLANT_init_x = ltiFit(modFitPath, 'SS02').sys;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Ts = 40e-6;
TsTicks = 1600;
Ki_x = 0.01;

% Hold the scan rate, in microns per sec constant to what we used before in
% the 5micron scan, because that seemed good for the z-axis BW.

% scan_velocity = 5/0.5; % 5 microns in 0.5 secs, for 1Hz raster.
scan_velocity = 5/0.1; % 5 microns in 0.1 secs, for 5Hz (T=.2) raster.
image_side = 10; % micro-meters.
raster_period = 2*(image_side/scan_velocity)
raster_freq = 1/raster_period; % Hz.


raster_amplitude = image_side/2; 
volts2mu = 5;
mu2volts = 1/volts2mu;


% resolution, ie, how many lines?
%  (eventually, this should be the same as as pixels)

number_of_pixels = 512;
square_num_lines = number_of_pixels;
y_height = (1/square_num_lines)*image_side



[x_rasterdata, truefreq, points_per_line] = raster(raster_freq, Ts, raster_period-Ts, 'coerce', 1,...
                            'shift', -raster_period/4);

% we want to start at 0, not -1.
x_rasterdata.Data = (x_rasterdata.Data +1)*raster_amplitude;
y_rasterdata = timeseries(linspace(0, y_height, length(x_rasterdata.Time))',...
                x_rasterdata.Time);



D_x = tf(Ki_x, [1 -1], Ts);
H_x = feedback(D_x*PLANT_init_x, 1);
[y, t] = lsim(H_x, x_rasterdata.Data*mu2volts, x_rasterdata.Time);

figure(1); clf; hold on
plot(x_rasterdata.Time, x_rasterdata.Data, t, y*volts2mu, '--')
plot(y_rasterdata.Time, y_rasterdata.Data);


%%
% We have to interleave the x & y data like
% [x(1), y(1), x(2), y(2), ....]

xy_data = zeros(2*length(x_rasterdata.Time), 1);
j = 1;
for k=1:2:length(xy_data)
    xy_data(k) = x_rasterdata.Data(j)*mu2volts;
    xy_data(k+1) = y_rasterdata.Data(j)*mu2volts;
    j = j+1;
end

% write it to a .csv file
%%

data_name = sprintf('raster_scan__%dpix_%dmic_%.2dHz.csv',number_of_pixels,image_side, raster_freq)
%%
data_in_path = fullfile(raster_root, data_name);

csvwrite(data_in_path, xy_data);






