% Generate x-y input waveforms.
%
% Use 0.2 hz x-dir triangle wave.
clear
clc
addpath('functions')
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
% lines_sec = 01;
% sec_line = 1/lines_sec;
% Trace is one line at sec_line. The whole period is trace and re-trace.


raster_freq = 1; % Hz.
raster_period = 1/raster_freq;

image_side = 5; % micro-meters.
raster_amplitude = image_side/2; 
% volts2mu = 5;
% mu2volts = 1/AFM.volts2mic_xy;


% resolution, ie, how many lines?
%  (eventually, this should be the same as as pixels)
npix = 128;
square_num_lines = npix;
y_height = (1/square_num_lines)*image_side;




[x_rasterdata, truefreq, points_per_line] = raster(raster_freq, Ts, raster_period-Ts, 'coerce', 1,...
                            'shift', -raster_period/4);

% we want to start at 0, not -1.
x_rasterdata.Data = (x_rasterdata.Data+1)*raster_amplitude;
y_rasterdata = timeseries(linspace(0, y_height, length(x_rasterdata.Time))',...
                x_rasterdata.Time);



D_x = tf(Ki_x, [1 -1], Ts);
H_x = feedback(D_x*PLANT_init_x, 1);
[y, t] = lsim(H_x, x_rasterdata.Data*AFM.mic2volt_xy, x_rasterdata.Time);

figure(1); clf; hold on
plot(x_rasterdata.Time, x_rasterdata.Data, t, y*AFM.volts2mic_xy, '--')
plot(y_rasterdata.Time, y_rasterdata.Data);



% We have to interleave the x & y data like
% [x(1), y(1), x(2), y(2), ....]

xy_data = zeros(2*length(x_rasterdata.Time), 1);

j = 1;
for k=1:2:length(xy_data)
    xy_data(k) = x_rasterdata.Data(j)*AFM.mic2volt_xy;
    xy_data(k+1) = y_rasterdata.Data(j)*AFM.mic2volt_xy;
    j = j+1;
end

%%
% scan_type: 0=raster, 1=CS

            
% write it to a .csv file

data_name = sprintf('raster_scan_%dpix_%dmic_%.2dHz.csv',npix,image_side, raster_freq)

target_dir = sprintf('%dmicrons', image_side)
%
data_root = PATHS.raster_image_data('5microns', 'parents')
% fullfile(getdataroot, 'raster', target_dir);
% if exist(data_root, 'file') ~=2
%     mkdir(fullfile(getdataroot, 'raster'), target_dir)
% end

data_in_path = fullfile(data_root, data_name)
meta_path = strrep(data_in_path, '.csv', '.json');

csvwrite(data_in_path, xy_data);
opts.FileName = meta_path;
            
savejson('', struct('raster_freq', raster_freq, 'npix', npix,...
              'width', image_side,...
              'points_per_line', int64(length(x_rasterdata.Time)/2),...
              'points_per_period', int64(length(x_rasterdata.Time)),...
              'total_num_points', int64(npix*length(x_rasterdata.Time)),...
              'scan_type', 0,  'fpga_input', xy_data(:)'), opts);







%%
% OLD METHOD, READING .MAT FILE IN LABVIEW
data_in_path = fullfile(data_root, data_name)
meta_path = strrep(data_in_path, '.csv', '.mat');

csvwrite(data_in_path, xy_data);
save(meta_path, 'meta');





