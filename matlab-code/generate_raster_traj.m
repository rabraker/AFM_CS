% Generate x-y input waveforms.
%
% Use 0.2 hz x-dir triangle wave.
clear
clc
addpath('functions')
addpath(fullfile(getMatPath(), 'dependencies', 'jsonlab'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Where to save raster data

addpath('functions/state_space_x')

plants = CanonPlants.plants_ns14(9, 1);
PLANT_init_x = plants.PLANT;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Trace is one line at sec_line. The whole period is trace and re-trace.
raster_freq = 6; % Hz.
image_side = 5; % micro-meters.
npix = 512;
xy_start_mic = [-image_side/2, -image_side/2]*0;


Ts = AFM.Ts;
% TsTicks = 1600;
Ki_x = 0.01;

raster_period = 1/raster_freq;
raster_amplitude = image_side/2; 

square_num_lines = npix;
y_height = (1/square_num_lines)*image_side;


[x_rasterdata, truefreq, points_per_line] = raster(raster_freq, Ts, raster_period-Ts, 'coerce', 1,...
                            'shift', -raster_period/4);

% we want to start at 0, not -1.
x_rasterdata.Data = (x_rasterdata.Data+1)*raster_amplitude;
x_rasterdata.Data = x_rasterdata.Data + xy_start_mic(1);

y_rasterdata = timeseries(linspace(0, y_height, length(x_rasterdata.Time))',...
                x_rasterdata.Time);
y_rasterdata.Data = y_rasterdata.Data + xy_start_mic(2);

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


% scan_type: 0=raster, 1=CS
% write it to a .json file

data_name = sprintf('raster_scan_%dpix_%dmic_%.2dHz.csv',npix,image_side, raster_freq);
target_dir = sprintf('%dmicrons', image_side);
data_root = PATHS.raster_image_data('5microns', 'parents');

data_in_path = fullfile(data_root, data_name)
meta_path = strrep(data_in_path, '.csv', '.json');
% csvwrite(data_in_path, xy_data);
opts.FileName = meta_path;

dat= struct('raster_freq', raster_freq,...
            'npix', npix,...
            'width', image_side,...
            'points_per_line', int64(length(x_rasterdata.Time)/2),...
            'points_per_period', int64(length(x_rasterdata.Time)),...
            'total_num_points', int64(npix*length(x_rasterdata.Time)),...
            'number_of_scans', 1,...
            'scan_type', 0,...
            'mu_length', image_side*2*npix,...
            'overscan_samples', 0,...
            'pre_scan_samples', 0,...
            'tip_velocity', 0,...
            'actual_sub_samble_perc', 100,...
            'fpga_input', xy_data(:)',...
            'pix_idx', int64([0,0])...
          );

savejson('', dat,  opts);


% 
% import scanning_v0.*
% meta = struct('raster_freq', raster_freq, 'npix', npix,...
%               'width', image_side);
% % OLD METHOD, READING .MAT FILE IN LABVIEW
% data_name = sprintf('raster_scan_%dpix_%dmic_%.2dHz.csv',npix,image_side, raster_freq);
% 
% target_dir = sprintf('%dmicrons', image_side);
% data_root = PATHS.raster_image_data('5microns', 'parents_csv');
% 
% data_in_path = fullfile(data_root, data_name)
% meta_path = strrep(data_in_path, '.csv', '.json');
% data_in_path = fullfile(data_root, data_name)
% meta_path = strrep(data_in_path, '.csv', '.mat');
% 
% csvwrite(data_in_path, xy_data);
% save(meta_path, 'meta');





