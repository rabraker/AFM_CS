% Generate x-y input waveforms.
%
% Use 0.2 hz x-dir triangle wave.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

lines_sec = 0.2;
sec_line = 1/lines_sec;
image_side = 5; % micro-meters.
raster_amplitude = image_side/2; 
volts2mu = 5;
mu2volts = 1/volts2mu;

% Trace is one line at sec_line. The whole period is trace and re-trace.
x_period = 2*sec_line;  
x_freq = 1/x_period;

% resolution, ie, how many lines?
%  (eventually, this should be the same as as pixels)
number_of_lines = 20;
number_of_pixels = 256;
square_num_lines = number_of_pixels;
y_height = (1/square_num_lines)*image_side



[x_rasterdata, truefreq, points_per_line] = raster(x_freq, Ts, x_period-Ts, 'coerce', 1,...
                            'shift', -x_period/4);

% we want to start at 0, not -1.
x_rasterdata.Data = (x_rasterdata.Data +1)*raster_amplitude;
y_rasterdata = timeseries(linspace(0, y_height, length(x_rasterdata.Time))',...
                x_rasterdata.Time);

% FeedForward Gains
dc_gainx = 0.65;

x_ffK = mu2volts/dc_gainx;
y_ffK = mu2volts/dc_gainx;




[y, t] = lsim(PLANT_init_x, x_rasterdata.Data*x_ffK, x_rasterdata.Time);

figure(1); clf; hold on
plot(x_rasterdata.Time, x_rasterdata.Data, t, y*volts2mu, '--')
plot(y_rasterdata.Time, y_rasterdata.Data);


%%
% We have to interleave the x & y data like
% [x(1), y(1), x(2), y(2), ....]

xy_data = zeros(2*length(x_rasterdata.Time), 1);
j = 1;
for k=1:2:length(xy_data)
    xy_data(k) = x_rasterdata.Data(j);
    xy_data(k+1) = y_rasterdata.Data(j);
    j = j+1;
end

% write it to a .csv file
%%
data_in_path = 'C:\Users\arnold\Documents\labview\afm_imaging\data\data-in-singleperiod.csv';

csvwrite(data_in_path, xy_data);


ln



