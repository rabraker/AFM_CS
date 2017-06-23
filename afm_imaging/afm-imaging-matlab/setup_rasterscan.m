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

lines_sec = 0.1;
sec_line = 1/lines_sec;
image_side = 15; % micro-meters.
raster_amplitude = image_side/2; 
volts2mu = 5;
mu2volts = 1/volts2mu;

% Trace is one line at sec_line. The whole period is trace and re-trace.
x_period = 2*sec_line;  
x_freq = 1/x_period;

% resolution, ie, how many lines?
pixel_size = 2;

time_total = pixel_size*x_period;

[x_rasterdata, truefreq] = raster(x_freq, Ts, time_total, 'coerce', 1);

% we want to start at 0, not -1.
x_rasterdata.Data = (x_rasterdata.Data)*raster_amplitude;
y_rasterdata = timeseries(linspace(0, image_side, length(x_rasterdata.Time))',...
                x_rasterdata.Time);

% FeedForward Gains
dc_gain = 0.65;

x_ffK = mu2volts/dc_gainx;
y_ffK = mu2volts/dc_gainx;




[y, t] = lsim(PLANT_init_x, x_rasterdata.Data*x_ffK, x_rasterdata.Time);

figure(1); clf; hold on
plot(x_rasterdata.Time, x_rasterdata.Data, t, y*volts2mu, '--')
plot(y_rasterdata.Time, y_rasterdata.Data);


%%









