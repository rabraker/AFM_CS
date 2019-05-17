% This script is used to generate the information labView needs to run the
% swept sines experiment. It outputs a .csv file which will be read in by
% labview. 
%
%   .csv File Outputs:
%   ------------------
%   Ts:   digital sample period
%   Nsettle: the number of samples will throw away in post processing. Ie.,
%   how long (approximately) we figure it will take the system to reach
%   steady state.
%   Amp: The amplitude of the reference sinusoid.
%   Freqs: A list of frequencies at which we will generate sin waves.
%   NumAve_s: The number of averages to be taken at a specific frequency.
%   MT_s: The number of sin periods in each average element.
%   NumSamp: the TOTAL number of samples each sin wave should be run for. 
%
% The script will save all of this data into a .JSON file (that will be read
% into labview), which will have the following format:
%
%

clc
clear
addpath('functions')

% File name
InputFileName = 'x-axis_sines_info_first_res.json';



root = PATHS.sysid; 
metapath = fullfile(root, InputFileName);

% Signal & control parameters.
Ts = 40e-6;  % Sampling Time
Fs = 1/Ts; % Sampling frequency 25 Khz
Ts_ticks = 40e6/Fs;

Amp = 0.20; % Sin wave amplitude.


% -------------------------------------------------------------------------
% Define the frequencies we will perturb the system at. 
% from Jeff:
% load(fullfile(pwd, '../data/input_freqs.mat')); 
% input_freqs_rad
% input_freqs_Hz = input_freqs_rad/2/pi;
% kk = find(input_freqs_Hz > 5500, 1, 'first');
% input_freqs_Hz = input_freqs_Hz(1:kk);
% freq0 = linspace(0.1, 10, 20)';

% freq0 = [];
freq1 = linspace(260, 480, 40)';


freqs_Hz = [freq1];

% You shouldn't have to change any of this. 
Nsettle = floor(0.5/Ts); % This should probably be sufficient for your AFM. 
Tsettle = Nsettle*Ts; % Amount of time to let transients die. Make it an integer multiple

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --------------- Nothing below here should need to be changed ------------
M_s = zeros(length(freqs_Hz), 1);
numAve_min = 1; 
M_min = 2;    %The minimum number of sin periods in each average.
tmin   = 0.2;  % minimum amount of time to generate each sine freq.
M_max = 400;

numAve_s = ones(length(freqs_Hz),1);
clc
for i = 1:length(freqs_Hz)
    freq_i = freqs_Hz(i);
    Tsin_i = 1/freq_i;
    
    M_i = M_s(i);
    tmax_i = max(tmin, Tsin_i*M_min*numAve_s(i));
    % want tmax = numAve*M*Tsin_i;
    
    if freqs_Hz(i) >2000
      M_s(i) = 10*200;
    else
      M_s(i) =min(M_max, ceil(tmax_i/(Tsin_i*numAve_s(i))));
    end
    
    
end

[freqs_Hz, omegas, periods, Nsamps] = SSUtils.coerce_freqs_hz(freqs_Hz, M_s, Ts);
% We may have (probably will) have repeated frequencies. Delete those, and
% update the associated data.
[freqs_Hz, idx_unique] = unique(freqs_Hz); 
omegas = omegas(idx_unique);
periods = periods(idx_unique);
Nsamps = Nsamps(idx_unique);
M_s = M_s(idx_unique)
figure(2); clf
semilogx(freqs_Hz, 0*freqs_Hz, 'x')
fprintf("number of frequencies: %d\n", length(freqs_Hz))


% error checking. Make sure we didn't miss anything.
if length(M_s) ~= length(freqs_Hz)
    error(['You screwed up somewhere because the length of MT is not the',  ...
           'same as length of freq_s. Need an MT defined for each frequency']);
end

% Now, compute the number of samples required for each sin wave.
nmax = 0;
NTsWomax = 0;
Nsettle_s = zeros(length(freqs_Hz), 1);
Ncollect_s = zeros(length(freqs_Hz), 1);
Nsamp_total = zeros(length(freqs_Hz), 1);

for i = 1:1:length(freqs_Hz)
    fsin = freqs_Hz(i);
    Tsin = periods(i);
    
    % Need total time = Tsettle + T*MT*NumAve
    NsampsMT_i = (Tsin/Ts)*M_s(i);
    assert( abs(Nsamps(i) - round(NsampsMT_i, 0)) < 1e-11) 
    
%     Nsettle_i = ceil(Nsettle*Ts/Tsin);
    Nsettle_i = Nsettle;
    % Pad with extra samples, because this prevents index errors and
    % doesn't hurt anything.
    Ncollect_i = (Nsamps(i) + 5)*numAve_s(i);

    Nsamp_total(i) = Nsettle_i + Ncollect_i;
    Nsettle_s(i) = Nsettle_i;
    Ncollect_s(i) = Ncollect_i;
    
    if Nsamp_total(i) > nmax
        nmax = Nsamp_total(i);
    end
    NTsWo_i = Nsamp_total(i)*Ts*2*pi*fsin;
    if NTsWo_i > NTsWomax
       NTsWomax = NTsWo_i; 
    end
end
fprintf('max number of samples: %d\n', nmax);
fprintf('Max of sin arg: %d\n', NTsWomax);

% ------------------------------------------------------------------------%
%                      Create MetaData file
% ------------------------------------------------------------------------%

opt.FileName = metapath;
opt.FloatFormat = '%.12f';
meta_data = struct('freqs_Hz', freqs_Hz(:)', 'numAve_s', numAve_s(:)', 'num_periodsMT', M_s(:)',...
  'Nsamp_total', Nsamp_total(:)', 'Nsettle_s', Nsettle_s(:)', 'Ncollect_s', Ncollect_s(:)', ...
  'Ts', Ts, 'Amp', Amp);
savejson('', meta_data, opt)
fprintf('Copy and paste the following into LabView:\n\n')
fprintf('Input File Name: %s\n', InputFileName);
fprintf('Input Folder: %s\n', root);


