% This script generates input data for the Swept Sines system
% identification method. A
%
% THIS SCRIPT IS SPECIFIC TO THE Z-AXIS.
%
% All data and some meta data is saved to two
% different files, which are generated off of the base name 'expName'
% variable. 
clc
clear
close all


% Some implementation notes: From an efficiency standpoint, we would like
% to be able to define, for each frequency, the number of averages and
% number of periods per average. If we do this, it leads to different
% lengths of input vectors for each frequency. Thus, there are a couple
% options:
% 1). We could forget that, figure out the MINIMUM number of samples and
% averages to take and just get a lot more data for all the other
% frequencies. This is probably the simplest route and would let me use the
% existing Labview Code.
%
% 2). We could zero pad everything. This doesn't seem like it will give
% much benifit, size wise. Maybe? Is it that big of a deal?
%
% 3). We can write each input vector line by line. This is probably the
% most efficient and perhaps the most difficult. Of course, I'm going to
% choose this. What is wrong with me?


%================================
% Control / experiment parameters
clc
clear
freq1 = logspace(log10(10), log10(180), 25)';
freq2 = logspace(log10(180), log10(300), 25)';
freq2_a = linspace(208, 230, 25)';
freq3 = logspace(log10(230), log10(10000), 200)';
input_freqs_Hz = [freq1; freq2; freq2_a; freq3];
input_freqs_Hz = forceMono(sort(input_freqs_Hz));

%%
% File parameters/experimant name
root = 'C:\Users\arnold\Documents\labview\afm_imaging\data\sys_id';
% expName = ['25-Aug-2017_exp02_ZAxis'];

sinFileName = 'z-axis_sines_in_long.csv';
sinpath = fullfile(root, sinFileName);
metaFile = fullfile(root, strrep(sinFileName,'.csv', '_metaData.csv'));

% Signal & control parameters.
Amp = .5; % Sin wave amplitude.

% Sampling frequency
Fs = 25e3; %25 Khz
Ts = 1/Fs; % Sample Time

% For labeview: PI gain (since we do this in closed loop and maximum signal
% values to keep things safe.
% Kp = 400*Ts; % PI gain.
Kp = .001;
z1 = 0.9;
ymax = 5;
umax = 8;

num = Kp;
den = [1 -1];
DD = tf(num, den, Ts);
[num, den] = tfdata(DD, 'v');


% We don't ever want to stop
Nsettle = floor(0.4/Ts);
Tsettle = Nsettle*Ts; % Amount of time to let transients die. Make it an integer multiple


NumAve         = 100;
tmin           = 0.1;
MT_min         = 2;

for i = 1:length(input_freqs_Hz)
    freq_i = input_freqs_Hz(i);
    Tsin_i = 1/freq_i;
    
    tmax    = max([tmin, NumAve*MT_min*Tsin_i]);
    if freq_i < 2000    
        MT_s(i) = ceil(tmax/NumAve/Tsin_i);
    else
        MT_s(i) = 20;
    end
end



% error checking. Make sure we didn't miss anything.
if length(MT_s) ~= length(input_freqs_Hz)
    error(['You screwed up somewhere because the length of MT is not the',  ...
           'same as length of freq_s. Need an MT defined for each frequency']);
end

% Delete the old file, since we are just appending to it. Otherwise, we
% will end up with two sets of data, one after the other inside one file.
fclose('all');
delete(sinpath)
% if this fails, use fopen('all') to see fids that haven't been closed.
%GO
h_sin = fopen(sinpath, 'w');

for i = 1:1:length(input_freqs_Hz)
    fsin = input_freqs_Hz(i);

    Tsin = 1/fsin;
    % Need total time = Tsettle + T*MT*NumAve
    wsin = 2*pi*fsin;
    
    Nsettle    = ceil(Tsettle/Ts);
    NsampsMT_i = ceil(Tsin*MT_s(i)/Ts)+1; %Number of samples per MT.
    
    N_i        = Nsettle + NsampsMT_i*NumAve;
    t          = (0:Ts:N_i*Ts)';
    
    xdk_i = Amp*sin(wsin*t);
    % write each signal to the dlmfile as a row.        
    %dlmwrite(sinpath, xdk_i', '-append', 'precision', 9);
    sf = repmat(['%f,'], 1, length(xdk_i)-1);
    sf = [sf, '%f\n'];
    
    fprintf(h_sin, sf, xdk_i);
    
    
    
    Nsamp_s(i) = length(xdk_i);
end

fclose(h_sin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create MetaData file


meta_data = zeros(length(input_freqs_Hz), 7);


meta_data(:,1) = input_freqs_Hz;
meta_data(:,2) = NumAve;
meta_data(:,3) = MT_s;
meta_data(:,4) = Nsamp_s;

TForder = length(den)-1;


metaH = fopen(metaFile,'w');
strn1  = 'Ts:, Nsettle:, ymax:, umax:, TfOrder\n';
fprintf(metaH, strn1);

strn2  = 'Freqs (Hz):, NumAve_s:,  MT_s:, NumSamp:\n'; 
fprintf(metaH, '%f, %f, %f, %f, %f\n', Ts, Nsettle, ymax, umax, TForder);

formatter = repmat('%f,', 1, TForder);
formatter = [formatter, '%f\n'];

fprintf(metaH, formatter, num);
fprintf(metaH, formatter, den);


fprintf(metaH, strn2);

for i = 1:length(input_freqs_Hz)
    fprintf(metaH, '%f, %f, %f, %f\n', meta_data(i,1), meta_data(i,2), meta_data(i,3), meta_data(i,4));
end

fclose(metaH);

% fprintf('experiment name: %s\n', expName);
