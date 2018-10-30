% A quick script to check out the vibration isolation chamber I built. The
% first section computes the a rough estimate of the spring constant and
% natural frequency. The second section plots the new and old data with the
% cantilevar just sitting on the surface. THis is a great success, at least
% when the xy-stage is turned off. 
%
% Oddly, when the xy-stage is turned on, we get a very visible sinusoid at
% about 360 hz. I don't know what to make of this. 
% 

W = 47.4; % lbs
l = 11; % inches

L = l * 2.54 / 100;
 
lbs2newtons = 4.4822162825086;
lbs2kilos = 0.45359237; 
M_kilos = W * lbs2kilos;

F_newtons = W * lbs2newtons;

% F = k * x
K = F_newtons / L;

wo = sqrt( K / M)
fo = wo/2/pi
%%

% dat1 = csvread('/media/labserver/afm-cs/z-scope/data-out_KI_p00_nobungee_10-26-2018.csv');
dat1 = csvread('/media/labserver/afm-cs/z-scope/data-out_KI_p00.csv');
dat2 = csvread('/media/labserver/afm-cs/z-scope/data-out_KI_p00_withbungee04_10-26-2018.csv');

ze1 = detrend(dat1(:,1));
N = length(ze1);
ze2 = detrend(dat2(end-N:end,1));
length(ze2)

% N = min(length(ze1), length(ze2));

figure(8); clf

plot(ze1(1:N))
hold on;
plot(ze2(1:N))

[Z1, freqs] = power_spectrum(ze1(1:N), Ts);
[Z2, freqs] = power_spectrum(ze2(1:N), Ts);


figure(9); clf

semilogx(freqs, Z1);
hold on
semilogx(freqs, Z2)



