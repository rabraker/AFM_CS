clear
clc
dat = csvread('C:\Users\arnold\Documents\labview\afm_imaging\data\zscope-data.csv');
dat_ol = csvread('C:\Users\arnold\Documents\labview\afm_imaging\data\zscope-data-ol.csv');
dat_ol = detrend(dat_ol);


plot(dat)

%%
load('models/modfit1.mat')
Di = Di*1;

figure(1)
margin(Di*G)
H = feedback(Di*G,1);
figure(2)
step(H)

Ts = G.Ts;
ref = 0;
t = [0:1:length(dat_ol)-1]'*Ts;
trun = t(end)
dist_dat = timeseries(dat_ol, t);

sim('models/sim_dist.slx')

%%
% plot(dat_ol)

dat_ft = fft(dat_ol);
L = length(dat_ft);
P2 = abs(dat_ft/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);



f = 25e3*(0:(L/2))/L;

semilogx(f, P1)