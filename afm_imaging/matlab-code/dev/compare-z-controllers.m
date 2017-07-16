
close all
% data from agilent PID:
f_agilent = 'C:\Users\arnold\Documents\labview\afm_imaging\data\analyze-z-axis-PIon.csv';

dat_agilent = csvread(f_agilent);

% plot(dat_agilent)
% hold on


f_mine = 'C:\Users\arnold\Documents\labview\afm_imaging\data\analyze-z-axis-myPIon.csv';
dat_mine = csvread(f_mine);


t = [0:1:length(dat_mine)-1]';

ind = 65000;
plot(t(ind:end), dat_mine(ind:end))
hold on
plot(t(end-length(dat_agilent)+1:end), dat_agilent)

legend('my pid z-axis', 'agilents pid')