
% dat2 = csvread('C:\Users\arnold\Documents\labview\afm_imaging\data\force-map-04.csv');

dat_s{1} = csvread(fullfile(getdataroot, 'force-map', 'force-map-short_cant-01.csv'));
dat_s{2} = csvread(fullfile(getdataroot, 'force-map', 'force-map-short_cant-02.csv'));
dat_s{3} = csvread(fullfile(getdataroot, 'force-map', 'force-map-short_cant-03.csv'));
%%
% clf
k1 = 5000;
k2 = 14906;
dat = dat_s{2};
err = dat(k1:k2,1);
uz = dat(k1:k2,2);
uz = uz - uz(1);

po = err(1);
p = err - po;


volt2nm = (7/20)*1000;
figure(3)
% hold on
close all
% ---------------------- Plotting ----------------------------------
% --------------------------- Build Figure ------------------------------
figwidth = 3.4;
figheight = 3;
F1 = figure(1); clf;
set(F1, 'Units', 'Inches', 'Position', [0,0, figwidth, figheight],...
    'PaperUnits', 'Inches', 'PaperSize', [figwidth, figheight])
set(F1, 'Color', 'w');


lft = .157;
ht1 = .78;
bt1 = .15
wd = 1-lft - .01;
ax1 = axes('Position', [lft, bt1, wd, ht1])

% xlabel('tip-sample separation [nm]', 'interpreter', 'latex')
% ylabel('interaction force F(r) [N]', 'interpreter', 'latex')



plot(uz*volt2nm, p)

xlabel('$u_z$ [$\approx$ nm]', 'interpreter', 'latex')
ylabel('p = err - po [v]', 'interpreter', 'latex')



export_fig(F1, fullfile(getfigroot(), 'deflection-distance.pdf'), '-q101')


