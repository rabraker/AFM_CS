
% dat2 = csvread('C:\Users\arnold\Documents\labview\afm_imaging\data\force-map-04.csv');

% dat_s{1} = csvread(fullfile(getdataroot, 'force-map', 'force-map-short_cant-01.csv'));
% dat_s{2} = csvread(fullfile(getdataroot, 'force-map', 'force-map-short_cant-02.csv'));
% dat_s{3} = csvread(fullfile(getdataroot, 'force-map', 'force-map-short_cant-03.csv'));

dat = csvread(fullfile(getdataroot, 'force-map', 'force-map-02.csv'));
%
% clf
savefig = 1;  %binary flag to save figure.
k1 = 1;
k2 = length(dat)

% dat = dat_s{2};
err_full = dat(k1:k2,1);
uz_full = dat(k1:k2,2);
uz_full = uz_full - uz_full(1);

po = err_full(1);
p_full = err_full - po; 
%%

% easy, stupid way to downsample. Kills two birds, since I want to find the
% point where uz changes directions. Cant directly use crossing(diff(uz))
% because there are thousands of points where diff(uz) == 0, because in
% fixed point, we don't actually increment every timestep (at 16 bits).
ind = find(diff(uz_full)~=0);
uz_down = uz_full(ind); 
k3 = crossing(diff(uz_down))
err_down = err_full(ind);
p_down = p_full(ind);
% plot(uz_down, err_down)





volt2nm = (7/20)*1000;
figure(3)
% hold on
close all
% ---------------------- Plotting ----------------------------------
% --------------------------- Build Figure ------------------------------
figwidth = 3.4;
figheight = 2.5;
F1 = figure(1); clf;
set(F1, 'Units', 'Inches', 'Position', [0,0, figwidth, figheight],...
    'PaperUnits', 'Inches', 'PaperSize', [figwidth, figheight])
set(F1, 'Color', 'w');


lft = .157;
ht1 = .78;
bt1 = .15;
wd = 1-lft - .01;
ax1 = axes('Position', [lft, bt1, wd, ht1]);
hold on
grid on
% xlabel('tip-sample separation [nm]', 'interpreter', 'latex')
% ylabel('interaction force F(r) [N]', 'interpreter', 'latex')

% k3 = 6328;
% k3 = 2653330;
h1 = plot(uz_down(1:k3)*volt2nm, p_down(1:k3), 'k');
h2 = plot(uz_down(k3:end)*volt2nm, p_down(k3:end), 'r');
h1.DisplayName = 'Descent';
h2.DisplayName = 'Retraction';
leg1 = legend([h1, h2]);
set(leg1, 'interpreter', 'latex', 'Position', [0.2243 0.2515 0.2961 0.1387])
ylim([1.1*min(p_full), 1.1*max(p_full)]);

xlabel('$u_z$ [$\approx$ nm]', 'interpreter', 'latex')
ylabel('p = err - po [v]', 'interpreter', 'latex')


if savefig
    export_fig(F1, fullfile(getfigroot(), 'deflection-distance.pdf'), '-q101')
end

