% ###############################################
% Cantelevar parameters
scl = 1e-9;  %meters to nm
addpath('force-modeling')
addpath('functions')

savefig = 1;  % Binary flag to save figure.
ro = 0.5e-9;
sig = ro*(30)^(1/6);

eps = 3.79e-22;
rho1 = 5e28;
rho2 = 5e28;
R = 50*scl;
nu1 = 0.5;
nu2 = 0.5;
E1 = 179e9;
E2 = 179e9;

k = 0.12;
Q = 100;
wo = 2*pi*20e3;
fo = R*(2/3)*pi*pi*eps*rho1*rho2*sig^4;
params = struct('sigma', sig, 'k', k, 'Q', 100, 'wo', wo,...
                'm', k/wo^2, 'z1', 100*scl, 'fo', fo);




go_den = 3*pi*( (1-nu1^2)/(pi*E1) + (1-nu2^2)/(pi*E2));
params.go = 8*sqrt(2)*sqrt(R)/(go_den);

params.zo = sig/30^(1/6);

ell = 1*scl;
params.x0 = [0; 0];
u0 = -0*scl;


% Plot the force curve with these parameters.
% close all
r_s = [0.4:.001:5]'*scl;
F = r_s*0;

params.fo = .5*scl;

for k=1:length(r_s)
    F(k) = force2(r_s(k), 0, params);
end
%%
% ---------------------- Plotting ----------------------------------
% --------------------------- Build Figure ------------------------------
clc
figwidth = 3.4;
figheight = 1.76;
F1 = figure(1); clf;
set(F1, 'Units', 'Inches', 'Position', [-10,0, figwidth, figheight],...
    'PaperUnits', 'Inches', 'PaperSize', [figwidth, figheight])
set(F1, 'Color', 'w');


lft = .157;

bt1 = .23;
ht1 = 1-bt1-.02;
wd = 1-lft - .01;

F_nm = F*1e9;
ax1 = axes('Position', [lft, bt1, wd, ht1]);
plot(r_s/scl, F_nm)
ylm = ylim;
ylim([min(F_nm)*1.3, max(F_nm)]);
xlabel('tip-sample separation [nm]', 'interpreter', 'latex');
hylab = ylabel('interaction force F(r) [nN]', 'interpreter', 'latex');
set(hylab, 'Units', 'normalized', 'Position', [-0.0666 0.36 0])
grid
%%
if savefig
%     export_fig(F1, fullfile(getfigroot(), 'force-curve.pdf'), '-q101');
    saveEps(F1, fullfile(getfigroot(), 'force-curve.eps'))
end
