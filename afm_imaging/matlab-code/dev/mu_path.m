clear; clc

matroot = getMatPath();
notepath = '/home/arnold/gradschool/proposals/nsf_afm_sept2017/latex/figs';
model_path = fullfile(matroot, 'publications', 'afmMPC', 'mimo_modelData.mat');
load(model_path);

PLANT_preT_x.InputDelay = 0;
G_x = PLANT_preT_x;
Xss  = (eye(size(G_x.B, 1)) - G_x.A)\G_x.B;
Ts = 40e-6;
Ki_x = 0.01;
D_x = ss(tf([1]*Ki_x, [1 -1], Ts));

% H = feedback(D_x*G_x, 1);
% step(H)

G_y = G_x; D_y = D_x;


% simple triangle
T=.2
%10 hz scan means 0 to 10 mu-m in 0.05 sec
mu2volts = 1/5;
volts2mu = 1/mu2volts;
tri_peak = 10; % mu-m.
Ts = G_x.Ts;
t1 = [0:Ts:T/2]';
l1 = t1*(tri_peak/t1(end));  % rise over run.
t2 = [t1(end)+Ts:Ts:T]';

l2 = l1(end)-(t2-t1(end))*(tri_peak/t1(end));

xr = [l1;l2];
t = [t1; t2];

x_ref.time = t;
x_ref.signals.values = xr;
y_ref.time = t;
y_ref.signals.values = t*(2/512)*mu2volts;

trun = t(end)

Dx_0 = 0
Dy_0 = 0;
x0 = Xss*Dx_0;
y0 = Xss*Dy_0;
FFCLI_x = zpk([], [], 1, Ts);
FFCLI_y = zpk([], [], 1, Ts);
M = 30;
rx = 0;
ry = 0;

sim('afm_tracker')

figure(10); clf
subplot(2,1,1)
plot(x_out.Time, x_out.Data*(1/mu2volts), x_ref.time, x_ref.signals.values)
legend('stage', 'x-ref')

ylabel('x [\mu m]')
subplot(2,1,2)
plot(ex_out_delayed.Time, ex_out_delayed.Data)
ylim([-0.05, 0.05])
ylabel('error [\mu m]')
xlabel ('t [s]')
grid

% saveEps(gcf(), fullfile(notepath, 'raster1.eps'))
%
% pixel rate is: micrometers/pixel

mpl = load('upathlocations.txt');

x = mpl(:,1) - mpl(1,1);
y= mpl(:,2) - mpl(1,2);;
lngs = mpl(:,3);
dirs = mpl(:,4);



pix2mu = 10/512;
% pix2mu = 1;
x = x*pix2mu;
y = y*pix2mu;
figure(2); clf
plot(x, y); hold on
xlabel('x [\mu m]')
ylabel('y [\mu m]')
title('continuous \mu-path')
% saveEps(gcf(), fullfile(notepath, 'mupath1.eps'))

figure(20); clf
plot(x, y); hold on
xlim([0,2])
ylim([0,2])
xlabel('x [\mu m]')
ylabel('y [\mu m]')
title('continuous \mu-path (zoomed in)')
% saveEps(gcf(), fullfile(notepath, 'mupath-zoom.eps'))


%%

% scanning at 10 microns per 0.05 sec, or 512 pixels per 0.05 sec.
secs_per_pix = (T/2)/512;

t_pix = linspace(0, secs_per_pix*length(x), length(x));
t_mu = [0:Ts:t_pix(end)]';

x_mu = interp1(t_pix, x, t_mu);
y_mu = interp1(t_pix, y, t_mu);

% plot(x_mu, y_mu, '--')

if 0
    figure(3); clf;
    plot(t_mu, x_mu)
end

x_ref.time = t_mu;
y_ref.time = t_mu;

x_ref.signals.values = x_mu;
y_ref.signals.values = y_mu;

trun = t_mu(end)

if 1
    H_x = feedback(D_x*G_x, 1);
    rx = order(H_x) - length(zero(H_x));
    FFCLI_x = (1/H_x)*zpk([], zeros(1,rx), 1, Ts);

    wlpf = 3000*2*pi;
    lpf = c2d(tf(wlpf, [1 wlpf])^7, Ts)
    FFCLI_x = FFCLI_x*lpf;

    %  -------------- Y ---------------------
    H_y = feedback(D_y*G_y, 1);
    ry = order(H_x) - length(zero(H_y));
    FFCLI_y = (1/H_y)*zpk([], zeros(1,ry), 1, Ts);

    FFCLI_y = FFCLI_y*lpf;

else
    FFCLI_x = zpk([], [], 1, Ts);    
    FFCLI_y = zpk([], [], 1, Ts);    
end



% The mu-path trajectory does not start (0,0). So we need an initial
% condition on the plant and we also need to preload the integrator with an
% initial condition to avoid a large initial transient. Easiest to convert
% the D_x and D_y to state space. 

% Dcgain of the systems 
DCG_x = dcgain(G_x);
DCG_y = dcgain(G_y);

% Initial steady state control inputs.
ux_0 = mu2volts*x(1)/DCG_x;
uy_0 = mu2volts*y(1)/DCG_y;

% Initial condition for the integrator which will result in ux_0 and uy_0.
Dx_0 = ux_0/D_x.C;
Dy_0 = uy_0/D_y.C;

% Initial condition for the plant.
x0  = Xss*x_ref.signals.values(1)*mu2volts/DCG_x;
y0  = Xss*y_ref.signals.values(1)*mu2volts/DCG_y;

M =11;
sim('afm_tracker')

% Plot Everything.
figure(4); clf; hold on;
plot(x_out.Data*volts2mu, y_out.Data*volts2mu);
plot(x_mu, y_mu, '--')

legend('stage', 'xy-ref')

plot(x_out.Data(1)*volts2mu, y_out.Data(1)*volts2mu, 'kx');
plot(x_mu(1), y_mu(1), 'rx')
ylabel('y-dir [\mu m]')
xlabel('x-dir [\mu m]')
ylim([0, 2])
xlim([0, 2])

% saveEps(gcf(), fullfile(notepath, 'mupath-follow.eps'))

kk = 2500;
figure(5)
subplot(2,1,1)
plot(x_out.Time(1:kk), x_out.Data(1:kk)*volts2mu,...
    x_ref.time(1:kk), x_ref.signals.values(1:kk))
ylabel('x-dir [\mu m]')
legend('stage x-dir', 'x-ref')
% subplot(3,1,2)
% plot(ux_out.Time(1:kk), ux_out.Data(1:kk))
% ylabel('u_x(t)')

subplot(2,1,2)
h1 = plot(ex_out.Time(1:kk), ex_out.Data(1:kk));
h1.DisplayName = 'Absolute error';
hold on

h2 = plot(ex_out_delayed.Time(1:kk), ex_out_delayed.Data(1:kk));
h2.DisplayName = sprintf('error w/ ref delay by %d samples', M);
legend([h1, h2])
ylim([-0.05, 0.05])
ylabel('error \mu m]')
% saveEps(gcf(), fullfile(notepath, 'mupath-time.eps'))

if 0
    figure(100);
    hold on
    plot(ux_out.Time, ux_out.Data);
end
%%
clf
kk = length(xr);
L = length(x_ref.signals.values(1:kk));
Y = fft(x_ref.signals.values(1:kk));

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = (1/Ts)*[0:(L/2)]'/L;

% df = pwelch(sin(50000*x_ref.time));
semilogx(f, 20*log10(P1))

L=length(xr);
Y = fft(xr);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = (1/Ts)*[0:(L/2)]'/L;

hold on
semilogx(f([2:2:length(f)]), 20*log10(P1([2:2:length(f)])) )


% gfrf = squeeze(freqresp(G_x, f*2*pi));
hfrf = squeeze(freqresp(feedback(G_x*D_x,1), f*2*pi));

gfrf = power_spec(impulse(G_x, [0:1:kk-1]*Ts)/2/pi, Ts);

semilogx(f, 20*log10(abs(gfrf)))
semilogx(f, 20*log10(abs(hfrf)))


%%
pix2mic = 10/512;
clc
figure(200); clf
subplot(2,1,1)
plot(xr)
hold on

[xr_filt, xr_fft_filt] = ideal_lpf(xr, 5, Ts);
ax1 = gca;
yyaxis(ax1, 'left')
% yyaxis('left')

plot(xr_filt, '--k')
xlm = xlim;

yyaxis(ax1, 'right')
plot(xr - xr_filt)
plot([xlm(1), xlm(2)], [-pix2mic, -pix2mic], '--k')
plot([xlm(1), xlm(2)], [pix2mic, pix2mic], '--k')
xr_fft = fft(xr);
power_frac = norm(xr_fft_filt(2:end))/norm(xr_fft(2:end))
%%
% ------------------
subplot(2,1,2)
ax2 = gca();
yyaxis(ax2, 'left')
xr2 = x_ref.signals.values(1:length(xr));
plot(xr2)
hold on

[xr_filt, xr_fft_filt] = ideal_lpf(xr2, 300, Ts);

plot(xr_filt, '--k')

yyaxis(ax2, 'right')
plot(xr2 - xr_filt);
power_frac = norm(xr_fft_filt)/norm(fft(xr2))

xlm = xlim;
plot([xlm(1), xlm(2)], [-pix2mic, -pix2mic], '--k')
plot([xlm(1), xlm(2)], [pix2mic, pix2mic], '--k')




% xr_fft = fft(xr);
% 
% xr_ifft = ifft(xr_fft);
% hold on
% 
% % plot(xr_ifft, '--', 'linewidth', 4)
% 
% xr_fft_filt = xr_fft;
% 
% k = find(f > 1000, 1,'first')
% kend = find(xr_fft_filt == conj(xr_fft_filt(k)), 1, 'last')
% 
% xr_fft_filt(k:kend) = 0;
% 
% xr_ifft_filt = ifft(xr_fft_filt);
% 
% plot(xr_ifft_filt, '.-k');













