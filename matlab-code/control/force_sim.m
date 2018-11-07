% % F(z) = 

clear; clc
Ts = 40e-6;
rand_fname = fullfile(PATHS.sysid, 'rand_noise_zaxis_10-30-2018_01.mat');
models = load(rand_fname);
G0 = models.modelFit.G_zdir;
G_plant = G0;
% G_plant = zpk([], [0.7, 0.7], 1, Ts)
% G_plant = G_plant * dcgain(G0)/dcgain(G_plant);
KI = 0.02;
D_I = zpk([0], 1, KI, G_plant.Ts) ;
Dinv = models.modelFit.Dinv;
gdrift = models.modelFit.gdrift;
% figure, bode(G0, G_plant)
% models.modelFit.Dinv = 1;
% models.modelFit.gdrift = 1;
%%
scl = 1e-9

v2nm = (7/20)*1000;
nm2v = 0.05/20;

m2nm = 1e9
nm2m = 1e-9

% ###############################################
% Cantelevar parameters
% sig = 0.34*scl;
ro = 0.5e-9;
sig = ro*(30)^(1/6)

eps = 3.79e-22;
rho1 = 5e28;
rho2 = 5e28;
R = 50*scl;
nu1 = 0.5;
nu2 = 0.5;
E1 = 179e9;
E2 = 179e9;

k = 0.12;
% Q = 100;
Q = 2;
wo = 2*pi*5e3;
m = (k)/wo^2;

fo = R*(2/3)*pi*pi*eps*rho1*rho2*sig^4;
params = struct('sigma', sig, 'k', k, 'Q', Q, 'wo', wo,...
                'm', m, 'fo', fo, 'ro', ro);




go_den = 3*pi*( (1-nu1^2)/(pi*E1) + (1-nu2^2)/(pi*E2))
params.go = 8*sqrt(2)*sqrt(R)/(go_den);

fun = @(x) -k*x + force2(x, 0, params);
% fun = @(xz) force2(xz(1), xz(2), params);
x1_guess = ro;
x10 = fsolve(fun, x1_guess)
x20 = 0;
fun(x10)
%%
params.x0 = [-2.4220e-09; x20];

% Plot the force curve with these parameters.
% close all
r_s = [0.4:.001:5]'*scl; % m
F = r_s*0;

% params.fo = .5*scl;

for k=1:length(r_s)
    F(k) = force2(-r_s(k), 0, params);
% F(k) = fun(r_s(k));
end

figure(100)
plot(-r_s, F)
%%
% --------------- Closed Loop Sim ----------------------------------------

clc
% D_I = zpk([],[], 1, Ts);
KK = -20;
v2m = (7/20)*1e-6;
m2v = 1/v2m;
mic2nm = 1000;
m2nm = 1e9;
To = 1/10;
w_surf = 2*pi/To;

M = 4*1;

n = floor(M*To/Ts);
t = (0:n-1)'*Ts;
Amp = 10*scl;
surface = sign(sin(w_surf*t))*Amp + Amp;
surface(1) = []; t(1) = [];
figure(2)
plot(surface)

surface = timeseries(surface, t);

trun = surface.Time(end);

sim('tip_interaction_closedloop')

%%
figure(10); clf
subplot(3,1,1)
plot(uz.Time, -uz.Data*(7/20)*mic2nm)
hold on
plot(surface.Time, 20-surface.Data*m2nm);
ylabel('control p(t)')
xlabel('time [s]')
legend('control', 'surface')
ax1 = gca();
grid on;

subplot(3,1,2)
plot(simout.Time, (simout.Data(:,1))*m2nm)
ylabel('defl [nm]')
xlabel('time [s]')
% ylim([-100, 100])
ax2 = gca();
grid on;

subplot(3,1,3)
r = simout.Data(:,1) - ztot.Data;
plot(ztot.Time, r*m2nm)
title('separation r')
ax3 = gca()
grid on;

linkaxes([ax1, ax2, ax3], 'x')
% q = simout.Data(:,1);
% r = q + p_dist_control.Data;
% 
% figure(30)
% plot(r, simout.Data(:,3))
% ylabel('Force')
% xlabel('separation r [nm]')
% title('separation r')
%%
figure(31); clf
yyaxis right

% ylim([-20, 1]*1e-9)
% ylim([min(simout.Data(:,3)), max(simout.Data(:,3))]*1.1);

plot(surface.Time, surface.Data*m2nm, '-k')
hold on
plot(uz.Time, uz.Data*(7/20)*mic2nm)

yyaxis left
plot(dfl.Time, dfl.Data, 'g')


% ylim([min(surface.Data), max(surface.Data)]*1.1);
