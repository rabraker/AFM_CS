% % F(z) = 

clear; clc

scl = 1e-9



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

k = 0.12
Q = 100;
wo = 2*pi*20e3;
fo = R*(2/3)*pi*pi*eps*rho1*rho2*sig^4;
params = struct('sigma', sig, 'k', k, 'Q', 100, 'wo', wo,...
                'm', k/wo^2, 'z1', 100*scl, 'fo', fo)




go_den = 3*pi*( (1-nu1^2)/(pi*E1) + (1-nu2^2)/(pi*E2))
params.go = 8*sqrt(2)*sqrt(R)/(go_den);

params.zo = sig/30^(1/6);

ell = 1*scl;
params.x0 = [0; 0]
u0 = -0*scl;


% Plot the force curve with these parameters.
% close all
r_s = [0.4:.001:5]'*scl;
F = r_s*0;

params.fo = .5*scl

for k=1:length(r_s)
    F(k) = force2(r_s(k), 0, params);
end

plot(r_s, F)


volts_per_tick = -5e-4;
uo = 100*scl;
tend = 2*-uo/volts_per_tick

dt = 40e-6;
t = linspace(0, tend, floor(tend/dt));

u = uo + t*volts_per_tick;

uz = timeseries(u, t);
trun = t(end);

Ki = 0.01;

sim('tip_interaction')

figure(1)
plot(simout.Time, simout.Data(:,1)/scl)
ylabel('defl [nm]')
% ylim([-100, 100])

figure(2)
plot(simout.Time, simout.Data(:,2))


%%
% ###############################################
% Simulation

dt = 1e-6;
t = [0:dt:4.5e-3]';

x0 = [1*scl; 0]
f = @(t, x)cantelevar.dyn(t, x, 0)
[t, x] = ode45(f, t, x0);
figure(2); clf
subplot(2,1,1)
plot(t, x(:,1))

subplot(2,1,2)
plot(t, x(:,2))

