% % F(z) = 

clear; clc

scl = 1e-9


ro = 0.5e-9;
% sigma = ro*(30)^(1/6)

% AF = atomicForces(sigma,.5*scl);
% 
% rs = [.4:.001:5]'*scl;
% 
% fs = AF.force(rs);
% 
% figure(1); clf
% plot(rs, fs)
% grid on

% ###############################################
% Cantelevar parameters
sig = 0.34*scl;
eps = 3.79e-22;
rho1 = 5e28;
rho2 = 5e28;
R = 50*scl;
nu1 = 0.5;
nu2 = 0.5;
E1 = 179e9;
E2 = 179e9;

params.sigma = sig;
params.k = .12;;
params.Q = 100;
params.wo = 2*pi*20e3;
params.m = params.k/params.wo^2;
params.z1 = 1000000*scl;
fo = (2/3)*pi*pi*eps*rho1*rho2*sig^4;

go_den = 3*pi*( (1-nu1^2)/(pi*E1) + (1-nu2^2)/(pi*E2))
params.go = 8*sqrt(2)*sqrt(R)/(go_den);

params.fo = fo*R;
params.zo = sig/30^(1/6);
params.x0 = [1*params.z1;0]
% cantelevar = Cantelevar(wo, Q, k, m, z1, AF)
u0 = -0*scl;

trun = .05;
sim('tip_interaction')

plot(simout.Time, simout.Data/scl)
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

