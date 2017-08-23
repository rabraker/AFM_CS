



%          |  
% % F(z) = 

clear; clc

scl = 1e-9


ro = 0.5e-9;
sigma = ro*(30)^(1/6)

AF = atomicForces(sigma,.5*scl);

rs = [.4:.001:5]'*scl;

fs = AF.force(rs);

figure(1); clf
plot(rs, fs)
grid on

% ###############################################
% Cantelevar parameters

k = .12;;
Q = 100;
wo = 2*pi*20e3;
m = k/wo^2;
z1 = 100*scl;

cantelevar = Cantelevar(wo, Q, k, m, z1, AF)


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

