clear, clc
addpath('functions')
Ts = 40e-6;
Ki_s = [0.005, 0.01, 0.05];
LPF = zpk([], [0.7, 0.7], 1, Ts);
LPF = LPF/dcgain(LPF);
models = load(fullfile(PATHS.sysid, 'x-axis_sines_infoFourierCoef_10-21-2018-03.mat'));
G = balreal(-models.modelFit.G_zdir);
sys_og = balreal(-models.modelFit.G_zdir_red);

% sys = ss(zpk([], 0.95, 1-0.95, Ts));
% G = sys;
Ki_s = [0.005, 0.01, 0.05];


root = '/media/labserver/afm-cs/z-scope';


% fnames = {'data-out_KI_p005.csv',...
%   'data-out_KI_p01.csv',...
%   'data-out_KI_p1.csv',...
%   'data-out_KI_p00.csv',...
%   'data-out_KI_p05.csv'};
f2 = 'data-out_KI_p01.csv';
fname = 'data-out_KI_p1.csv';
dat_k = csvread(fullfile(root, f2));
dat_k = csvread(fullfile(root, fname));

ze_k = detrend(dat_k(:,1));
uz_k = detrend(dat_k(:,2));
t_k = (0:length(ze_k)-1)'*Ts;

z = zero(sys_og)
p = pole(sys_og)
p_eject = p(end);
z_eject = z(end);
g_eject = zpk(p_eject, z_eject, (1-z_eject)/(1-p_eject), Ts)

p_eject2 = p(end-2:end-1);
z_eject2 = z(end-2:end-1);
g_eject2 = zpk(p_eject2, z_eject2, 1, Ts);
g_eject2 = g_eject2/dcgain(g_eject2)

sys = minreal(g_eject2*g_eject*sys_og);
G = sys;
%%
[Z_fft, freqs] = power_spectrum(ze_k, Ts);
figure(150)
semilogx(freqs, Z_fft)


%%

dist = timeseries(ze_k, t_k);
n = length(t_k);
% surface = zeros(n,1); 
% surface(n/2:end) = 0.05*0;
w_surf = 5*2*pi;
surface = sign(sin(w_surf*t_k))*0.05;

dist.Data = surface + sin(29.8*2*pi*t_k)*0.01; %dist.Data;
ref = timeseries(surface, t_k);
wo = 29.8*2*pi;

% Aw = [0, 1; %, 0; 
%   -wo^2, 0]; %, 0;
%     %0,   0  0];
% Hw = [1, 0]; % 1]
% Gw = ss(Aw, [0;0], Hw, 0);
% Gwz = c2d(Gw, Ts);
% Cdd = [sys.c*0, 1,0];

Aw = [0, 1, 0; 
  -wo^2, 0, 0;
    0,   0  0];
Hw = [1, 0, 1]
Gw = ss(Aw, [0;0; 0], Hw, 0);
Gwz = c2d(Gw, Ts);
Cdd = [sys.c*0, 1,0,0];

% Aw = 1;
% Hw = 1;
% Bw = 0;
% Gwz = ss(Aw, Bw, Hw, 0, Ts);
% Cdd = [sys.c*0, 1];

ns1 = size(sys.b,1);
nsw = size(Gwz.b, 1);

a = [sys.a, sys.b*0*Gwz.c;
             zeros(nsw, ns1), Gwz.a];
b = [sys.b; 0*Gwz.b];
c = [sys.c, Gwz.c];
sys_obs = ss(a, b, c, 0, Ts);

% Lx = dlqr(sys.a', sys.c', eye(ns1), 10000)';
% p_int = [0.9876+0.002844j; 0.9876-0.002844j; 0.98];
% % p_int = 0.98;
% pdes = [eig(sys.a - Lx*sys.c); p_int];
% pdes = getCharDes(sys_obs, [1.5,1,1], 0.99, [0.8,.7, .7], [1, 1])
pdes = getCharDes(sys_obs, [1.,1], 0.9, [0.9, 0.7], [1]);
Lz = place(a', c', pdes)';
% Lz = dlqr(sys_obs.a', sys_obs.c', eye(ns1+nsw)*1, 10000)';
% p1 = eig(a - Lz2*c);
% p2 = eig(a - Lz*c);
% 
% figure(3), clf
% plot(real(p1), imag(p1), 'xk')
% hold on
% plot(real(p2), imag(p2), 'xr')
% zgrid

Dz = zpk([0], 1, 0.05, Ts);
I = eye(ns1);
K = dlqr(sys.a, sys.b, sys.c'*sys.c, 1);

Nbar = -1; %SSTools.getNbar(sys, K);

trun = t_k(end)

x0_obs = zeros(ns1+nsw,1);

M = [I, zeros(ns1, nsw)];
ens1 = zeros(1, ns1+nsw);
ens2 = zeros(1, ns1+nsw);
ens1(ns1+1) = 1;
ens2(end) = 1;



sim('afm_z_dist_sin_est')

figure(10); clf
subplot(3,2,1)
plot(dist)
hold on
plot(surf_est+sin_est, '--')
title('30hz sin estimate')
% legend('truth', 'estimate')
grid
subplot(3,2,3)
hold on
plot(y_hat1)
title('y-hat')
grid

subplot(3,2,5)
plot(surf_est.Time, surf_est.Data)
hold on
plot(uz.Time, surface)
title('surface')
grid

subplot(3,2,2)
plot(ze)
hold on
title('ze')

subplot(3,2,4)
plot(uz)
title('uz')

% figure(200), clf
% hold on
% plot(uz.Time, uz.Data, '--')
%%
Cc = Cdd*0; Cc(end) = 1;
[Sens, H_dc_d, H_ud] = loops(G, Gwz, Lz, K, Cc, sys_obs)
%
bode(H_dc_d)

% %%
% a_t = a - Lz*c;
% c_t = [sys.c*0, Gwz.c];
% H1 = ss(a_t, b, c_t, 0, Ts);
% H2 = ss(a_t, Lz, c_t, 0, Ts);
% 
% H_de = minreal(  (H2 -1)/( 1+ Dz*(G-H1-H2*G)));
% 
% figure(7)
% 
% bode(minreal(H_de*Dz))
% u = lsim(Dz*H_de, surface.Data, surface.Time);
% figure(200), clf
% plot(surface.Time, u)
% %%

function [Sens, H_dc_d, H_ud] = loops(G, Gw, LL, K, Cc, sys_obs)
  ns = size(G.a,1);
  nd = size(Gw.a,1);
  
  Lx = LL(1:ns);
  Ld = LL(ns+1:end);
  Nbar = SSTools.getNbar(G, K);
  [A, B, C] = ssdata(G);
  [Ad, Bd, Cd] = ssdata(Gw);
  
  At = [A - B*K - Lx*C,  -Lx*Cd - B*Nbar*Cd;
       -Ld*C,             Ad - Ld*Cd];
  Lt = [Lx; Ld];
  Kt = [K, Cd*Nbar];
  
  g_1 = ss(At, Lt, Kt, 0, G.Ts);
  Loop = G*g_1;
  
  Sens = minreal(1/(1 + Loop));
  
  H_ud = minreal(-g_1*Sens);
  
  
  g_2 = ss(At, Lt, Cc, 0, G.Ts);
  
  H_dc_d = minreal(g_2*Sens);
  
end





