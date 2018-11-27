
win = window(@hann, 2^14);

y = lsim(models.modelFit.Dinv, cs_exp.uz, cs_exp.t);

[tuz_ze, f] = tfestimate(cs_exp.uz, y, win, [1000], [], 25e3);

Fig = figure(111); clf
frfBode(tuz_ze, f, Fig, 'Hz')

%%

clc

G = models.modelFit.G_zdir;
Ts = G.Ts;
[z, p] = zpkdata(models.modelFit.Dinv);
Dinv = models.modelFit.Dinv;
Ddamp = zpk(z, [0.8, 0.8]*.0, 1, Ts);
Ddamp = Ddamp/dcgain(Ddamp);
% Ddamp = 1;
Ki = -0.02;
Dki = zpk(0, [1], Ki, Ts)



H2 = minreal(Dki/(1+Dki*Dinv*G));
H = minreal(Dki/(1+Dki*Ddamp*G));

figure(1); clf
hold on
step(H, H2)

%%
figure(2); clf
bode(Dki*Ddamp*G, Dki*Dinv*G)

hold on
bode(H, H2)
grid
%%

figure(3)
%%
freqs = logspace(log10(.1), log10(12.5e3), 1000)';

rand_fname = fullfile(PATHS.sysid, 'rand_noise_zaxis_10-30-2018_01.mat');
models = load(rand_fname)

Ki_s = [0.001, 0.01, 0.04];
Ki_s = 0.02;
G = minreal(ss(models.modelFit.G_zdir)); %/models.modelFit.gdrift);
Dinv = models.modelFit.Dinv;
p = pole(G);

D2 = zpk(p(end-1:end), [0,0], 1, Ts);
D2 = D2/dcgain(D2);

F5 = mkfig(5, 5, 4); clf

% Dinv = 1;
F6 = mkfig(6, 5, 4); clf
Hu = gobjects(length(Ki_s),1);
clrs = {'b', 'r', 'k'}
for k=1:length(Ki_s)
  
  Dz = zpk([0], [1], -Ki_s(k), Ts);
  
  H_de = -minreal(1/(1+Dz*Dinv*G));
  H_du = -minreal(Dz/(1+Dz*Dinv*G));
  H_dy = minreal(Dz*Dinv*G/(1+Dz*Dinv*G));b
  frf_bode_mag(H_de, freqs, F5, 'Hz', clrs{k}, 'LineWidth', 1.5);
  hold on, grid on
  
  Hu(k) = frf_bode_mag(H_du, freqs, F6, 'Hz', clrs{k}, 'LineWidth', 1.5);
  Hu(k).DisplayName = sprintf('$K_I = %.3f$', Ki_s(k));
end

figure(F5);
title('Disturbance to error')
ylim([-50, 5])
lower = [-75, -75];
upper = [5,5];
h = ciplot(lower, upper, [7, 30],'g');
alpha(h, '.25')

figure(F6)
leg = legend(Hu);
lower = [-75, -75];
upper = [5,5];
h = ciplot(lower, upper, [7, 30],'g');
alpha(h, '.25')
h.DisplayName = 'Building Vibration Band';
set(leg, 'location', 'SouthWest', 'FontSize', 14)
title('Disturbance to control')
ylim([-50, 5])

%%
clc
G2 = absorbDelay(ss(models.modelFit.G_reduce));
figure(1)
ax1 = subplot(2,2,1); cla
ax2 = subplot(2,2,3); cla

We = zpk(.8, [(1-eps*1e8)], 0.04, Ts);
% We = 0.1*tf([4, -0.9841], [1, -], Ts);
Wu = 0.1;
W3 = zpk([0.95, .95], [0, 0], 90, Ts);
W3 = Ddamp;

frf_bode_mag(We, freqs, ax1, 'Hz');
frf_bode_mag(W3, freqs, ax1, 'Hz');



[K, CL] = mixsyn(G, We, Wu, W3);
if isempty(K)
  error('Control synthesis failed.')
end

freqs = linspace(1, 12500, 1000)';
freqs2 = models.modelFit.frf.freqs_Hz;
K_frf = squeeze(freqresp(K, freqs2*2*pi));
G_frf = models.modelFit.frf.G_uz2stage;
H_frf = (K_frf.*G_frf)./(1 + K_frf.*G_frf);

H = feedback(G*K, 1);
S = feedback(1, G*K);
frf_bode_mag(H, freqs, ax2, 'Hz');
frf_bode_mag(H_dy, freqs, ax2, 'Hz');
frf_bode_mag(H_frf, freqs2, ax2, 'Hz');

title(ax2, '$H_{yr}$')

ax4 = subplot(2,2,4); cla
title('Sens ($H_{de}$')
frf_bode_mag(S, freqs, ax4, 'Hz');
frf_bode_mag(H_de, freqs, ax4, 'Hz');

ax3 = subplot(2,2,2); cla
frf_bode_mag(K, freqs, ax3, 'Hz');