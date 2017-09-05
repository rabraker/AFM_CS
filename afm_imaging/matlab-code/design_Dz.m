% Load the z-axis frequency response data.
load('C:\Users\arnold\Documents\MATLAB\AFM_SS\System_Identification\data\data_ZAxis\z-axis_sines_in_329_out_9-3-2017-01.mat');

P_frf =squeeze( modelFit.frf.G_frf);
freq_s = modelFit.frf.freq_s;
w_s = modelFit.frf.w_s;
Ts = modelFit.frf.Ts;
% Plot it.
F1 = figure(1); clf
frfBode(P_frf, freq_s, F1, 'r', 'Hz')

% Fit a TF to the resonance anti-resonance pair by hand.
wz = 214.5*2*pi
wp = 212*2*pi;
zz = .011;
zp = .008;

G1 = c2d((wp*wp/wz^2)*tf([1, 2*zz*wz, wz^2], [1, 2*zp*wp, wp^2]), Ts);

% This approximates the rolloff, but ignore it.
a = 4200*2*pi;
G2 = c2d(a*tf(1, [1, a]), Ts);

% G = 2.5*zpk(G1*G2*G2*G2*G2);
G = 2.5*zpk(G1)
G_frf = squeeze(freqresp(G, w_s));

frfBode(G_frf,  freq_s, F1, 'k', 'Hz');
% The I controller. 
K = 1000
Ki = .01;

Di = tf(Ki*[1 0], [1 -1], Ts)
D_inv = tf(zpk((1/G1)*dcgain(G1))); %make it have unity dc gain.
D = -Di*D_inv;
Dfrf =squeeze(freqresp(D, w_s));

subplot(2,1,2)
xlm = xlim;
plot(xlm, [-180, -180], 'k')
% and plot
F2 = figure(2); 
frfBode(P_frf.*Dfrf, freq_s, F1, 'k', 'Hz');
%%

% Nyquist:
figure(10)
plot(real(P_frf.*Dfrf), imag(P_frf.*Dfrf));
hold on
plot(real(P_frf.*Dfrf), -imag(P_frf.*Dfrf));
grid on
%%

kk = [0:1:5000]';
t = kk*Ts;
u = 2 - (5e-5)*kk;
figure(50)
plot(t, u)

y = lsim(D_inv, u, t);
plot(t, y)


