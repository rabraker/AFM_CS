% Load the z-axis frequency response data.
% load('C:\Users\arnold\Documents\MATLAB\AFM_SS\System_Identification\data\data_ZAxis\z-axis_sines_in_329_out_9-3-2017-01.mat');
clear
if ispc
    dataroot = 'C:\Users\arnold\Documents\labview\afm_imaging\data\sys_id';
else
    dataroot = '/home/arnold/gradschool/publications/afm-cs/afm_imaging/data/sys_id';
end
% load('C:\Users\arnold\Documents\MATLAB\AFM_SS\System_Identification\data\data_ZAxis\z-axis_sines_in_329_out_9-6-2017-01.mat');
load(fullfile(dataroot, 'z-axis_sines_in_long_out_9-8-2017-01.mat'));

P_frf =squeeze( modelFit.frf.G_frf);
w_s = modelFit.frf.w_s;
Ts = modelFit.frf.Ts;

[P_frf, w_s] = monotonicFRF(P_frf, w_s);
freq_s = w_s/2/pi;
% freq_s = modelFit.frf.freq_s;
pmag = abs(P_frf);
phase = unwrap(angle(P_frf));
phase = phase + pi;
P_frf = pmag.*exp(j*phase);

d = 7;
Del = zpk(zeros(d,1), [], 1, Ts)
Del_frf = squeeze(freqresp(Del, w_s));

P_frf_nodelay = P_frf.*Del_frf;

% Plot it.
F1 = figure(1+10); clf
frfBode(P_frf, freq_s, F1, 'Hz', 'r')



k1 = find(freq_s < 200, 1, 'last')
k2 = find(freq_s < 240, 1, 'last')
P_fit = P_frf_nodelay(k1:k2);
freq_fit = freq_s(k1:k2);
p_frd = frd(P_fit, freq_fit*2*pi, Ts);

F2 = figure(2+10); clf
frfBode(P_fit, freq_fit, F2,  'Hz', 'r')


% Make a transfer function for the initial guess.
wz = 220.*2*pi;
wp = 214*2*pi;
zz = .0071;
zp = .0072;
D = 2.8*c2d((wp*wp/wz^2)*tf([1, 2*zz*wz, wz^2], [1, 2*zp*wp, wp^2]), Ts);

theta0 = [D.Numerator{1}, D.Denominator{1}(2:3)];
logcost_ = @(theta) logcost(theta, P_fit, freq_fit*2*pi);

opts = optimset('MaxFunEvals', 20e3, 'MaxIter', 1000, 'TolX', 1e-6);
theta = fminsearch(logcost_, theta0, opts);

G1 = tf(theta(1:3), [1, theta(4:5)], Ts);

% frfBode(D, freq_fit, F2, '--k', 'Hz');
frfBode(G1, freq_fit, F2, 'Hz', '--g');

%%


% This approximates the rolloff, but ignore it.
a = 4200*2*pi;
G2 = c2d(a*tf(1, [1, a]), Ts);

G = 2.5*zpk(G1*G2*G2*G2*G2);

G_frf = squeeze(freqresp(G, w_s));

% frfBode(G_frf,  freq_s, F1, 'k', 'Hz');
% The I controller.
K = 1000
Ki = .01;

Di = -tf(Ki*[1 0], [1 -1], Ts)
D_inv = tf(zpk((1/G1)*dcgain(G1))); %make it have unity dc gain.
G1_frf =squeeze(freqresp(G1, w_s));
frfBode(G1_frf, freq_s, F1, 'Hz', 'k');
subplot(2,1,2)
xlm = xlim;
plot(xlm, [-180, -180], 'k')

D = D_inv*Di;
Dfrf =squeeze(freqresp(D, w_s));
% and plot
F2 = figure(2+20); clf
frfBode(P_frf.*Dfrf, freq_s, F2,  'Hz', 'k');
% frfBode(P_frf.*Dfrf*.1, freq_s, F1, '--k', 'Hz');
% frfBode(P_frf.*Dfrf*.5, freq_s, F1, '--m', 'Hz');
frfBode(P_frf.*Dfrf./(1+P_frf.*Dfrf), freq_s, F2, 'Hz', 'r');
format long
D_inv.Numerator{1}
D_inv.Denominator{1}

DD = tf(zpk(D));
DD.Numerator{1};
DD.Denominator{1};
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

