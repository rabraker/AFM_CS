
rand_fname = fullfile(PATHS.sysid, 'rand_noise_zaxis_10-30-2018_01.mat');
models = load(rand_fname);
G0 = models.modelFit.G_zdir;
p = pole(G0); %p(2) = [];
z = zero(G0); %z(end) = [];
G2 = zpk(z, p, -1, Ts);
G2 = G2 * abs(freqresp(G0, 2*pi*10))/abs(dcgain(G2));
G_plant = G2;


Dz = models.modelFit.Dinv;
% Dz = zpk([], [], 1, Ts);
gdrift = models.modelFit.gdrift;
% gdrift = -G_plant;

KI = -0.2;
D_I = zpk([0], 1, KI, G_plant.Ts) ;

H = -ss(minreal(feedback(D_I, Dz*G_plant)));

F5 = figure(5); clf

freqs = logspace(-1, log10(12500), 1000)';
frfBode(G_plant, freqs, F5, 'Hz', '-r');
frfBode(G2, freqs, F5, 'Hz', '-g');
frfBode(models.modelFit.frf.G_uz2stage, models.modelFit.frf.freqs_Hz, F5, 'Hz', '--k');


frfBode(H, freqs, F5, 'Hz', 'b');


Ts = 40e-6;

To = 1/10;
w_surf = 2*pi/To;

M = 10*1;

n = floor(M*To/Ts);
t = (0:n-1)'*Ts;
Amp = 0.05;
surface = sign(sin(w_surf*t))*Amp;

surface = timeseries(surface, t);

trun = surface.Time(end);
sim('afm_z_simple')


figure(2); clf

plot(surface)
hold on
plot(uz_surf, '-b')

plot(uz_drift, '--k')

