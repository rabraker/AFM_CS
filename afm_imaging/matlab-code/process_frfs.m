
frf_dat = csvread('C:\Users\arnold\Documents\labview\afm_imaging\frf_data\z_axis_frf_ontennis.csv', 1);

freqs_hz = frf_dat(2:end,1);
freqs_rad = freqs_hz/2/pi;

magDB = frf_dat(2:end,2);
mag = 10.^(magDB/20);

phaseDEG = frf_dat(2:end, 4);
phaserad = phaseDEG*pi/180;

frf = mag.*exp(1i*phaserad);
%%
Ki = 0.001;
Dz1 = zpk([0], [1], Ki, 40e-6);
[~,~,~,wcp] = margin(Dz1)

Dz2 = zpk([0], [1], 1, 1e-5);

mag_wcp = abs(freqresp(Dz2, wcp))

Ki2 = 1/mag_wcp
Dz = Ki2*Dz2;


figure(10);
h = bodeplot(Dz1, Dz, {10, 10e3*2*pi});
setoptions(h, 'FreqUnits', 'Hz')
grid on


%
dz_frf = squeeze(freqresp(Dz, freqs_rad));

g = dz_frf.*frf./(1 + dz_frf.*frf);

figure(1); clf

subplot(2,1,1)
% semilogx(freqs_hz, magDB)
% hold on
% semilogx(freqs_hz, 20*log10(abs(g)))
semilogx(freqs_hz, 20*log10(abs(dz_frf.*frf)))
grid on


subplot(2,1,2)
semilogx(freqs_hz, unwrap(phaseDEG))
hold on;
semilogx(freqs_hz, unwrap(angle(g)) );
semilogx(freqs_hz, unwrap(angle(dz_frf.*frf)))
grid on
