
frf_dat = csvread('C:\Users\arnold\Documents\labview\afm_imaging\frf_data\z_axis_frf_ontennis.csv', 1);

freqs_hz = frf_dat(2:end,1);
freqs_rad = freqs_hz/2/pi;

magDB = frf_dat(2:end,2);
mag = 10.^(magDB/20);

phaseDEG = frf_dat(2:end, 4);
phaserad = phaseDEG*pi/180;

frf = mag.*exp(1i*phaserad);

Ki = 0.01;
Dz = zpk([0], [1], Ki, 1e-5);
dz_frf = squeeze(freqresp(Dz, freqs_rad));

g = dz_frf.*frf./(1 + dz_frf.*frf);

figure(1); clf

subplot(2,1,1)
semilogx(freqs_hz, magDB)
hold on
semilogx(freqs_hz, 20*log10(abs(g)))
semilogx(freqs_hz, 20*log10(abs(dz_frf.*frf)))



subplot(2,1,2)
semilogx(freqs_hz, unwrap(phaseDEG))
hold on;
semilogx(freqs_hz, unwrap(angle(g)) );
semilogx(freqs_hz, unwrap(angle(dz_frf.*frf)))