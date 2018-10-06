clear

load('C:\Users\arnold\Documents\MATLAB\AFM_SS\System_Identification\data\data_multiAxis\13-Sep-2017_exp01_multiAxis.mat')

modelFit_nomass = modelFit;

load('C:\Users\arnold\Documents\MATLAB\AFM_SS\System_Identification\data\data_multiAxis\13-Sep-2017_exp02_multiAxis.mat')
modelFit_mass = modelFit;


G_frf_no_mass = modelFit_nomass.frf.G_frf;
freqs = modelFit_nomass.frf.freq_s;

G_frf_mass = modelFit_mass.frf.G_frf;
freqs_mass = modelFit_mass.frf.freq_s;


F1 = figure(1);
frfBode(squeeze(G_frf_no_mass(1,1,:)), freqs, F1, '-k', 'Hz');

frfBode(squeeze(G_frf_mass(1,1,:)), freqs_mass, F1, '--r', 'Hz');


F2 = figure(2);
frfBode(squeeze(G_frf_no_mass(2,2,:)), freqs, F2, '-k', 'Hz');

frfBode(squeeze(G_frf_mass(2,2,:)), freqs_mass, F2, '--r', 'Hz');


