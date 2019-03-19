clc
ss_fname = fullfile(PATHS.sysid, 'z-axis_sines_info_intsamps_quickFourierCoef_1-19-2019z-drive-02.json');
ss_data = SweptSinesOnline(ss_fname);
G_frf = ss_data.FC_s(:,2)./ss_data.FC_s(:,1);
freqs = ss_data.freq_s;

ss_fname2 = fullfile(PATHS.sysid, 'ALL-axis_sines_info_intsamps_quickFourierCoef_1-19-2019z-drive-01.json');
ss_data2 = SweptSinesOnline(ss_fname2);
G_frf2 = ss_data2.FC_s(:,2)./ss_data2.FC_s(:,1);
freqs = ss_data2.freq_s;


f1 = figure(2002);clf

h1 = frf_bode_mag(G_frf, freqs, f1, 'Hz', 'k');
h2 = frf_bode_mag(G_frf2, freqs, f1, 'Hz', 'r');

% dat = load('C:\Users\arnold\Documents\matlab-bk-5-13-2018-ecepc\AFM_SS\System_Identification\data\data_ZAxis\z-axis_sines_in_329_out_9-6-2017-01.mat')
% frf_bode_mag(dat.modelFit.frf.G_frf(:), dat.modelFit.frf.freq_s, f1, 'Hz', '--g')

h1.DisplayName = 'Cantilevar A';
h2.DisplayName = 'Cantilevar B';

leg = legend([h1, h2]);
leg.Position = [0.1450 0.8349 0.2028 0.0764];

f1.Color = 'w';


%%
clc
ss_fname = fullfile(PATHS.sysid, 'ALL-axis_sines_info_intsamps_quickFourierCoef_1-20-2019z-drive-05.json');
ss_data = SweptSinesOnline(ss_fname);
G_frf = ss_data.FC_s(:,2)./ss_data.FC_s(:,1);
freqs1 = ss_data.freq_s;

ss_fname2 = fullfile(PATHS.sysid, 'ALL-axis_sines_info_intsamps_quickFourierCoef_1-17-2019z-drive-03.json');
ss_data2 = SweptSinesOnline(ss_fname2);
G_frf2 = ss_data2.FC_s(:,4)./ss_data2.FC_s(:,1);
freqs2 = ss_data2.freq_s;


f1 = figure(2002);clf

h1 = frf_bode_mag(G_frf, freqs1, f1, 'Hz', 'k');
h2 = frf_bode_mag(G_frf2, freqs2, f1, 'Hz', 'r');

% dat = load('C:\Users\arnold\Documents\matlab-bk-5-13-2018-ecepc\AFM_SS\System_Identification\data\data_ZAxis\z-axis_sines_in_329_out_9-6-2017-01.mat')
% frf_bode_mag(dat.modelFit.frf.G_frf(:), dat.modelFit.frf.freq_s, f1, 'Hz', '--g')

h1.DisplayName = 'Cantilever A';
h2.DisplayName = 'Cantilever B';

leg = legend([h1, h2]);
leg.Position = [0.1450 0.8349 0.2028 0.0764];

f1.Color = 'w';

%%
root = '/media/labserver/afm-cs/sysID/z-axis-weird/';

fn_0p02 = fullfile(root, 'ALL-axis_sines_info_intsamps_quickFourierCoef_1-21-2019-tip1-amp0p02--01.json');
fn_0p05 = fullfile(root, 'ALL-axis_sines_info_intsamps_quickFourierCoef_1-21-2019-tip1-amp0p05.json');
fn_0p2 = fullfile(root, 'ALL-axis_sines_info_intsamps_quickFourierCoef_1-21-2019-tip1-amp0p2.json');

ss_0p02 = SweptSinesOnline(fn_0p02);
ss_0p05 = SweptSinesOnline(fn_0p05);
ss_0p2 = SweptSinesOnline(fn_0p2);


G_0p02 = ss_0p02.FRF_from_FC(1,2);
G_0p05 = ss_0p05.FRF_from_FC(1,2);
G_0p2 = ss_0p2.FRF_from_FC(1,2);

figure(100); clf

h1 = frf_bode_mag(G_0p02, [], gcf, 'Hz');
h2 = frf_bode_mag(G_0p05, [], gcf, 'Hz', 'k');
h3 = frf_bode_mag(G_0p2, [], gcf, 'Hz', 'r');


h1.DisplayName = 'Amp = 0.02 [v]';
h2.DisplayName = 'Amp = 0.05 [v]';
h3.DisplayName = 'Amp = 0.2 [v]';


legend([h1,h2,h3])

%%

frf_bode_mag(G_frf2, freqs1, gcf, 'Hz', '--m');


