clear

fname_p25 = '/media/labserver/afm-cs/sysID/z-axis_frf_mp25-offset.json';
fname_0 = '/media/labserver/afm-cs/sysID/z-axis_frf_zero-offset.json';

fname_ss = '/media/labserver/afm-cs/sysID/z-axis_sines_info_intsampsFourierCoef_11-17-2018-01.json';


ss_fc = SweptSinesOnline(fname_ss);

rnt_0 = RandNoiseFRF(fname_0);

rnt_p25 = RandNoiseFRF(fname_p25);


frfBode(rnt_0.G_frd, [], gcf);
frfBode(rnt_p25.G_frd, [], gcf, 'Hz', '--r');
%%

ss_fc.translate_channel_names(AFM.nametrans)
%%

G = ss_fc.FRF_from_FC(1,2);

frfBode(G, [], gcf, [], '--g');

%%
fname_sst = '/media/labserver/afm-cs/sysID/z-axis_sines_info_intsamps_out_11-17-2018-01.json'

ss_t = SweptSinesOffline(fname_sst);

%%
ss_t.translate_channel_names(AFM.nametrans);

ss_t.compute_FC_all_periods();

%%
G_ss = ss_t.FRF_from_FC(1,2);

f1 = figure(10);
f2 = figure(11);
ax1 = subplot(2,1,1);
ax2 = subplot(2,1,2);

fun = @(freq_idx) ss_t.data_view_callback(freq_idx, [ax1, ax2]);

frf_bode_dataview2({squeeze(G_ss.ResponseData)}, G_ss.Frequency,f1, 'Hz', {'r'}, fun)