% Load time data and compute PSDs, looking at what happens when the fan is
% connected, or not.
clear
clc
fig_root = '~/matlab/afm-cs/matlab-code/notes/figures';
root = '/media/labserver/afm-cs/z-scope/psd-test-fan';
root_nb = '/media/labserver/afm-cs/z-scope/psd-test-fan/9-18-2018';

fname_ext_psd = 'xz_engaged_frf_externalpower.json';
fname_int_psd = 'xz_engaged_frf_internalpower.json';
fname_nobung_psd = 'z_zengaged_xyon_nobungee_frd.json';


fpath_int = fullfile(root, fname_int_psd);
fpath_ext = fullfile(root, fname_ext_psd);
fpath_nb = fullfile(root_nb, fname_nobung_psd);

psd_ext = LvPSD(fpath_ext);
psd_int = LvPSD(fpath_int);
psd_nobung = LvPSD(fpath_nb);

% Load the time domain data
fname_ext_tm = 'xz_engaged_time_externalpower.json';
fname_int_tm = 'xz_engaged_time_internalpower.json';
fname_nobung_tm = 'z_zengaged_xyon_nobungee_frd_time.json';

fpath_int_tm = fullfile(root, fname_int_tm);
fpath_ext_tm = fullfile(root, fname_ext_tm);
fpath_nobung_tm = fullfile(root_nb, fname_nobung_tm);


rnt_ext = LvPSDTime(fpath_ext_tm);
rnt_int = LvPSDTime(fpath_int_tm);
rnt_nb = LvPSDTime(fpath_nobung_tm);
%%
% -------------------- Step 0 ---------------------------------------------
% figure(10)
clc
freqs = psd_ext.freqs;

real_idx = 30;
end_idx = 21799;
z_idx = 3;
fan_idx = 4;
idx1 = 2486:3143;

tvec = (0:1:rnt_ext.num_points-1)'*rnt_ext.Ts;
% tend = tvec(end_idx);
% fan_v_int = detrend(rnt_int.data(fan_idx, :, real_idx));
% fan_v_ext = detrend(rnt_ext.data(fan_idx, :, real_idx));

F10 = mkfig(10, 5, 5.35, true); clf
ax10 = axes('Units', 'inches', 'Position', [.6, .5, 4.25 4.5]);
ax10.Color = fig_color;
ht1 = plot(tvec, detrend(rnt_nb.data(z_idx, :, real_idx)));
hold on
ht2 = plot(tvec, detrend(rnt_int.data(z_idx, :, real_idx)));

ht3 = plot(tvec, detrend(rnt_ext.data(z_idx, :, real_idx)));
title('deflection', 'FontSize', 14)
ylim([-0.07, 0.04])
xlim([tvec(1), tvec(end_idx)])
hold on
xlabel('time [s]', 'FontSize', 14)
ylabel('deflection [v]', 'FontSize', 14)
ylm = ylim;
grid on

F11 = mkfig(11, 5, 5.35, true); clf
% ax33 = axes( 'Units', 'inches', 'position', [4.8479 3.6354 3.4496 2.6146]);
ax11 = axes('Units', 'inches', 'Position', [.6, .5, 4.35 4.5]);

% hold on, grid on, ylim(ylm); xlim(xlm);

h1 = semilogx(freqs, psd_nobung.psd(:, z_idx));
hold on, grid on,
h2 = semilogx(freqs, psd_int.psd(:, z_idx));

h3 = semilogx(freqs, psd_ext.psd(:, z_idx)); 
ylim([-110, -40])
xlim([4.4729, 3000])
h1.DisplayName = 'Step 0: un-modified';
h2.DisplayName = 'Step 1: isolate from building';
h3.DisplayName = 'Step 1: + electrical isolation';

ht1.DisplayName = 'Step 0: un-modified';
ht2.DisplayName = 'Step 1: isolate from building';
ht3.DisplayName = 'Step 1: + electrical isolation';


ax11.Color=fig_color;
title('z-deflection')
lower = [-75, -75];
upper = [-40, -40];
h = ciplot(lower, upper, [7, 75],'g');
alpha(h, '.25')

xlabel('Frequency [Hz]', 'FontSize', 14);
ylabel('power spectral density')

leg1 = legend([h1, h2, h3]);
set(leg1, 'FontSize', 14, 'Position', [0.1361 0.1586 0.5938 0.1504], 'box', 'off')

leg2 = legend([ht1, ht2, ht3]);
set(leg2, 'FontSize', 14, 'Position', [0.1361 0.1197 0.5938 0.1504], 'box', 'off')

% ax11_2 = axes('Units', 'inches', 'Position', [3.0208 2.8854 1.8979 2.0833],...
%   'XTick', [], 'YTick', []);

% h2.Visible = 'off';
% ht2.Visible = 'off';
% h3.Visible = 'off';
% ht3.Visible = 'off';
% % h.Visible = 'off';
% leg1.Visible = 'off';


% h2.Visible = 'on';
% ht2.Visible = 'on';
% h.Visible = 'off';
leg1 = legend([h1, h2, h3]);
set(leg1, 'FontSize', 14, 'Position', [0.1361 0.0963 0.5938 0.1504], 'box', 'off')

save_fig(F10, fullfile(fig_root, 'dfl_time'), false);
save_fig(F11, fullfile(fig_root, 'dfl_psd'), false);

%%
lower = [-65, -65];
upper = [-40, -40];
h = ciplot(lower, upper, [300, 450],'g');
alpha(h, '.25')

leg1 = legend([h1, h2]);
set(leg1, 'FontSize', 14, 'Position', [0.1361 0.1586 0.5938 0.1504], 'box', 'off')

print(F10, '-dpdf', '~/gradschool/rob_ss_prez/figures/initial_dfl_time2.pdf')
print(F11, '-dpdf', '~/gradschool/rob_ss_prez/figures/initial_dfl_psd2.pdf')

%%
h3.Visible = 'on';
ht3.Visible = 'on';
h.Visible = 'off';

leg1 = legend([h1, h2, h3]);
set(leg1, 'FontSize', 14, 'Position', [0.1361 0.1586 0.5938 0.1504], 'box', 'off')
lower = [-65, -65];
upper = [-40, -40];
h = ciplot(lower, upper, [300, 450],'g');
alpha(h, '.25')
leg1 = legend([h1, h2, h3]);
set(leg1, 'FontSize', 14, 'Position', [0.1361 0.1586 0.5938 0.1504], 'box', 'off')

print(F10, '-dpdf', '~/gradschool/rob_ss_prez/figures/initial_dfl_time3.pdf')
print(F11, '-dpdf', '~/gradschool/rob_ss_prez/figures/initial_dfl_psd3.pdf')

%%
figure(1); clf

power_line_freqs = 60*[1:100];

subplot(2,1,1)


freqs = psd_ext.freqs;

semilogx(psd_ext.freqs, (psd_ext.psd))
ylm = [-120, -40];
ylim(ylm)
xlim([1, 12500])
hold on
grid on
legend('x-sensor', 'y-sensor', 'z-deflection', 'nPoint Fan voltage supply')
subplot(2,1,2)
% semilogx([power_line_freqs; power_line_freqs], ylm, 'LineStyle', '--', 'Color', [.5,.5,.5])
% hold on
semilogx(psd_int.freqs, (psd_int.psd) )
ylim(ylm)
xlim([1, 12500])

grid
legend('x-sensor', 'y-sensor', 'z-deflection', 'nPoint Fan voltage supply')
%%
freqs = psd_ext.freqs;

idx_x = 1;
idx_y = 2;
idx_z = 3;
idx_fan = 4;

ylm = [-120, -40];
xlm = [10, 12500];

% figure(3); 
F3 = mkfig(3, 8.5, 6.5); clf
% subplot(2,2,1)
ax11 = axes('Units', 'inches', 'Position', [0.7000 3.6354 3.4496 2.6146]);
h1 = semilogx(freqs, psd_int.psd(:, idx_x));
hold on, grid on, ylim(ylm); xlim(xlm);
h2 = semilogx(freqs, psd_ext.psd(:, idx_x));
% h3 = semilogx(freqs, psd_nobung.psd(:, idx_x));

title('x-sensor');
ylabel('[dB]');
h1.DisplayName = 'with nPoint fan voltage';
h2.DisplayName = 'with external fan voltage';
% h3.DisplayName = 'no-bungee';

lower = [-120, -120];
upper = [-40, -40];
h = ciplot(lower, upper, [327, 437],'g');
alpha(h, '.25')



% subplot(2,2,3)
ax22 = axes('Units', 'inches', 'Position', [0.7000 0.5521 3.4496 2.6146]);
semilogx(freqs, psd_int.psd(:, idx_y));
hold on, grid on, ylim(ylm); xlim(xlm);
semilogx(freqs, psd_ext.psd(:, idx_y));
% semilogx(freqs, psd_nobung.psd(:, idx_y));
lower = [-120, -120];
upper = [-40, -40];
h = ciplot(lower, upper, [327, 437],'g');
alpha(h, '.25')

title('y-sensor');
xlabel('frequency [Hz]');
ylabel('[dB]');

% subplot(2,2,2)
ax33 = axes( 'Units', 'inches', 'position', [4.8479 3.6354 3.4496 2.6146]);
semilogx(freqs, psd_int.psd(:, idx_z))
hold on, grid on, ylim(ylm); xlim(xlm);
semilogx(freqs, psd_ext.psd(:, idx_z));
% semilogx(freqs, psd_nobung.psd(:, idx_z));
title('z-deflection')
lower = [-120, -120];;
upper = [-40, -40];
hz = ciplot(lower, upper, [327, 437],'g');
alpha(hz, '.25')

ax55 = axes( 'Units', 'inches', 'position', [7.02083 5.01250 1.2812 1.2500]);
semilogx(freqs, psd_int.psd(:, idx_z))
hold on, grid on, ylim(ylm); xlim(xlm);
semilogx(freqs, psd_ext.psd(:, idx_z));
xlim([327, 437]);
ylim([-90, -40]);
set(ax55, 'XTickLabel', [], 'YTickLabel', []);
lower = [-120, -120];
upper = [-40, -40];
h = ciplot(lower, upper, [327, 437],'g');
alpha(h, '.25')


% subplot(2,2,4)
ax44 = axes('Units', 'inches', 'position', [4.8479 0.5521 3.4496 2.6146]);
semilogx(freqs, psd_int.psd(:, idx_fan));
hold on, grid on, ylim(ylm); xlim(xlm);
semilogx(freqs, psd_ext.psd(:, idx_fan));
xlabel('frequency [Hz]');
title('nPoint fan voltage');
lower = [-95, -95];;
upper = [-40, -40];
h = ciplot(lower, upper, [327, 437],'g');
alpha(h, '.25')


leg1 = legend([h1, h2]);
set(leg1, 'FontSize', 12)

print(F3, '-dpdf', '~/matlab/website/fan_noise_psd.pdf')
print(F3, '-dpdf', '~/gradschool/rob_ss_prez/figures/fan_noise_psd.pdf')
%%


%%
clc
real_idx = 30;
end_idx = 21799;
z_idx = 3;
fan_idx = 4;
idx1 = 2486:3143;
tvec = (0:1:rnt_ext.num_points-1)'*rnt_ext.Ts;
tend = tvec(end_idx);
fan_v_int = detrend(rnt_int.data(fan_idx, :, real_idx));
fan_v_ext = detrend(rnt_ext.data(fan_idx, :, real_idx));

f2 = mkfig(2, 8, 6); clf
% ax1 = subplot(2,1,1);
wd = 3.25;
ht1 = 3.25;
ax1 = axes('Units', 'inches', 'Position', [0.6458 2.5 wd ht1]);
            
h1 = plot(tvec, fan_v_int);
xlim([0, tvec(end)])
hold on
h2 = plot(tvec, fan_v_ext);
h1.DisplayName = 'with internal fan voltage source';
h2.DisplayName = 'with external fan voltage source';
xlim([0, tvec(end_idx)])
ylim([-0.035, 0.015])
title('fan voltage (deviation from mean)')
ylabel('fan voltage')
xlabel('time [s]')
ylm = ylim;
lower = [ylm(1), ylm(1)];
upper = [ylm(2), ylm(2)];

h = ciplot(lower, upper, [.1, .125],'g');
alpha(h, '.25')
leg1 = legend([h1, h2]);
set(leg1, 'FontSize', 12, 'Position', [0.1291 0.8845 0.3583 0.0745])


ax3 = axes('Units', 'inches', 'Position', [0.6458    .2564    wd    1.6]);
plot(tvec(idx1), fan_v_int(idx1))
hold on
plot(tvec(idx1), fan_v_ext(idx1))
title(ax3, 'fan voltage detail')
xlim([tvec(idx1(1)), tvec(idx1(end))])
ylm = ylim;
lower = [ylm(1), ylm(1)];
upper = [ylm(2), ylm(2)];
h = ciplot(lower, upper, [tvec(idx1(1)), tvec(idx1(end))],'g');
alpha(h, '.25')
% set(ax3, 'XTickLabel', [])

a1 = annotation('doublearrow', [0.1393 0.1719], [0.0937 0.0937]);
s1 = sprintf('To = %.5f s, f=%.1f Hz', 0.0026, 1/0.0026)
text([0.103], [-0.035], s1)



ax2 = axes('Units', 'inches', 'Position', [4.7 2.5 wd ht1]);
% h0 = plot(tvec, detrend(rnt_nb.data(z_idx, :, real_idx)));
hold on
h3 = plot(tvec, detrend(rnt_int.data(z_idx, :, real_idx)));

h4 = plot(tvec, detrend(rnt_ext.data(z_idx, :, real_idx)));
title('deflection')
ylim([-0.025, 0.025])
xlim([tvec(1), tvec(end_idx)])
hold on
xlabel('time [s]')
ylabel('deflection [v]')
ylm = ylim;
lower = [ylm(1), ylm(1)];
upper = [ylm(2), ylm(2)];

h = ciplot(lower, upper, [.1, .125],'g');
alpha(h, '.25')


ax4 = axes('Units', 'inches', 'Position', [4.7 .2564, wd, 1.6]);
dat1 = detrend(rnt_int.data(z_idx, :, real_idx));
dat2 = detrend(rnt_ext.data(z_idx, :, real_idx));
h3 = plot(tvec, dat1);
hold on
h4 = plot(tvec, dat2);
title('deflection')
ylim([-0.025, 0.025])

xlim([tvec(idx1(1)), tvec(idx1(end))])
ylm = ylim;
lower = [ylm(1), ylm(1)];
upper = [ylm(2), ylm(2)];
h = ciplot(lower, upper, [tvec(idx1(1)), tvec(idx1(end))],'g');
alpha(h, '.25')
set(ax3, 'XTickLabel', [])
title('deflection detail')


print(f2, '-dpdf', '~/matlab/website/fan_noise_time.pdf')

%%















%%
clc
Ts = rnt_ext.Ts;
N = rnt_ext.num_realizations;
Npts = rnt_ext.num_points;

% psd_fan = zeros(ceil(Npts/2), rnt_fan.number_outputs);
psd_fan = 0;
psd_nofan = 0;
psd_batteryfan = 0;

for k=1:N
  df_k = rnt_ext.data(3, :, k)';
  df_k = df_k - mean(df_k);
  win = window(@hann, Npts);
  [psd_k, freqs] = power_spectrum(df_k, Ts, win);
  

%   [psd_k, freqs] = pwelch(df_k, [], 1, length(f)*2, 1/Ts, [], f, 'onesided', 'psd');
% [psd_k, freqs] = pwelch(df_k, Npts, [], f, 1/Ts, 'onesided', 'psd');

%   dnf_k = rnt_nofan.data(3, :, k)';
%   dnf_k = dnf_k - mean(dnf_k);
%   [psd_nofan_k, freqs] = pwelch(dnf_k, Npts, [], [], 1/Ts, 'onesided', 'Power');
% 
%   
%   dbtf_k = rnt_batteryfan.data(3, :, k)';
%   dbtf_k = dbtf_k - mean(dbtf_k);
%   [psd_batteryfan_k, freqs] = pwelch(dbtf_k, Npts, [], [], 1/Ts, 'onesided', 'Power');
% 
  psd_fan = psd_fan + psd_k;
%   psd_nofan = psd_nofan + psd_nofan_k;
%   psd_batteryfan = psd_batteryfan + psd_batteryfan_k;
%   
%   
end

freqs = freqs(2:end);
psd_est = (psd_fan(2:end)/N);

figure(2); clf;

semilogx(psd_ext.freqs, (psd_ext.psd(:,3)))
hold on

semilogx(freqs, 10*(log10(psd_est)) )

% semilogx(freqs, 10*log10(psd_nofan/N))
% semilogx(freqs, 10*log10(psd_batteryfan/N))

% dat = csvread(fullfile(root, 'lv_psd_data.csv'),1);
%%

x = squeeze(rnt_ext.data(3, :, 1))';
% win = rectwin(length(x));
window(@hann, Npts);
x_dft1 = power_spectrum(x, 1, win);

x_dft2 = periodogram(x, win, length(x), 1);

figure(2); clf;
plot(10*log10(x_dft1))
hold on
plot(10*log10(x_dft2))

%%

figure
plot(squeeze(rnt_ext.data(3, :, 1)))

hold on
plot(squeeze(rnt_int.data(3, :, 1)))


%%
psd_nofan = [];
psd_batteryfan = [];
start_idx = 3;
for k=1:rnt_ext.number_outputs
  
  datf_chan_k = squeeze(rnt_ext.data(k, :, :));
  datnf_chan_k = squeeze(rnt_int.data(k, :, :));
  datbtf_chan_k = squeeze(rnt_batteryfan.data(k, :, :));
  
  [psd_fan_k, ~] = pwelch(dat_chan_k, Npts, [], Npts, 1/Ts, 'onesided', 'Power');
  psd_fan(:,k) = mean(psd_fan_k(start_idx:end,:), 2);
  
  [psd_nofan_k, ~] = pwelch(rnt_int.data(:, :, k)', Npts, [], Npts, 1/Ts, 'onesided', 'Power');
  psd_nofan(:,k) = mean(psd_nofan_k(start_idx:end,:), 2);
  
  [psd_batteryfan_k, f] = pwelch(rnt_batteryfan.data(:, :, k)', Npts, [], Npts, 1/Ts, 'onesided', 'Power');
  psd_batteryfan(:,k) = mean(psd_batteryfan_k(start_idx:end,:), 2);
  
  
end
freqs = f(start_idx:end);


%%

figure(1), clf;
for k=1:3
subplot(3,1,k)

h_fan = semilogx(freqs, 10*log10(psd_fan(:,k)) );

hold on
h_nofan = semilogx(freqs, 10*log10(psd_nofan(:,k)) );
h_batteryfan = semilogx(freqs, 10*log10(psd_batteryfan(:, k)));
grid on
end

h_fan.DisplayName = 'With standard npoint fan';
h_nofan.DisplayName = 'Fan unplugged';
h_batteryfan.DisplayName = 'Fan with 9v battery';

subplot(3,1,1)
legend([h_fan, h_nofan, h_batteryfan])
