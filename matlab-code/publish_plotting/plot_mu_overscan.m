clear
rmpath ~/gradschool/thesis/afm_mpc_journal/functions/canon/
rmpath ~/gradschool/thesis/afm_mpc_journal/functions
addpath functions
xdirControl = get_xdir_standard_control();
Hyr = xdirControl.Hyr;
%%
Ts = Hyr.Ts;

width =5;  % microns
raster_freq = 2;

% raster_period = 1/raster_freq;
microns_per_second = 2*width * raster_freq;
x_rate = microns_per_second * AFM.mic2volt_xy() * AFM.Ts  % Ts = seconds per sample


mu_micron = 0.5
mu_volts = mu_micron * AFM.mic2volt_xy();

mu_Nsamples = ceil(mu_volts / x_rate);

%%

% TF from reference to error.
Her = minreal(1 - Hyr); %%feedback(1, D*G);
Int_z = zpk([], [1], 1, Ts);
% Compute the steady state error to a ramp
ess = dcgain(minreal(Her*Int_z))*x_rate;

% sometimes, the integrator gets perturbed to a really slow pole, so dc-gain
% is infinitate. So evaluate low-freq gain instead.
if isinf(ess)
  ess = abs(freqresp(Her*Int_z, 1))*x_rate;
end


N_extra = ceil(ess/x_rate)


% figure(100); clf
fig0 = mkfig(100, 5, 4); clf
ha = tight_subplot(1, 1, 0, [.1, .02], [.32, .05])

% 1. Niave, no pre-scan, no post scan
x_N = (mu_Nsamples)*x_rate;
r_vec = linspace(0, x_N, mu_Nsamples);
t_vec = [0:1:mu_Nsamples-1]'*Ts;
y_vec = lsim(Hyr, r_vec, t_vec);

plot(t_vec, r_vec, 'b', 'LineWidth', 1.5)
hold(ha, 'on')
plot(t_vec, y_vec, 'r', 'LineWidth', 1.5)
remove_ticks(ha)
xlim(ha, [0, 1.25*t_vec(end)])
ylim(ha, [0, 1.25*r_vec(end)])

plot(ha, [t_vec(end), t_vec(end)], [0, r_vec(end)], '--k', 'LineWidth', 1.1)
plot(ha, [0, t_vec(end)], [r_vec(end), r_vec(end)], '--k', 'LineWidth', 1.1)


ylabel('y(k)', 'FontSize', 14)
xlabel('time', 'FontSize', 14)
%
N = mu_Nsamples + N_extra;
x_N = (N)*x_rate;
r_vec_post = linspace(0, x_N, N);
t_vec_post = [0:1:N-1]'*Ts;
y_vec_post = lsim(Hyr, r_vec_post, t_vec_post);

plot(t_vec_post, r_vec_post, 'b', 'LineWidth', 1.5)
hold(ha, 'on')
plot(t_vec_post, y_vec_post, 'r', 'LineWidth', 1.5)

plot(ha, [t_vec_post(end), t_vec_post(end)], [0, r_vec_post(end)], '--k', 'LineWidth', 1.1)
plot(ha, [0, t_vec_post(end)], [r_vec(end), r_vec(end)], '--k')
plot(ha, [0, t_vec_post(end)], [r_vec_post(end), r_vec_post(end)], '--k', 'LineWidth', 1.1)


text(0, r_vec(end), '$y(t_{\textrm{scan}})$', 'HorizontalAlignment', 'right', 'FontSize', 12);
t1 = text(t_vec(end)-t_vec(end)/2.5, r_vec(end)-r_vec(end)/10,...
    '$e_{ss}=H_{er}(1)v$', 'HorizontalAlignment', 'right', 'FontSize', 12);

a1 = annotation('arrow', [0.6417 0.82], [0.7328 0.7240])


text(mu_Nsamples*Ts, 0.002, '$t_{scan}$',...
    'HorizontalAlignment', 'left', 'VerticalAlignment',...
    'bottom','FontSize', 18, 'rot', 90);

s=sprintf('$t_{scan}+ t_{post}$');

text(t_vec_post(end), +0.001, s,...
    'HorizontalAlignment', 'left', 'VerticalAlignment',...
    'bottom','FontSize', 18, 'rot', 90);
s2 = sprintf('$y(t_{\\textrm{scan}} + t_{\\textrm{post}})$');
text(0, r_vec_post(end), s2, 'HorizontalAlignment', 'right', 'FontSize', 12);
%
save_fig(fig0, fullfile(PATHS.defense_fig(), 'overscan'), true)
%%
% When we do the prescan, we reach quasi steady state by the end, typically.
% This means that we dont necessarily need the all the overscan samples,
% because what we are really interested in is getting a scan that spans at
% mu_len microns.
clc
N_prescan = N_extra;
N = mu_Nsamples + N_prescan;
x_N = (N)*x_rate;
r_vec = linspace(0, x_N, N);
t_vec = [0:1:N-1]'*Ts;
y_vec = lsim(Hyr, r_vec, t_vec);

% y_vec(end) - y_vec(N_prescan)
% mu_volts
% ess



% figure(100); clf
fig1 = mkfig(101, 5, 4); clf
% ha = tight_subplot(1, 1, 0, [.1, .02], [.26, .02])
ha = tight_subplot(1, 1, 0, [.1, .02], [.32, .05])
hold on




plot(t_vec, y_vec, 'r', 'LineWidth', 1.5);
plot(t_vec, r_vec, 'b', 'LineWidth', 1.5);
% plot(t_vec, y_vec, '--')
remove_ticks(ha)
xlim(ha, [0, 1.1*t_vec(end)])
ylim(ha, [0, 1.1*r_vec(end)])

% grid
xlabel('time ', 'FontSize', 14)
ylabel('y(k)', 'FontSize', 14)


N_pre_fake = 2*N_prescan;

plot([N_pre_fake, N_pre_fake]*Ts, [0, r_vec(N_pre_fake)], '--k', 'LineWidth', 1.1)

plot([N, N]*Ts, [0, r_vec(N_prescan+mu_Nsamples)], '--k', 'LineWidth', 1.1)
plot([0, N*Ts], [1, 1]*y_vec(end), '--k')

% X0 is the start of the mu-path scan, which happens after the last prescan
% sample.

x0 = y_vec(N_prescan);
x0_fake = y_vec(N_pre_fake+1);
% X1 is the mu_path
mu_len_volts = mu_Nsamples*x_rate;
x1 = x0+mu_len_volts;

% x1_idx = find(y_vec >= x1, 1, 'first')
% plot([0, x1_idx]*Ts, [x1, x1], '--k', 'LineWidth', 1.5)

plot([0, N_pre_fake+1]*Ts, [x0_fake, x0_fake], '--k', 'LineWidth', 1.5)


% set(gca(), 'YTick', [x0, x1])
set(gca(), 'YTickLabel', [], 'XTickLabel', [])

txt_opts_xl = {'Units', 'Data', 'HorizontalAlignment', 'center', 'FontSize', 12};
txt_opts_yl = {'Units', 'Data', 'HorizontalAlignment', 'right', 'FontSize', 12};

text(-20*Ts, x1, '$y(t_{prescan}+t_{scan})$', txt_opts_yl{:})
text(-20*Ts, x0_fake, '$y(t_{prescan})$', txt_opts_yl{:})
% text(-20*Ts, 0, '$0$', txt_opts_yl{:})

yly=-.005;
text(N_pre_fake*Ts, .002, '$t_{prescan}$', txt_opts_xl{:},...
    'HorizontalAlignment', 'left', 'VerticalAlignment',...
    'top','FontSize', 18, 'rot', 90);

s=sprintf('$t_{prescan} + t_{scan}$')
text(N*Ts, .002, s,...
    'HorizontalAlignment', 'left', 'VerticalAlignment',...
    'bottom','FontSize', 18, 'rot', 90);

a2 = annotation('DoubleArrow', [0.3646 0.3646], [0.2120 0.8])
st3 = text(93*Ts, 0.0335, '$\mu$-path length', 'FontSize', 18, 'rot', 90);


save_fig(fig1, fullfile(PATHS.defense_fig(), 'pre_scan'), true)


