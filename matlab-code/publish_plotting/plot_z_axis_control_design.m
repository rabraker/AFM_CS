addpath functions
close all
root = fullfile(PATHS.sysid, 'z_axis_evolution', 'batch-1');
file = 'first_res_fit-3-17-2019-1.json';



Fig = mkfig(1, 5, 4); clf
[ha, pos] = tight_subplot(1, 1, [.02, .01 ], [.1, 0.03], [.1, .05]);

if k==1
  prefix='(A) ';
else
  prefix = '(B) ';
end
dat = loadjson(fullfile(root, file));

gz_k = dat.SOS_frf_data.resp_real + dat.SOS_frf_data.resp_imag*1i;
h(k) = frf_bode_mag(gz_k, dat.SOS_frf_data.omegas/2/pi, Fig, 'Hz');
h(k).DisplayName = [prefix name];


ss_fname = fullfile(root, 'ALL-axis_sines_info_intsamps_quickFourierCoef_3-17-2019-01.json');
ss_data = SweptSinesOnline(ss_fname);

G_frf = frd(ss_data.FC_s(:,2)./ss_data.FC_s(:,1), ss_data.freq_s, AFM.Ts, 'FrequencyUnit', 'Hz');


KI = -0.03;

D_inv = frd(tf(dat.Dinv_Num, dat.Dinv_den, AFM.Ts),  ss_data.freq_s,  'FrequencyUnit', 'Hz');

D_ki = zpk([], [1], KI, AFM.Ts);

DD = D_ki * D_inv;

D_ki_frf = frd(D_ki, ss_data.freq_s,  'FrequencyUnit', 'Hz');
DD_frf = frd(DD, ss_data.freq_s, 'FrequencyUnit', 'Hz');
Dinv_frf = frd(D_inv, ss_data.freq_s, 'FrequencyUnit', 'Hz');

Loop = DD_frf * G_frf;
Loop_noinv = D_ki_frf * G_frf;


Huz_d = feedback(D_ki_frf, D_inv*G_frf);
Hyr = feedback(D_ki_frf*D_inv*G_frf, 1);
Hyr_noinv = feedback(D_ki_frf,  G_frf);


F2 = mkfig(2, 7, 4); clf
[ha, pos] = tight_subplot(1, 1, [.02, .01 ], [.1, 0.03], [.1, .05]);

h1 = frf_bode_mag(Loop, ss_data.freq_s, ha, 'Hz', '-g');
h2 = frf_bode_mag(G_frf, ss_data.freq_s, ha, 'Hz');
h3 = frf_bode_mag(Hyr_noinv, ss_data.freq_s, ha, 'Hz', '-r');
h4 = frf_bode_mag(Huz_d, ss_data.freq_s, ha, 'Hz', '-k');
ylim(ha(1), [-60, 20]);

h1.DisplayName = 'Loop gain (with inverse)';
h2.DisplayName = '$G_{d_{fl},u_z}$';
h3.DisplayName = '$H_{u_{Z},s}$ (no inverse)';
h4.DisplayName = '$H_{u_{Z},s}$ (with inverse)';
leg = legend([h1, h2, h3, h4]);

set(leg, 'Location', 'southwest', 'FontSize', 12);

save_fig(F2, fullfile(PATHS.thesis_fig_final, 'z_control_design'));
%%
bandwidth(Huz_d)/2/pi
bandwidth(Hyr)/2/pi

%%
clc
addpath('system_id')
F3 = mkfig(3, 7, 4); clf
[ha, pos] = tight_subplot(1, 1, [.02, .01 ], [.1, 0.03], [.1, .05]);

Ts = AFM.Ts;
ss_opts = frf2ss_opts('Ts', Ts, 'r', 500, 's', 300);
Nd2 = 4;
ns2 = 7;
Kc = 1;
freqs = ss_data.freq_s(:);
f2ss = frf2ss(Huz_d.ResponseData*Kc, freqs*2*pi, Nd2, ss_opts); % 12
sys_z_cl = f2ss.realize(ns2); % 12
k_estmax = floor(length(freqs)/2)
LGopts = optimoptions(@lsqnonlin, 'Display', 'iter',...
    'FunctionTolerance', 1e-9, 'MaxIter', 5000,'MaxFunctionEvaluations', 5000,...
    'StepTolerance', 1e-9, 'Jacobian','on', 'CheckGradients', false);

sos_fos = SosFos(zpk(sys_z_cl), 'iodelay', sys_z_cl.InputDelay);
LG = LogCostZPK(Huz_d.ResponseData(1:k_estmax)*Kc, freqs(1:k_estmax)*2*pi, sos_fos);
LG.solve_lsq(2, LGopts)
[sys_z_log, p] = LG.sos_fos.realize();
sys_z_log.InputDelay = 1; %max(round(p, 0), 0);
fprintf('LG says delay = %.2f\n', p);
sys_z_cl = sys_z_log;


h4 = frfBode(Huz_d*Kc, ss_data.freq_s, F3, 'Hz', '-k');
frfBode(sys_z_cl, ss_data.freq_s(:), F3, 'Hz')
% h4 = frf_bode_mag(Hyr*Kc, ss_data.freq_s, ha, 'Hz', '-k');
% frf_bode_mag(sys_z_cl, ss_data.freq_s(:), ha, 'Hz')


width = 5;

A = 0.07;
Ts = AFM.Ts;
rast_freq = 6; %Hz
period = 1/rast_freq;
holes_per_line = 10;
holes_per_period = 2*holes_per_line
holes_per_sec =  holes_per_period / period;

t_half_pitch = 0.5 * (1/holes_per_sec);

N_half = floor(t_half_pitch / Ts)

one_hole = [zeros(N_half, 1); zeros(N_half,1)+A];
N = 10
N_hole = repmat(one_hole, N, 1);
% t_one_hole = (0:1:2*N_half-1)'*Ts;
t = (0:length(N_hole)-1)'*Ts;

% plot(t_one_hole, one_hole)

sys_ze = (1/D_ki)*sys_z_cl;
figure(6); clf
yy = lsim(sys_z_cl, N_hole, t);
lsim(sys_z_cl, N_hole, t);


[y, t] = lsim(sys_ze, N_hole, t);

%
figure(38); clf
hold on
hye = plot(t+22.0955, y-0.3)

damage_metric(y)
%%
y_=[]
for k=1:4:floor(length(y))
  if k+4 > length(y)
    break
  end
  y_ = [y_; mean(yy(k:k+4))];
  
end

figure, plot(y_)



function damage = damage_metric(y)
  % Compute a damage metric based on the deflection signals positivity.
  % This is computed as the power of the positive values of the negative
  % error signal. The motivation is that, for a given setpoint, we do not
  % care, from a damage perspective, if the error dips below the setpoint
  % (though that will affect image quality), because this corresponds to the
  % tip parachiting off a ledge. Rather from a damage perspective, what we
  % care about is events where (ze - ref) signal becomes positive.
  
  ref = 0; %y.meta_exp.z_axis_params.setpoint_scan;
  % rather than subtracting mean, subtract the reference value.
  err = y - ref;  % shift to zero.
  err_pos = err(err>0);
  damage = sum(err_pos.^2)/length(err_pos)/AFM.Ts;
end


