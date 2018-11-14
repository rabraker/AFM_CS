% Implement augmented state space integral control with shift register to
% 'estimate' the delay states.


% clear all
clc
% Options
figbase  = 50;
verbose = 0;
controlParamName = 'LinControls01.csv';
refTrajName      = 'ref_traj_track.csv';
outputDataName = 'exp01outputBOTH.csv';
% Build data paths

addpath(fullfile(getCsRoot(), 'matlab-code', 'functions'));
addpath(fullfile(getCsRoot(), 'matlab-code', 'functions', 'state_space_x'));

full_mod_path = fullfile(PATHS.sysid, 'xy-axis_sines_info_intsamps_quickFourierCoef_11-11-2018xydrive-01.mat');
mf_full = load(full_mod_path);
% ---- Paths for shuffling data to labview and back. ------
%labview reads data here
controlDataPath = fullfile(PATHS.step_exp, controlParamName);
% labview saves experimental results/data here
data_out_path    = fullfile(PATHS.step_exp, outputDataName);
% labview reads desired trajectory here
traj_path     = fullfile(PATHS.step_exp, refTrajName);
% location of the vi which runs the experiment.


TOL = 14/512; % max volts by pixels
% TOL = .01;
tol_mode = 'abs';
% which simulations to run

do_sim_hyst = false;
do_inv_hyst = false;
do_drift = false;
do_invdrift = false;

plotstate = false;
plotpoles = false;
plotdriftpoles = false;
plotdriftbode = false;
saveon = true;

md = 1;
gam_rob = 46.4;
gam_s = sort([gam_rob]);
exp_id_str = 'const-sig';

% ------- Load Plants -----
[plants, frf_data] = CanonPlants.plants_ns14(9);
Ts  = plants.SYS.Ts;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                  Design reference "trajectory"                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get a ref trajectory to track.
N  = 800;
r1 =.4;

% step_ref = StepRef([0, r1], N);
% yref = step_ref.yref;
t = (0:floor(1.5/Ts))'*Ts;
N = length(t);
n1 = floor(N/3);
n2 = floor(N*(2/3));
u1 = zeros(n1,1);
u2 = linspace(0, 1, n2-n1)';
u3 = linspace(u2(end), 0, N-n2)';

yref = timeseries([u1; u2; u3], t);
step_ref = yref;


yref.Data = yref.Data*1;
step_descr = 'single_step';

dist_traj = yref;
dist_trajy = yref;
sig = 0.006;
dist_traj.Data = randn(size(yref.Data,1), 1)*sig;
dist_trajy.Data = randn(size(yref.Data,1), 1)*sig;
thenoise = dist_traj;
thenoise.Data = thenoise.Data*0;

% F_yudu = figure(60); clf
% subplot(3,1,1)
% hold on, grid on;
% step_ref.plot(F_yudu, '-k', 'LineWidth', 0.5);
% F_y = figure(61); clf
% hold on, grid on
% if max(abs(yref.Data)) > 0
%   step_ref.plot(F_y);
%   step_ref.plot_settle_boundary(F_y, TOL, tol_mode);
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Design control/estimator gains. This gains are what we actually         %
% use in both the simulation and experimentant system.                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3). LQR generation gain.
% -------------------------------------------------------------------------
% Adjust the du_max to account for the gain of gdrift_inv.
du_max_orig = StageParams.du_max;
du_max = du_max_orig/norm(plants.gdrift_inv, Inf);

cmplx_rad = 0.9;
% [Q1, R0, S1, P_x] = build_control_constsigma(plants.sys_recyc, cmplx_rad);
can_cntrl = CanonCntrlParamsChoozeZeta();
[Q1, R0, S1, P_x] = build_control_choosezeta(plants.sys_recyc, can_cntrl);

% ------------------------- Observer Gain ---------------------------------
can_obs_params = CanonObsParams_01();
can_obs_params.beta = 50;
[sys_obsDist, L_dist] = build_obs(plants.SYS, can_obs_params);

% 2). Design FeedForward gains.
[Nx, Nu] = SSTools.getNxNu(plants.sys_recyc);

% -------------------------------------------------------------------
% -------------------- Setup Fixed Point stuff -----------------------------
A_obs_cl = sys_obsDist.a - L_dist*sys_obsDist.c;



gam_iter = gam_s(1);
fprintf('===========================================================\n');

K_lqr = dlqr(plants.sys_recyc.a, plants.sys_recyc.b, Q1, R0+gam_iter, S1);
Qp = dare(plants.sys_recyc.a, plants.sys_recyc.b, Q1, R0+gam_iter, S1);

if 0
  verbose = 0;
  analyze_margins(plants, sys_obsDist, K_lqr, L_dist, verbose);
end


% --------------------  Fixed Linear stuff -----------------------------
%%
clc
uss_0 = 0;
PLANT = absorbDelay(mf_full.modelFit.G_all_axes);
sim_obj.sys_obs = sys_obsDist;
sim_obj.Nx = Nx;
sim_obj.du_max = du_max;
sim_obj.L = L_dist;
sim_obj.x0_obs = zeros(size(sys_obsDist.a, 1), 1);

x0 = zeros(size(PLANT.a, 1),1); %SSTools.getNxNu(PLANT(1,1))*0;

[ ndist, Ns_obs] = size(sim_obj.sys_obs.c);
Ident_obs = eye(Ns_obs);
 % the last row
 C_ydist = Ident_obs(end-ndist+1:end, :);
 % all rows but the last 1
 Ident_obs = Ident_obs(1:end-ndist, :);
 ref_traj = yref;
 trun = ref_traj.Time(end); 
sim('AFMss_fp_obshas_uk_threeaxis') 
figure(1)
plot(Y)
%%
dat=csvread('/media/labserver/afm-cs/imaging/raster/5microns/raster_scan_512pix_5mic_01Hz_out_11-12-2018mica_-01.csv');
%%
parent_dat = csvread('/media/labserver/afm-cs/imaging/raster/5microns/raster_scan_512pix_5mic_01Hz.csv');
xyref = reshape(parent_dat', 2, [])';
xref = xyref(:,1);

% figure(1); plot(datmat(:,1));
% ax1 = gca
% figure(2); plot(datmat(:,3));
% ax3 = gca;
% figure(3); plot(datmat(:,4));
% ax4 = gca;
% linkaxes([ax1, ax3, ax4])
Ts = 40e-6;
samps_per_period = size(parent_dat,1)/2 % twice as many in here for x & y.
samps_per_line = samps_per_period/2

nperiods = 512;
pix = nperiods;
dat = dat([1:nperiods*samps_per_period], :);
figure(2)
subplot(2,1,1)
ax1 = gca();
N1 = 2;
N2 = N1+6;
idx = N1*samps_per_line+1:N2*samps_per_line;
t = idx*Ts;
plot(t,dat(idx,3))
subplot(2,1,2)
plot(t, dat(idx, 4))
ax2 = gca();

linkaxes([ax1, ax2], 'x')
%%
clc
D = zpk([0], [1], -0.02, Ts);
figure(10); clf;
subplot(2,1,1)
subplot(2,1,2)
UXX = 0;

ZZ = 0;
%  idx_trace = get_trace_indeces(nperiods, samps_per_period);
%  idx_retrace = get_retrace_indeces(nperiods, samps_per_period);
 start_idx = 500;
for k=0:512-1
  idx = k*(samps_per_line*2)+1:(k+1)*(samps_per_line*2);
  ux = dat(idx(start_idx:end-99), 4);
  ux = ux - mean(ux);
  [UX, freqs] = power_spectrum(ux, Ts);
  UXX = UXX +UX;

  z = dat(idx(start_idx:end-99), 3);
  z = z-mean(z);
  [Z] = power_spectrum(z, Ts);
  ZZ = ZZ + Z;
  
  figure(10);
  subplot(2,1,1), cla
  semilogx(freqs, log10(ZZ/k))
  grid on
  
  subplot(2,1,2), cla
  semilogx(freqs, log10(UXX/k))
  hold on
  semilogx(freqs, log10( (ZZ/k) .* abs( squeeze(freqresp(D, freqs*2*pi)))), '--r')
  grid on
  drawnow();
  
  figure(11)
  subplot(2,1,1)
  plot(z)
  subplot(2,1,2)
  plot(ux)
  
%   keyboard

end
%%
figure(10)
subplot(2,1,1)
title('ux')
subplot(212)
title('ze')

%%
% sims_fpl = SimAFM(absorbDelay(mf_full.modelFit.G_all_axes(1,1)), K_lqr, Nx, sys_obsDist, L_dist, du_max,...
%   false, 'thenoise', thenoise);
% 
% 
% [, U_full_fxpl, U_nom_fxpl, dU_fxpl, Xhat_fxpl] = sims_fpl.sim(yref, dist_traj);
% name = sprintf('FXP lin Sim. (%s)', exp_id_str);
% 
% fxpl_Opts = stepExpDuOpts('pstyle', '-r', 'TOL', TOL, 'step_ref', step_ref,...
%   'controller', K_lqr, 'name',  name);
% sim_exp_fxpl = stepExpDu(y_fxpl, U_full_fxpl, dU_fxpl, fxpl_Opts);
% Ts_vec_fxpl = sim_exp_fxpl.settle_time(TOL, tol_mode, 0);
% 
% h2 = plot(sim_exp_fxpl, F_yudu, 'umode', 'both');
% legend([h2(1)])
% figure(F_y)
% h22 = sim_exp_fxpl.ploty(F_y);
% legend([h22]);



%%
function handle_vec =  make_legend_vec(varargin)
  
%   handle_vec = gobjects();
handle_vec = [];
  for k = 1:length(varargin)
    hand = varargin{k};
    if isvalid(hand)
      handle_vec = [handle_vec, hand];
    end
  end
  
end

function analyze_margins(plants, sys_obsDist, K_lqr, L_dist, verbose)
  [Sens, Hyd, Hyr, Hyeta, Loop] = ss_loops_delta_dist(plants.SYS, plants.sys_recyc,...
    sys_obsDist, K_lqr, L_dist);
  
  F_clbode = figure(25);clf; hold on, grid on
  Hbode_sens = bodeplot(Sens);
  setoptions(Hbode_sens, 'FreqUnits', 'Hz')
  legend('S: linear')
  grid on, hold on;
  
  [Gm_lin, Pm_lin] = margin(Loop);
  
  fprintf('-------- MARGINS ------------------\n')
  fprintf('Linear: GM = %.2f [], PM = %.2f [deg]\n', Gm_lin, Pm_lin)
  
  if verbose >=2  
    figure(101)
    rlocus(Loop);
    title('Klin')
  end
  
  if verbose >=1
    figure(104)
    nyquist(Loop)
    title('Klin')
    xlim([-2, 0.1])
    ylim([-1, 1])
  end
  
  
end


