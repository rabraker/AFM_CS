% Implement augmented state space integral control with shift register to
% 'estimate' the delay states.


clear all
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


% exp_id_str = 'const-sig';
exp_id_str = 'choose-zeta';

if strcmp(exp_id_str, 'const-sig')
  gam_rob = 46.4;
elseif strcmp(exp_id_str, 'choose-zeta')
  gam_rob = 25;
else
  error('unrecognized exp_id_str: %s', exp_id_str)
end
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
r1 =1.37;

step_ref = StepRef([r1], N);
yref = step_ref.yref;
yref.Data = yref.Data*1;
step_descr = 'single_step';

dist_traj = yref;
dist_traj.Data = dist_traj.Data*0;
thenoise = dist_traj;


F_yudu = figure(60); clf
subplot(3,1,1)
hold on, grid on;
step_ref.plot(F_yudu, '-k', 'LineWidth', 0.5);

F_y = figure(61); clf
hold on, grid on
if max(abs(yref.Data)) > 0
  step_ref.plot(F_y);
  step_ref.plot_settle_boundary(F_y, TOL, tol_mode);
end

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

if strcmp(exp_id_str, 'const-sig')
  cmplx_rad = 0.9;
  [Q1, R0, S1, P_x] = build_control_constsigma(plants.sys_recyc, cmplx_rad);
elseif strcmp(exp_id_str, 'choose-zeta')
  can_cntrl = CanonCntrlParamsChoozeZeta();
  [Q1, R0, S1, P_x] = build_control_choosezeta(plants.sys_recyc, can_cntrl);
end

% ------------------------- Observer Gain ---------------------------------
can_obs_params = CanonObsParams_01();
can_obs_params.beta = 50;
[sys_obsDist, L_dist] = build_obs(plants.SYS, can_obs_params);

% 2). Design FeedForward gains.
[Nx, Nu] = SSTools.getNxNu(plants.sys_recyc);


fprintf('===========================================================\n');

K_lqr = dlqr(plants.sys_recyc.a, plants.sys_recyc.b, Q1, R0+gam_rob, S1);
Qp = dare(plants.sys_recyc.a, plants.sys_recyc.b, Q1, R0+gam_rob, S1);

if 1
  verbose = 0;
  analyze_margins(plants, sys_obsDist, K_lqr, L_dist, verbose);
end

%%
% -------------------------------------------------------------------
% -------------------- Setup Fixed Point stuff -----------------------------
A_obs_cl = sys_obsDist.a - L_dist*sys_obsDist.c;

nw = 32;
nf = 26;

du_max_fxp = fi(du_max, 1, 32, 26);
Nx_fxp = fi(Nx, 1, 32, 30);
L_fxp = fi(L_dist, 1, 32, 30);

sys_obs_fxp.a = fi(sys_obsDist.a -L_dist*sys_obsDist.c, 1, nw, nw-7);
sys_obs_fxp.b = fi(sys_obsDist.b, 1, nw, 29);
sys_obs_fxp.c = fi(sys_obsDist.c, 1, nw, 28);

K_fxp = fi(K_lqr, 1, nw,32-10);
sims_fxpl = SimAFM(plants.PLANT, K_fxp, Nx_fxp, sys_obs_fxp, L_fxp, du_max_fxp,...
  true, 'nw', nw, 'nf', nf, 'thenoise', thenoise);
if 1
  sims_fxpl.r = plants.hyst_sat.r;
  sims_fxpl.w = plants.hyst_sat.w;
  sims_fxpl.rp = fi(plants.hyst_sat.rp, 1, 16, 11);
  sims_fxpl.wp = fi(plants.hyst_sat.wp, 1, 16, 11);
  sims_fxpl.d = plants.hyst_sat.d;
  sims_fxpl.ws = plants.hyst_sat.ws;
  sims_fxpl.dp = fi(plants.hyst_sat.dp, 1, 16, 11);
  sims_fxpl.wsp = fi(plants.hyst_sat.wsp, 1, 16, 11);
end
if 1
  sims_fxpl.gdrift_inv = plants.gdrift_inv;
  sims_fxpl.gdrift = plants.gdrift;
end

[y_fxpl, U_full_fxpl, U_nom_fxpl, dU_fxpl, Xhat_fxpl] = sims_fxpl.sim(yref, dist_traj);
name = sprintf('FXP lin Sim. (%s)', exp_id_str);

fxpl_Opts = stepExpDuOpts('pstyle', '-r', 'TOL', TOL, 'step_ref', step_ref,...
  'controller', K_lqr, 'name',  name);
sim_exp_fxpl = stepExpDu(y_fxpl, U_full_fxpl, dU_fxpl, fxpl_Opts);
Ts_vec_fxpl = sim_exp_fxpl.settle_time(TOL, tol_mode, 0);

h2 = plot(sim_exp_fxpl, F_yudu, 'umode', 'both');
legend([h2(1)])
figure(F_y)
h22 = sim_exp_fxpl.ploty(F_y);
legend([h22]);

%%
fprintf('===========================================================\n');
fprintf('Writing control data...\n');
fprintf('===========================================================\n');
clc

sims_fxpl.sys_obs_fp = sys_obsDist;
sims_fxpl.sys_obs_fp.a = sys_obsDist.a - L_dist*sys_obsDist.c;

sims_fxpl.write_control_data(controlDataPath, yref, traj_path)
control_path = fullfile(PATHS.step_exp, sprintf('LinControls-%s.json', exp_id_str))
sims_fxpl.write_control_data_json(control_path)

%%


% column spec: [gamma, lin_exp, mpc_exp]

%--------------------------------------------------------------------------
% --------------------------- LINEAR Experiment ---------------------------

fprintf('===========================================================\n');
fprintf('Starting Experiments for gamma = %\n', gam_rob);
fprintf('===========================================================\n');
% Build the u-reset.
if 1
  dry_run = false;
  reset_piezo('t1', 15, 't_final', 25, 'umax', 9, 'k1', 0.55,...
    'verbose', false, 'dry_run', dry_run)
  fprintf('...finished running piezo-reset.\n')
end

SettleTicks = 20000;
Iters = length(yref.Data)-1;

% create and pack data. Then save it.
[num, den] = tfdata(plants.gdrift_inv);
num = num{1};
den = den{1};

umax = 11;
ymax = max(yref.Data)*1.3;
%%
clear e;
clear vi;
% -----------------------RUN THE Experiment--------------------------------
vipath =['C:\Users\arnold\Documents\matlab\afm_mpc_journal\',...
  'labview\fixed-point-host\play_FXP_AFMss_LinearDistEst_singleAxis.vi'];

[e, vi] = setup_VI(vipath, false, 'SettleTicks', SettleTicks, 'Iters', Iters,...
  'num', num, 'den', den, 'TF Order', 1*(length(den)-1),...
  'r_s', plants.hyst_sat.rp, 'w_s', plants.hyst_sat.wp, 'N_hyst', 1*length(plants.hyst_sat.rp),...
  'sat_ds', plants.hyst_sat.dp, 'sat_ws', plants.hyst_sat.wsp, 'N_sat', 1*length(plants.hyst_sat.dp),...
  'du_max', du_max,'dry_run', false,...
  'read_file', true, 'umax', umax, 'ymax', ymax, 'outputDataPath', data_out_path,...
  'traj_path', traj_path, 'control_data_path', fxplin_dat_path);

vi.Run


% Now, read in data, and save to structure, and plot.
[y_exp, ~, du_exp, ufull_exp, Ipow_exp, xhat_exp, yy] =  load_exp_data(data_out_path, sys_obsDist);

expOpts = stepExpDuOpts('TOL', TOL, 'step_ref', step_ref, 'controller', K_lqr,...
  'pstyle', '-b', 'name', sprintf('AFM Stage (Linear) (%s)', exp_id_str));

afm_exp_lin = stepExpDu(y_exp, ufull_exp, du_exp, expOpts);
afm_exp_lin.Ipow = Ipow_exp;

Ts_vec_afm_lin = afm_exp_lin.settle_time(TOL, tol_mode, 0);
fprintf('Total AFM lin FXP settle-time = %.3f [ms]\n', sum(Ts_vec_afm_lin)*1000);

H_linexp = plot(afm_exp_lin, F_yudu, 'umode', 'both');
subplot(3,1,1)

legend(make_legend_vec(h1(1), h2(1), h3(1), H_linexp(1), H_mpcexp(1) )); 

plot(y_exp.Time, yy, ':k')

H_linexp2 = afm_exp_lin.ploty(F_y);

legend(make_legend_vec(h12, h22, h32, H_linexp2 ));



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


function [y, u, du, ufull, Ipow, xhat, yhat] =  load_exp_data(data_out_path, sys_obsDist)
  Ts = sys_obsDist.Ts;
  AFMdata = csvread(data_out_path);

  t_exp = (0:size(AFMdata,1)-1)'*Ts;
  y = timeseries(AFMdata(:,1), t_exp);
  u = timeseries(AFMdata(:, 2), t_exp);
  du = timeseries(AFMdata(:,3), t_exp);
  ufull = timeseries(AFMdata(:,4), t_exp);
  
  Ipow = timeseries(AFMdata(:,5), t_exp);
  xhat = timeseries(AFMdata(:,6:end), t_exp);
  yhat = xhat.Data*sys_obsDist.c';
  
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


