function xdirControl = get_xdir_standard_control()
addpath(fullfile(getCsRoot(), 'matlab-code', 'functions'));
addpath(fullfile(getCsRoot(), 'matlab-code', 'functions', 'state_space_x'));

md = 1;
gam_rob = 46.4;


% ------- Load Plants -----
[plants, frf_data] = CanonPlants.plants_ns14(9);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Design control/estimator gains. This gains are what we actually         %
% use in both the simulation and experimentant system.                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3). LQR generation gain.
% -------------------------------------------------------------------------

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
K_lqr = dlqr(plants.sys_recyc.a, plants.sys_recyc.b, Q1, R0+gam_rob, S1);

if 0
  verbose = 0;
  analyze_margins(plants, sys_obsDist, K_lqr, L_dist, verbose);
end

[Sens, Hyd, Hyr, Hyeta, Loop] = ss_loops_delta_dist(plants.SYS, plants.sys_recyc,...
  sys_obsDist, K_lqr, L_dist);
xdirControl = struct('Sens', Sens,...
  'Hyd', Hyd,...
  'Hyr', Hyr,...
  'Loop', Loop,...
  'plants', plants,...
  'sys_obsDist', sys_obsDist,...
  'Nx', Nx,...
  'K_lqr', K_lqr);
  


end