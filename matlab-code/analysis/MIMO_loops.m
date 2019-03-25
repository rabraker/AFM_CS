
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

Gxy = mf_full.modelFit.G_all_axes(1:3,:);

% ------- Load Plants -----
[plants, frf_data] = CanonPlants.plants_ns14(9);
Ts  = plants.SYS.Ts;
Dy = zpk([0], [1], 0.01, AFM.Ts);

du_max_orig = StageParams.du_max;
du_max = du_max_orig/norm(plants.gdrift_inv, Inf);

cmplx_rad = 0.9;
[Q2, R2, S2] = build_control_constsigma(plants.sys_recyc, cmplx_rad);
can_cntrl = CanonCntrlParamsChoozeZeta();
[Q1, R0, S1] = build_control_choosezeta(plants.sys_recyc, can_cntrl);

% ------------------------- Observer Gain ---------------------------------
can_obs_params = CanonObsParams_01();
can_obs_params.beta = 50;
[sys_obsDist, L_dist] = build_obs(plants.SYS, can_obs_params);

% 2). Design FeedForward gains.
[Nx, Nu] = SSTools.getNxNu(plants.sys_recyc);


gam_rob1 = 25;
K_lqr1 = dlqr(plants.sys_recyc.a, plants.sys_recyc.b, Q1, R0+gam_rob1, S1);

gam_rob2 = 100;
K_lqr2 = dlqr(plants.sys_recyc.a, plants.sys_recyc.b, Q2, R2+gam_rob2, S2);


% Qp = dare(plants.sys_recyc.a, plants.sys_recyc.b, Q1, R0+gam_rob, S1);

[~, ~, Hyr] = ss_loops_delta_dist(plants.SYS, plants.sys_recyc, sys_obsDist, K_lqr1, L_dist);


[HH1] = close_two_axis(plants.SYS, plants.sys_recyc, sys_obsDist, K_lqr1, L_dist, Gxy, Dy, 0);
HH1 = HH1(:, 1:2);
[HH2] = close_two_axis(plants.SYS, plants.sys_recyc, sys_obsDist, K_lqr2, L_dist, Gxy, Dy, 0);
HH2 = HH2(:, 1:2);



F11 = mkfig(10, 7, 9, true); clf
[ha, pos] = tight_subplot(3, 2, .01, .04, [.06, .02]);
ha = reshape(ha', 2, [])';
freq_s = logspace(log10(1), log10(12500), 350)';

mimo_bode_mag(HH1, freq_s, ha, '-k');
mimo_bode_mag(HH2, freq_s, ha, '--r');




function mimo_bode_mag(HH, freq_s, ha, varargin)
  ylm = [-95, 30];
  coords = {'X', 'Y', 'Z'};
  for ny = 1:3
    for nu=1:2
      frf_bode_mag(HH(ny,nu), freq_s, ha(ny,nu), 'Hz', varargin{:})
      set(ha(ny,nu), 'YLim', ylm);
      if nu >1
        set(ha(ny,nu), 'YLabel', [], 'YTickLabel', []);
      else
        ylabel(ha(ny,nu), sprintf('Output %s (Mag [dB] )', coords{ny}));
      end
      if ny ==1
        title(ha(ny,nu), sprintf('Input %s', coords{nu}));
      end
      if ny <3
        set(ha(ny,nu),  'XLabel', [], 'XTickLabel', []);
      else
        set(ha(ny,nu), 'XTick', [10,100,1000,10000]);
      end
    end
  end
  
end

function [HH] = close_two_axis(sys, sys_recyc, sys_obs, KxKu, LxLd, Gxy, Dy, Dz)
  Ts = sys.Ts;
  Kx = KxKu(1:end-1);
  Ku = KxKu(end);
  Nbar = SSTools.getNbar(sys_recyc, KxKu);

  Ld = LxLd(end);
  Lx = LxLd(1:end-1);
  Cd = sys_obs.c(end);

  A_tilde = sys.a - sys.b*Kx - Lx*sys.c;

  AA_ = [A_tilde, sys.b-sys.b*Ku, -(sys.b*Nbar+Lx*Cd);
       -Kx,     1-Ku,             -Nbar;
       -Ld*sys.c,    0,              1-Ld*Cd];

  BB_ = [sys.b*Nbar; Nbar; 0];
  LL_ = [Lx; 0; Ld];

  K_c = [KxKu, Nbar];
  K_c(end-1) = K_c(end-1)-1;
  D1 = ss(AA_, LL_, K_c, 0, Ts);
  D2 = ss(AA_, BB_, -K_c, Nbar, Ts);

  % % Estimated state feedback loop gain:
  Loop = sys*D1;
  Sens = 1/(1+Loop);

  Hyr = ((sys*D2*Sens));
  Hyd = minreal( sys*Sens);
  
 
  DD2 = [D2, 0, 0;
    0,  Dy, 0];
%     0, 0, Dz];
  DD1 = [D1, 0, 0;
         0, Dy, 0];
%     0, 0, Dz];
  I = eye(3);
  
  HH = Gxy*DD2/(I + Gxy*DD1);
  
  HH = minreal(HH);

end
