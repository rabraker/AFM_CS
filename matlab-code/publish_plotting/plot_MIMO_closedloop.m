
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

full_mod_path = fullfile(PATHS.sysid, 'xy-axis_sines_info_intsamps_quickFourierCoef_11-11-2018xydrive-01.mat');
mf_full = load(full_mod_path);

%%
Gxyz = mf_full.modelFit.frf.all_axes; 

% ------- Load Plants -----
[plants, frf_data] = CanonPlants.plants_ns14(9);
Ts  = plants.SYS.Ts;
Dy = zpk([0], [1], 0.01, AFM.Ts);

du_max_orig = StageParams.du_max;
du_max = du_max_orig/norm(plants.gdrift_inv, Inf);


can_cntrl = CanonCntrlParamsChoozeZeta();
[Q1, R0, S1] = build_control_choosezeta(plants.sys_recyc, can_cntrl);

cmplx_rad = 0.9;
[Q2, R2, S2] = build_control_constsigma(plants.sys_recyc, cmplx_rad);

% ------------------------- Observer Gain ---------------------------------
can_obs_params = CanonObsParams_01();
can_obs_params.beta = 50;
[sys_obsDist, L_dist] = build_obs(plants.SYS, can_obs_params);

gam_rob1 = 25;
K_lqr1 = dlqr(plants.sys_recyc.a, plants.sys_recyc.b, Q1, R0+gam_rob1, S1);

gam_rob2 = 46.4;
K_lqr2 = dlqr(plants.sys_recyc.a, plants.sys_recyc.b, Q2, R2+gam_rob2, S2);


% [~, ~, Hyr] = ss_loops_delta_dist(plants.SYS, plants.sys_recyc, sys_obsDist, K_lqr1, L_dist);

Gxyz_frd = get_H();
%%
% Gxyz_frd = [plants.gdrift_inv, 0, 0; 0, 1, 0; 0 0 1]* Gxyz_frd;


Gzz = Gxyz_frd(3,3);
[DDz, Dki, Dinv] = get_gz_dz(Gzz);


% Gz_interp = interp1(Gz_frf.Frequency, Gz_frf.ResponseData(:), mf_full.modelFit.frf.freqs_Hz);
% Gxyz(3,3,:) = Gz_interp;
% Gxyz(2,3,:) = 0;
% Gxyz(1,3,:) = 0;

% Gxyz_frd = frd(Gxyz, mf_full.modelFit.frf.freqs_Hz, 'FrequencyUnit', 'Hz', 'Ts', AFM.Ts);


[HH1] = close_three_axis(plants.SYS, plants.sys_recyc, sys_obsDist, K_lqr1, L_dist, Gxyz_frd, Dy, Dki, Dinv);

[HH2] = close_three_axis(plants.SYS, plants.sys_recyc, sys_obsDist, K_lqr2, L_dist, Gxyz_frd, Dy, Dki, Dinv);



 
F11 = mkfig(10, 7, 9, true); clf
[ha, pos] = tight_subplot(3, 3, .01, [.045, 0.03], [.065, .02]);
ha = reshape(ha', 3, [])';
freq_s = logspace(log10(1), log10(12500), 350)';

h1 = mimo_bode_mag(HH1, freq_s, ha, '-k');
h2 = mimo_bode_mag(HH2, freq_s, ha, '--r');

% fprintf('Bandwidth CZ: %f\n', bandwidth(HH1(1,1))/2/pi)
% fprintf('Bandwidth CR: %f\n', bandwidth(HH2(1,1))/2/pi)
%%
analyze_margins(plants, sys_obsDist, K_lqr1, L_dist, ha(1,1), '-m');
analyze_margins(plants, sys_obsDist, K_lqr2, L_dist, ha(1,1), '-b');

% frf_bode_mag(HH(ny,nu), freq_s, ha(ny,nu), 'Hz', varargin{:});
% frf_bode_mag(HH(ny,nu), freq_s, ha(ny,nu), 'Hz', varargin{:});

h1(1,1).DisplayName = 'choose-$\zeta$';
h2(1,1).DisplayName = 'constant-$\rho$';
leg = legend([h1(1,1), h2(1,1)]);

set(leg, 'location', 'southwest')

%%
save_fig(F11, fullfile(PATHS.thesis_root, 'plots_mpc_slf/figures/MIMO_CL'))



function hands = mimo_bode_mag(HH, freq_s, ha, varargin)
  hands = gobjects(3, 3);
  ylm = [-95, 30];
  coords = {'X', 'Y', 'Z'};
  for ny = 1:3
    for nu=1:3
      hands(ny, nu) = frf_bode_mag(HH(ny,nu), freq_s, ha(ny,nu), 'Hz', varargin{:});
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
      end
        set(ha(ny,nu), 'XTick', [10,100,1000,10000]);

    end
  end
  
end

function H_frd = get_H()
  %
  sysx_fn = fullfile(PATHS.sysid(), 'ALL-axis_sines_info_intsamps_quickFourierCoef_1-17-2019x-drive-01.json');
  sysy_fn = fullfile(PATHS.sysid(), 'ALL-axis_sines_info_intsamps_quickFourierCoef_1-17-2019y-drive-01.json');
  sysz_fn = fullfile(PATHS.sysid(), 'ALL-axis_sines_info_intsamps_quickFourierCoef_1-17-2019z-drive-03.json');
  sysx_ = SweptSinesOnline(sysx_fn);
  sysy_ = SweptSinesOnline(sysy_fn);
  sysz_ = SweptSinesOnline(sysz_fn);
  Gux = sysx_.FRF_from_FC(1, [2,3,4]);
  Guy = sysy_.FRF_from_FC(1, [2,3,4]);
  Guz = sysz_.FRF_from_FC(1, [2,3,4]);
  
  H_frd = [Gux, Guy, Guz];
  
end

function [DD, D_ki, D_inv] = get_gz_dz(Gzz)
  root = fullfile(PATHS.sysid, 'z_axis_evolution', 'batch-1');
  file = 'first_res_fit-3-17-2019-1.json';
  dat = loadjson(fullfile(root, file));
  
  ss_fname = fullfile(root, 'ALL-axis_sines_info_intsamps_quickFourierCoef_3-17-2019-01.json');
  ss_data = SweptSinesOnline(ss_fname);
  freqs = ss_data.freq_s;
  
%   G_frf = frd(ss_data.FC_s(:,2)./ss_data.FC_s(:,1), ss_data.freq_s, AFM.Ts, 'FrequencyUnit', 'Hz');
  
  idx1 = find(Gzz.Frequency < 180, 1, 'last');
  idx2 = find(Gzz.Frequency > 230, 1, 'first');
  
  
  gd0 = tf(dat.Dinv_den, dat.Dinv_Num, AFM.Ts)/dat.K
  sos_fos = SosFos(gd0);
  lg = LogCostZPK(Gzz.ResponseData(idx1:idx2), Gzz.Frequency(idx1:idx2)*2*pi, sos_fos);
  lg.solve_lsq(1)
  D = lg.sos_fos.realize();
  
  D_inv = 1/D;
  
  KI = -0.05;
  
%   D_inv = frd(tf(dat.Dinv_Num, dat.Dinv_den, AFM.Ts),  freqs,  'FrequencyUnit', 'Hz');
%   D_inv = tf(dat.Dinv_Num, dat.Dinv_den, AFM.Ts)
  D_ki = zpk([], [1], KI, AFM.Ts);
  
  DD = D_ki * D_inv;
  
  
  
%   Loop = DD_frf * G_frf;
%   Loop_noinv = D_ki_frf * G_frf;

end

function [HH] = close_three_axis(sys, sys_recyc, sys_obs, KxKu, LxLd, Gxyz_frf, Dy, Dki, Dinv)
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
  D2x = ss(AA_, LL_, K_c, 0, Ts);     % The feedback
  Mx = ss(AA_, BB_, -K_c, Nbar, Ts); % The feedforward

  % % Estimated state feedback loop gain:
  Loop = sys*D2x;
  Sens = 1/(1+Loop);

  Hyr = ((sys*Mx*Sens));
  Hyd = minreal( sys*Sens);
  
%    R-->[ M ]--->O---->[ D1 ]-->[ G ]--+---> y
%                 |                     |           
%                 +-----[ D2 ]<---------+
%
% y =        G*D1
%      --------------- * M * R
%      I + D2*G*G1


  if 1
      M = [Mx, 0, 0;
          0,  1, 0;
          0, 0, 1];
      
      DD2 = [D2x, 0, 0;
          0,   1, 0;
          0,   0, 1];
      DD1 = [1, 0, 0;
          0, Dy, 0;
          0, 0,  Dki*Dinv];
      
  else
      load('shaped_loop_x.mat');
      Dx = D_x.Dx;
      M = eye(3);
      
      DD2 = eye(3);
      
      DD1 = [Dx*0.5, 0, 0;
          0, Dy, 0;
          0, 0,  Dki*Dinv];
      
  end
  I = eye(3);

  
  HH = (Gxyz_frf*DD1 / (I + DD2*Gxyz_frf*DD1)) * M;
  
% % %   DD2 = [D2, 0, 0;
% % %          0,  Dy, 0;
% % %          0, 0, 0];
% % %   DD1 = [D1, 0, 0;
% % %          0, Dy, 0;
% % %           0, 0, Dki*Dinv];
% % %   I = eye(3);
% % %   
% % %   M1 = (Gxyz_frf*DD2);
% % %   M1(3,3) = Dki;
% % %   
% % %   HH = M1 / (I + Gxyz_frf*DD1);
% % % %   HH = Gxyz_frf / (I + Gxyz_frf*DD1);
% % %   


end

function analyze_margins(plants, sys_obsDist, K_lqr, L_dist, ax, varargin)
  [Sens, Hyd, Hyr, Hyeta, Loop] = ss_loops_delta_dist(plants.SYS, plants.sys_recyc,...
    sys_obsDist, K_lqr, L_dist);
  
  F_clbode = figure(25);clf; hold on, grid on
  Hbode_sens = bodeplot(Sens);
  setoptions(Hbode_sens, 'FreqUnits', 'Hz')
  legend('S: linear')
  grid on, hold on;
  
  [Gm_lin, Pm_lin] = margin(Loop);
  
  fprintf('-------- MARGINS ------------------\n')
  fprintf('Linear: GM = %.2f [dB], PM = %.2f [deg]\n', 20*log10(Gm_lin), Pm_lin)
  fprintf('Bandwidth: %f [Hz]\n', bandwidth(minreal(Hyr))/2/pi);
  
  if nargin > 4
      ws = logspace(log10(1), log10(12500), 300);
      frf_bode_mag(Hyr, ws, ax, 'Hz', varargin{:});
  end

  
end
