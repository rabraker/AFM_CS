
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

% ------- Load Plants -----
[plants, frf_data] = CanonPlants.plants_ns14(9, '5micron');

Ts  = plants.SYS.Ts;
Dy = zpk([0], [1], 0.01, AFM.Ts);

du_max_orig = StageParams.du_max;
du_max = du_max_orig/norm(plants.gdrift_inv, Inf);

xdirControl = get_xdir_standard_control('const-sig')

%%

cmplx_rad = 0.9;
[Q2, R2, S2] = build_control_constsigma(plants.sys_recyc, cmplx_rad);

% ------------------------- Observer Gain ---------------------------------
can_obs_params = CanonObsParams_01();
can_obs_params.beta = 50;
[sys_obsDist, L_dist] = build_obs(plants.SYS, can_obs_params);

gam_rob2 = 46.4;
K_lqr2 = dlqr(plants.sys_recyc.a, plants.sys_recyc.b, Q2, R2+gam_rob2, S2);

%%
Gxyz_frd = get_H();
xdir_cntrl = get_xdir_standard_control('const-sig');

w = 300*2*pi;
d = c2d(tf(w, [1 w]), AFM.Ts, 'matched');

Dx_ff = make_dx_ff()*d;
Dy_ff = make_dy_ff()*d;

Gzz = Gxyz_frd(3,3);
[DDz, Dki, Dinv] = get_gz_dz(Gzz);

[HH1] = close_three_axis(Gxyz_frd, xdir_cntrl, Dy, Dki, Dinv, 1, 1);
[HH2] = close_three_axis(Gxyz_frd, xdir_cntrl, Dy, Dki, Dinv, Dx_ff, Dy_ff);
 
F11 = mkfig(10, 7, 5, true); clf
% [ha, pos] = tight_subplot(3, 3, .01, [.061, 0.03], [.065, .02]);
[ha, pos] = tight_subplot(3, 2, .01, [.061, 0.03], [.065, .02]);
ha = reshape(ha', [], 3)';
freq_s = logspace(log10(1), log10(12500), 350)';

h1 = mimo_bode_mag(HH1(:,1:2), freq_s, ha, '--r');
h2 = mimo_bode_mag(HH2(:, 1:2), freq_s, ha, '-b');
h3 = mimo_bode_mag(Gxyz_frd(:, 1:2), freq_s, ha, ':k');

h4 = frf_bode_mag(Dx_ff, freq_s, ha(1,1), 'Hz', 'm');
h5 = frf_bode_mag(Dy_ff, freq_s, ha(2,2), 'Hz', 'm');
ylabel(ha(2,2), '')


% fprintf('Bandwidth CZ: %f\n', bandwidth(HH1(1,1))/2/pi)
% analyze_margins(plants, sys_obsDist, K_lqr2, L_dist, ha(1,1), '-b');
% frf_bode_mag(HH(ny,nu), freq_s, ha(ny,nu), 'Hz', varargin{:});
% frf_bode_mag(HH(ny,nu), freq_s, ha(ny,nu), 'Hz', varargin{:});

h1(1,1).DisplayName = 'C.L w/o F.F';
h2(1,1).DisplayName = 'C.L w/ F.F';
h3(1,1).DisplayName = 'O.L. Plant';
h4.DisplayName = '$D_{x,ff}$';
h5.DisplayName = '$D_{y,ff}$';
%%
leg = legend([h1(1,1), h2(1,1), h3(1,1), h4], 'FontSize', 10);
legy = legend(h5, 'FontSize', 12);
set(leg, 'location', 'southwest')
set(legy, 'location', 'southwest')

%%
save_fig(F11, fullfile(PATHS.tmech_fig(), 'MIMO_CL_uxuy'))


function Dy_ff = make_dy_ff()
    wz1 = 214 * 2 * pi;
    wz2 = 505 * 2 * pi;
    zz1 = 0.01;
    zz2 = 0.04;
    
    g = tf([1, 2*zz1*wz1, wz1^2], conv([1, 300*2*pi],  [1, 300*2*pi]));
    g2 = tf([1, 2*zz2*wz2, wz2^2], conv([1, 450*2*pi],  [1, 450*2*pi]));
    
    g = c2d(g, AFM.Ts, 'matched') * c2d(g2, AFM.Ts, 'matched');
    Dy_ff = tf(g) * (1/dcgain(g))
end

function Dx_ff = make_dx_ff()
    wz1 = 214 * 2 * pi;
    wz2 = 505 * 2 * pi;
    zz = 0.001;
    
    
    g = tf([1, 2*zz*wz1, wz1^2], conv([1, 300*2*pi],  [1, 300*2*pi]));
    g2 = tf([1, 2*zz*wz2, wz2^2], conv([1, 450*2*pi],  [1, 450*2*pi]));
    
    g = c2d(g, AFM.Ts, 'matched') * c2d(g2, AFM.Ts, 'matched');
    Dx_ff = tf(g) * (1/dcgain(g))
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


function hands = mimo_bode_mag(HH, freq_s, ha, varargin)
  hands = gobjects(3, 3);
  ylm = [-95, 30];
  coords = {'X', 'Y', 'Z'};
  for ny = 1:size(HH, 1)
    for nu=1:size(HH, 2)
        try
      hands(ny, nu) = frf_bode_mag(HH(ny,nu), freq_s, ha(ny,nu), 'Hz', varargin{:});
        catch
            keyboard
        end
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
    % Load initial guess.
    root = fullfile(PATHS.sysid, 'z_axis_evolution', 'batch-1');
    file = 'first_res_fit-3-17-2019-1.json';
    dat = loadjson(fullfile(root, file));
    
    idx1 = find(Gzz.Frequency < 180, 1, 'last');
    idx2 = find(Gzz.Frequency > 230, 1, 'first');
    % Initial guess
    gd0 = tf(dat.Dinv_den, dat.Dinv_Num, AFM.Ts)/dat.K;
    sos_fos = SosFos(gd0);
    lg = LogCostZPK(Gzz.ResponseData(idx1:idx2), Gzz.Frequency(idx1:idx2)*2*pi, sos_fos);
    lg.solve_lsq(1)
    D = lg.sos_fos.realize();
    
    D_inv = 1/D;
    
    KI = -0.05;
    
    D_ki = zpk([], [1], KI, AFM.Ts);
    
    DD = D_ki * D_inv;
    
end