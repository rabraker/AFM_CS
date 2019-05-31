clear all
% clc
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
Gxyz_frd = get_H_frd();
size(Gxyz_frd)
% Gxyz_frd = Gxyz_frd(:,:, 1:end-1)


%%

ydir_cntrl = get_ydir_standard_control();


%%
[plants, frf_data] = CanonPlants.plants_ns14(9);
newest_xfrf_path = fullfile(PATHS.sysid, 'x-axis_sines_info_first_resFourierCoef_5-9-2019-01.json');
sysx_new = SweptSinesOnline(newest_xfrf_path);
%%
[z, p, k] = zpkdata(plants.Gvib, 'v');
p1 = p(end-1:end);
z1 = z(end-3:end-2);
p(end-1:end) = [];
z(end-3:end-2) = [];

g_tmp = zpk(z, p, k, plants.Gvib.Ts)*plants.gdrift;
g_tmp = g_tmp*dcgain(plants.Gvib*plants.gdrift)/dcgain(g_tmp);

Gx_bend_ = sysx_new.FRF_from_FC(1, [2]);
Gx_bend = Gx_bend_/g_tmp;

figure(100)
opts = bodeoptions();
opts.FreqUnits = 'Hz';
bodeplot(Gx_bend, opts)
grid

gx_bend0 = zpk(z1, p1, 1, plants.Gvib.Ts);
%%
sos_fos_xbend = SosFos(gx_bend0, 'iodelay', 9)
lg = LogCostZPK(squeeze(Gx_bend.Response), Gx_bend.Frequency*2*pi, sos_fos_xbend);
lg.solve_lsq(1);
gx_bend_lg = lg.sos_fos.realize();

hold on
bodeplot(gx_bend_lg, Gx_bend.Frequency*2*pi)

%%
Gxyz = mf_full.modelFit.G_all_axes(1:3,:);

% ------- Load Plants -----
[plants, frf_data] = CanonPlants.plants_ns14(9);

[z, p, k] = zpkdata(plants.Gvib, 'v');

figure(1); clf
opts = bodeoptions();
opts.FreqUnits = 'Hz';
h = bodeplot(plants.Gvib, {1*2*pi, 12.5e3*2*pi}, opts);
grid on

% figure(2);
% pzplot(plants.Gvib)

p1 = p(end-1:end);
z1 = z(end-3:end-2);

D_notch= zpk([0.979 + 1i*0.163, 0.979  - 1i*0.163], [0.8, 0.8], 1, AFM.Ts);
D_notch = D_notch/dcgain(D_notch);

Dinv = zpk(p1, z1, 1, AFM.Ts);
Dki = zpk(0, 1, 0.03, AFM.Ts);

Loop = plants.Gvib*Dinv*Dki*D_notch;

figure(1)
hold on
bodeplot(Loop, opts)
bodeplot(1/(1+Loop), opts)
grid on


figure(3)
H = feedback(ss(Loop), 1);

bodeplot(H, {1*2*pi, 12.5e3*2*pi}, opts)
grid on

figure(2)
step(ss(H))

H_ruz = ss(Dinv*Dki*D_notch)/ss(1+(ss(Loop)));


figure(6)
step(ss(H_ruz))

[gm, pm] = margin(Loop);
gm = 20*log10(gm);

fprintf('bandwidth: %.2f[Hz], PM: %.2f [deg], GM: %.2f [dB]\n', bandwidth(H)/2/pi, pm, gm)
%
u = raster(15, AFM.Ts, 3);
y = lsim(ss(H), u.Data, u.Time);

figure(5);clf
plot(u.Time, u.Data)
hold on
plot(u.Time, y, '--')
%%
figure(60); clf;
[ha, pos] = tight_subplot(3, 3, .01, [.045, 0.03], [.065, .02]);

ha = reshape(ha', 3, [])';

[Dz, Dz_ki, Dz_inv] = get_gz_dz(Gxyz_frd(3,3));
% Dz = Dz
Dx = Dinv*Dki*D_notch;
H3 = close_three_axis(Gxyz_frd, Dx, Dy, Dz);

hands = mimo_bode_mag(H3, H3.Frequency, ha, '-r')
%%

    


function [HH] = close_three_axis(Gxyz_frf, Dx, Dy, Dz)

  
%    R-->[ M ]--->O---->[ D1 ]-->[ G ]--+---> y
%                 |                     |           
%                 +-----[ D2 ]<---------+
%
% y =        G*D1
%      --------------- * M * R
%      I + D2*G*G1

M = eye(3);

DD2 = eye(3);

DD1 = [Dx*1, 0, 0;
    0, Dy, 0;
    0, 0,  Dz];
  
  I = eye(3);
  HH = (Gxyz_frf*DD1 / (I + DD2*Gxyz_frf*DD1)) * M;
  
end


function [DD, D_ki, D_inv] = get_gz_dz(Gzz)
  root = fullfile(PATHS.sysid, 'z_axis_evolution', 'batch-1');
  file = 'first_res_fit-3-17-2019-1.json';
  dat = loadjson(fullfile(root, file));
  
  ss_fname = fullfile(root, 'ALL-axis_sines_info_intsamps_quickFourierCoef_3-17-2019-01.json');
  ss_data = SweptSinesOnline(ss_fname);

  idx1 = find(Gzz.Frequency < 180, 1, 'last');
  idx2 = find(Gzz.Frequency > 230, 1, 'first');
  
  gd0 = tf(dat.Dinv_den, dat.Dinv_Num, AFM.Ts)/dat.K
  sos_fos = SosFos(gd0);
  lg = LogCostZPK(Gzz.ResponseData(idx1:idx2), Gzz.Frequency(idx1:idx2)*2*pi, sos_fos);
  lg.solve_lsq(1)
  D = lg.sos_fos.realize();
  
  D_inv = 1/D;
  
  KI = -0.06;
  
  D_ki = zpk([], [1], KI, AFM.Ts);
  
  DD = D_ki * D_inv;
  
end

function H_frd = get_H_frd()
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

function hands = mimo_bode_mag(HH, freq_s, ha, varargin)
  hands = gobjects(3, 3);
  ylm = [-95, 30];
  coords = {'X', 'Y', 'Z'};
  for ny = 1:3
    for nu=1:3
      hands(ny, nu) = frf_bode_mag(HH(ny,nu), freq_s, ha(ny,nu), 'Hz', varargin{:});
      set(ha(ny,nu), 'YLim', ylm);
      hold(ha(ny, nu), 'on')
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