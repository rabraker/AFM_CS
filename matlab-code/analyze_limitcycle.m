clc
clear
pix = 256;
width = 5;
microns_per_volt = 50/10;
pix_per_volt = (pix/width)*microns_per_volt;
Ts = 40e-6;
addpath('functions')
% data_root = 'C:\Users\arnold\Documents\labview\afm_imaging\data\';
% data_root = fullfile(getdataroot(), 'cs-data');
data_root = fullfile(getdataroot(), 'cs-data\20microns\9-23-2017');
data_root = '/media/labserver/afm-cs'
% ---------------------------------------------------


% f = 'cs-traj-z-bounce2.csv';
% fpath = fullfile(data_root, f);
% dat = csvread(fpath);



cs_exp_data_name = 'cs-traj-z-bounce_out_10-8-2018-13.csv'; % illustrates problem
% cs_exp_data_name = 'cs-traj-z-bounce_out_10-8-2018-22.csv';
cs_exp_meta_name = strrep(cs_exp_data_name, '.csv', '-meta.mat');

cs_data_path = fullfile(data_root, cs_exp_data_name);
cs_meta_path = fullfile(data_root, cs_exp_meta_name);


dat = csvread(cs_data_path); 

load(cs_meta_path);  % Provides ExpMetaData

fprintf('loading done...\n')

indc = {'k',        'r',       [0, .75, .75], 'm', [.93 .69 .13], 'b';
       'xy-move', 'tip down', 'tip settle',  'na', 'tip up', '$\mu$-path scan'};

xy = false;
figbase = 100;

x = dat(:,1);
y=dat(:,2);
t = [0:1:length(x)-1]'*Ts;

z_err = dat(:,3);
uz = dat(:,4);
met_ind = dat(:,5);



Fig_uz = figure(20+figbase); clf
% plot(uz)
ax2 = gca
plotbyindex(ax2, t, uz, met_ind, indc);
title('uz')

Fig_ze = figure(30+figbase); clf
% plot(z_err)
ax3 = gca();
hp = plotbyindex(ax3, t, z_err, met_ind, indc);
title('z-err')
hold on
plot([t(1), t(end)], [.05, .05])
plot([t(1), t(end)], -[.05, .05])
legend(hp(2:6))

if xy
 
  figure(10+figbase), clf
  % plot(x)
  title('x')
  ax1 = gca;
  plotbyindex(ax1, t, x, met_ind, indc);
  
  figure(40+figbase), clf
  ax4=gca();
  plotbyindex(ax4, t, y, met_ind, indc);
  title('y')
  linkaxes([ax1, ax2, ax3, ax4], 'x')
else
  linkaxes([ax2, ax3], 'x')
end



% figure(40)
% plot(min(met_ind, 3))
% title('meta index')
% ax4 = gca();
% linkaxes([ax1, ax2, ax3, ax4], 'x')



clc
Ts = 40e-6;
tdown_idx = find(met_ind == -5);

indc = {'k',        'r',       [0, .75, .75], 'm', [.93 .69 .13], 'b';
       'xy-move', 'tip down', 'tip settle',  'na', 'tip up', '$\mu$-path scan'};

met_ind0 = met_ind(2:end);
met_ind0(met_ind0 >0) = -4;
met_ind0(met_ind0 ==0) = -4; % for data corruption.

%%
figure(32), clf
% yyaxis right
plot(met_ind0, 'g')
grid on, hold on
%
clc
typ = -5;
met_ind_temp = abs(met_ind0);
zerr_s = {};
t_s = {};
%          -1       -2       -3       -4      -5
names = {'move', 'tdown', 'tsettle', 'scan', 'tup'};
names_idx = [1:5];
idx_state_s.scan = {};
idx_state_s.move = {};
idx_state_s.tdown = {};
idx_state_s.tsettle = {};
idx_state_s.tup = {};

% figure(31);
% yyaxis right
clrs = {'k', 'r', 'g', 'b', 'm'};
idx_end = 0;

break_out = false;
for k=1:50
  
  for j = 1:length(names)
    idx_end_prev = idx_end;
    j_next = max(mod(j+1, 6), 1); % wrap back to 1 when j=5.
    idx_start = find(met_ind_temp(idx_end_prev+1:end) == names_idx(j), 1, 'first') + idx_end_prev;
    
     plot(idx_start, met_ind0(idx_start), 'o', 'MarkerSize', 10, 'Color', clrs{j})
    
    idx_end   = find(met_ind_temp(idx_end_prev+1:end) == names_idx(j_next), 1, 'first') + idx_end_prev-1;
    
    if isempty(idx_end)
      idx_end   = find(met_ind_temp(idx_end_prev+1:end) == names_idx(j), 1, 'last') + idx_end_prev;
    end
    
    plot(idx_end, met_ind0(idx_end), 'x', 'MarkerSize', 10, 'Color', clrs{j})
    
    if isempty(idx_end) || idx_end >=length(met_ind_temp)
      break_out = true;
      break
    end
    idx_state_s.(names{j}){k} = idx_start:idx_end;
   
    
  end
 
  if break_out
    break
  end
  
end

%
% ---------------- Visualize uz engagement point. -------------------------

for k=1:length(idx_state_s.tdown)
  
  ze_idx = idx_state_s.tdown{k};
  t_ = t(ze_idx);
  ze = z_err(ze_idx);
  uz_ = uz(ze_idx);
  
  [~, idx_min] = min(ze);
  figure(Fig_ze);
  plot(t_(idx_min), ze(idx_min), 'x')
  
  figure(Fig_uz);
  plot(t_(idx_min), uz_(idx_min), 'x')
  
  
  
end


legend(hp(2:6))



%%
% ---------------------------- Test out variance detector -----------------
figure(300), clf
hold on, grid on

ax3 = gca();
hp = plotbyindex(ax3, t, z_err, met_ind, indc);
title('z-err')
hold on
plot([t(1), t(end)], [.05, .05])
plot([t(1), t(end)], -[.05, .05])

ze_idx = idx_state_s.tup{1};

t_ = t(ze_idx);
ze = z_err(ze_idx);
var_std = sqrt(var(ze(end-100:end)))*2


TOL = var_std;

N = 32;


found = false;
for j=1:length(idx_state_s.tup)
  ze_idx = idx_state_s.tup{j};
  t_ = t(ze_idx);
  ze = z_err(ze_idx);

%   plot(t_, ze)
  for k = N:length(t_)
    
    ze_k = ze(k-N+1:k);
    
    std_k = sqrt(var(ze_k));
    
    std_s(k-N+1) = std_k;
    
    
    mu = mean(ze_k);
    if std_k < TOL
      plot([t_(k-N+1), t_(k)], [mu, mu], 'b')
      
      errorbar(t_(k-N+1), mu, std_k)
      errorbar(t_(k), mu, std_k)
      found = true;
      break
    else
      found = false;
    end
  end

  if j==5
%       keyboard
  end

  if ~found
    [~, idx_min] = min(ze);
    
    plot(t_(idx_min), ze(idx_min), 'x')
  end

end

legend(hp(2:6))







