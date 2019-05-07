addpath functions
close all
root = fullfile(PATHS.sysid, 'z_axis_evolution', 'batch-1');
files = {...
'first_res_fit-3-4-2019.json',...
'first_res_fit-3-9-2019-2.json',...
'first_res_fit-3-17-2019-1.json',...
};
files = fullfile(root, files);

file_nclr = 'first_res_fit-4-01-2019-mica-nclr.json';
files{end+1} = fullfile(PATHS.sysid, 'z_axis_evolution', file_nclr);

% 'first_res_fit-3-11-2019-1.json',...
% 'first_res_fit-3-5-2019.json',...
% 'first_res_fit-3-7-2019.json',...
names = {'(A) SICON, L=450~$\mu$m',...
  '(B: Mar. 9) SICON, L=450~$\mu$m',...
  '(B: Mar. 17) SICON, L=450~$\mu$m',...
  '(C) NCLR, L=225~$\mu$m'};


Fig = mkfig(1, 5, 4); clf
[ha, pos] = tight_subplot(1, 1, [.02, .01 ], [.14, 0.03], [.14, .05]);

h = gobjects(length(files), 1);
for k=1:length(files)
  if k==1
    prefix='(A) ';
  else
    prefix = '(B) ';
  end
  name = strrep(files{k}, 'first_res_fit-', '');
  name = strrep(name, '-1.json', '');
  name = strrep(name, '-2.json', '');
  name = strrep(name, '.json', '');
  dat = loadjson(files{k});

  
  gz_k = dat.SOS_frf_data.resp_real + dat.SOS_frf_data.resp_imag*1i;
  h(k) = frf_bode_mag(gz_k, dat.SOS_frf_data.omegas/2/pi, Fig, 'Hz');
  h(k).DisplayName = names{k};
end

leg0 = legend(h);
set(leg0, 'Position', [0.1493 0.1554 0.4665 0.1643]); 


save_fig(Fig, fullfile(PATHS.thesis_fig_final, 'z_cant_evolution'))
%%
names = {'(A) L=450~$\mu$m',...
  '(B: Mar. 9) L=450~$\mu$m',...
  '(B: Mar. 17) L=450~$\mu$m',...
  '(C) L=225~$\mu$m'};

for k=1:length(names)
    h(k).DisplayName = names{k};
end
ha.FontSize = 14;
leg0.FontSize = 12;
leg0.Position = [0.1394 0.1370 0.4607 0.2143];
%%

save_fig(Fig, fullfile(PATHS.defense_fig(), 'z_cant_evolution'), true)







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


