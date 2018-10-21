clear, clc

dry_run = false;
saveon = false;

Ts = 40e-6;
if 0
  u_max = 1.75;
  n_space = 8000;
  n_up = 5;
  step_sz = u_max/n_up;
  
  imps = [0;step_sz*ones(n_up, 1); -step_sz*ones(2*n_up-1,1)];
  
  
  for k = 2*n_up-2:-1:1
    if mod(k,2) == 1
      sgn = -1;
    else
      sgn = 1;
    end
    imps = [imps; sgn * step_sz * ones(k, 1)];
    
  end
  % imps = [imps; -sum(imps)]
  % ref_s = cumsum(imps)
  N_imp = length(imps)
  
  
  impulse_idx= (1:n_space:N_imp*n_space)';
  u_vec = zeros((N_imp)*(n_space), 1);
  u_vec(impulse_idx) = imps;
  u_vec = cumsum(u_vec);
  
  t_vec = (0:length(u_vec)-1)'*Ts;
  
% lpf
  w1 = 100*2*pi;
  F = tf(w1, [1, w1])
  F = c2d(F*F, Ts);
  u_vec = lsim(F, u_vec, t_vec);
  
  figure(1); clf
  plot(t_vec, u_vec);
  grid on

else
  u_max = 1.75;
  n_space = 16000;
  n_up = 5;
  
  To = Ts*n_space;
  wo = 2*pi/To;
  t_vec = (0:n_up*n_space-1)'*Ts;
  u_vec = sign(sin(wo*t_vec))*.75 + .5;

end

  figure(1); clf
  plot(t_vec, u_vec);
  grid on

%%

% reset_piezo('t1', 15, 't_final', 25, 'umax', 10, 'k1', 0.55,...
%             'verbose', true, 'dry_run', dry_run)
          

umax = 10;
clear vi;
vipath_reset = 'C:\Users\arnold\Documents\afm-cs\labview\play_nPoint_id_slew_OL_FIFO.vi';
[e, vi] = setupVI(vipath_reset, 'Abort', 0,...
  'umax', umax, 'TsTicks', 1600, 'u_in', u_vec, 'DAC-select', 2);
vi.Run;
stage_dat = vi.GetControlValue('stage_data_out');

u_exp = stage_dat(:,1);
yx_exp = stage_dat(:,4); % - dat(1,2);
t_exp = (0:length(u_exp)-1)'*Ts;

figure(201); clf
plot(t_exp, -u_exp)
hold on
plot(t_exp, yx_exp - yx_exp(1))
grid on

%%
if saveon +1
  hystData.t_exp = t_exp;
  hystData.u_exp = u_exp;
  hystData.y_exp = yx_exp;
  hystData.umax = umax;
  hystData.impulse_idx = impulse_idx;
  % hystData.u_reset = u_reset;
  data_path = fullfile(PATHS.sysid, ['hysteresis/hystID_data_', date(), '_03.mat']);
  fprintf('\n\n--------------------------\n');
  fprintf('data path: %s\n', data_path);
  save(data_path, 'hystData')
end

%%

idd = iddata(yx_exp - yx_exp(1), u_exp, Ts);

g = tfest(idd, 12, 'Ts', Ts)
%%
figure(201)
y_sim = lsim(g, u_exp, t_exp);

plot(t_exp, y_sim)

