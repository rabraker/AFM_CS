

function [N_extra] = mu_overscan(Hyr, x_rate, mu_Nsamples, verbose, N_prescan)
% [N_extra] = mu_overscan(Hyr, x_rate, mu_Nsamples, verbose, N_prescan)
% 
%  Hyr: LTI model of the closed loop.
%  x_rate: scan speed, in volts per sample period.
%  mu_Nsamples: 

  
  Ts = Hyr.Ts;
  
  % TF from reference to error.
  Her = minreal(1 - Hyr); %%feedback(1, D*G);
  Int_z = zpk([], [1], 1, Ts);
  % Compute the steady state error to a ramp 
  ess = dcgain(minreal(Her*Int_z))*x_rate;
  
  % sometimes, the integrator gets perturbed to a really slow pole, so dc-gain 
  % is infinitate. So evaluate low-freq gain instead.
  if isinf(ess)
    ess = abs(freqresp(Her*Int_z, 1))*x_rate;
  end
  
  N = mu_Nsamples+N_prescan;
  N_extra = floor(ess/x_rate);

  x_N = (N-1)*x_rate;
  r_vec = linspace(0, x_N, N);
  t_vec = [0:1:N-1]'*Ts;
  y_vec = lsim(Hyr, r_vec, t_vec);
  
  x_N_extra = (N+N_extra-1)*x_rate;
  r_vec_extra = linspace(0, x_N_extra, N+N_extra);
  t_vec_extra = [0:1:N+N_extra-1]'*Ts;
  y_vec_extra = lsim(Hyr, r_vec_extra, t_vec_extra);
  
  % When we do the prescan, we reach quasi steady state by the end, typically.
  % This means that we dont necessarily need the all the overscan samples,
  % because what we are really interested in is getting a scan that spans at
  % mu_len microns.

  % X0 is the start of the mu-path scan, which happens after the last prescan
  % sample.
  x0 = y_vec_extra(N_prescan+1);
  % X1 is the mu_path
  mu_len_volts = mu_Nsamples*x_rate;
  x1 = x0+mu_len_volts;
  
  x1_idx = find(y_vec_extra >= x1, 1, 'first');
  if isempty(x1_idx)
    error('x1_idx is empty. Fix this corner case.')
  end
  % Then the *actual* number of samples to overscan is 
  N_extra = x1_idx - (mu_Nsamples + N_prescan);
  
  
  if exist('verbose', 'var') && verbose
    % figure(100); clf
    fig = mkfig(100, 5, 4); clf
    ha = tight_subplot(1, 1, 0, [.06, .02], [.26, .02]);
    
    plot(N*Ts, r_vec(end), 'x')
    hold on
    plot(t_vec_extra, y_vec_extra, t_vec_extra, r_vec_extra)
    plot(t_vec, y_vec, '--')
    
    grid
    xlabel('time [s]')
    ylabel('y [v]')
    ylim([0, 1.2*max(r_vec)]);
    
    plot([N_prescan, N_prescan]*Ts, [0, r_vec(N_prescan)], '--k')
    plot([N_prescan+mu_Nsamples, N_prescan+mu_Nsamples]*Ts,...
      [0, r_vec(N_prescan+mu_Nsamples)], '--k')
    
    % X0 is the start of the mu-path scan, which happens after the last prescan
    % sample.
    x0 = y_vec_extra(N_prescan+1);
    % X1 is the mu_path
    mu_len_volts = mu_Nsamples*x_rate;
    x1 = x0+mu_len_volts;
    
    x1_idx = find(y_vec_extra >= x1, 1, 'first');
    
    plot([0, N_prescan+1]*Ts, [x0, x0], '--k')
    plot([0, x1_idx]*Ts, [x1, x1], '--k')
    
    % set(gca(), 'YTick', [x0, x1])
    set(gca(), 'YTickLabel', [], 'XTickLabel', [])
    
    txt_opts_xl = {'Units', 'Data', 'HorizontalAlignment', 'center', 'FontSize', 12};
    txt_opts_yl = {'Units', 'Data', 'HorizontalAlignment', 'right', 'FontSize', 12};
    
    text(-20*Ts, x1, '$x(t_{prescan}+t_{scan})$', txt_opts_yl{:})
    text(-20*Ts, x0, '$x(t_{prescan})$', txt_opts_yl{:})
    text(-20*Ts, 0, '$0$', txt_opts_yl{:})
    
    yly=-.005;
    text(N_prescan*Ts, yly, '$t_{prescan}$', txt_opts_xl{:})
    
    s=sprintf('$t_{prescan} +$ \n $t_{scan}$');
    text((N_prescan+mu_Nsamples)*Ts, yly-0.005, s, txt_opts_xl{:})
%     figure(100)
%     plot(t_vec, y_vec)
%     hold on
%     plot(mu_Nsamples*Ts, r_vec(end), 'x')
%     hold on
%     plot(t_vec_extra, y_vec_extra, t_vec_extra, r_vec_extra)
%     
%     grid
%     xlabel('time [s]')
%     ylabel('y [v]')
  end
  
  
end