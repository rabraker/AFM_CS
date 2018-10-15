

function [N_extra] = mu_overscan(G, D, x_rate, mu_Nsamples, verbose)
  
  
  Ts = G.Ts;
  Hyr = feedback(D*G, 1);
  
  % TF from reference to error.
  Her = feedback(1, D*G);
  Int_z = zpk([], [1], 1, G.Ts);
  % Compute the steady state error to a ramp 
  ess = dcgain(minreal(Her*Int_z))*x_rate;
  
  N = mu_Nsamples;

  
  
  N_extra = floor(ess/x_rate)


  if exist('verbose', 'var') && verbose
    x_N = (N-1)*x_rate;
    r_vec = linspace(0, x_N, N);
    t_vec = [0:1:N-1]'*Ts;
    y_vec = lsim(Hyr, r_vec, t_vec);
    
    x_N_extra = (N+N_extra-1)*x_rate;
    r_vec_extra = linspace(0, x_N_extra, N+N_extra);
    t_vec_extra = [0:1:N+N_extra-1]'*Ts;
    y_vec_extra = lsim(Hyr, r_vec_extra, t_vec_extra);
    
    figure(100)
    plot(t_vec, y_vec)
    hold on
    plot(mu_Nsamples*Ts, r_vec(end), 'x')
    hold on
    plot(t_vec_extra, y_vec_extra, t_vec_extra, r_vec_extra)
    
    grid
    xlabel('time [s]')
    ylabel('y [v]')
  end
  
  
end