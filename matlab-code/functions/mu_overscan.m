

function [N_extra] = mu_overscan(Hyr, x_rate, mu_Nsamples, N_prescan)
% [N_extra] = mu_overscan(Hyr, x_rate, mu_Nsamples, verbose, N_prescan)
% 
%  Hyr: LTI model of the closed loop.
%  x_rate: scan speed, in volts per sample period.
%  mu_Nsamples: 

  
% This is hueristic and will fail for super short mu-paths. But basically, the
% idea is to simulate a long enough ramp, and compute how many extra samples are
% needed to actually travel for the length of 1 mu-path.
    mu_len = mu_Nsamples * x_rate;
    N_ = 4*mu_Nsamples + N_prescan;
    
    t_vec = (0:N_-1)'*AFM.Ts;
    r_vec = (0:N_-1)*x_rate;
    y_vec = lsim(Hyr, r_vec, t_vec);
    
    y_scan = y_vec(N_prescan+1:end);
    y_scan = y_scan - y_scan(1);
    
    
    last_idx = find(y_scan > mu_len, 1, 'first');

    N_extra = last_idx - mu_Nsamples;

    return

  
end