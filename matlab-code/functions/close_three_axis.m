% [HH] = close_three_axis(Gxyz_frf, xdir_cntrl, Dy, Dzki, Dinv, Dx_ff, Dy_ff)  
%
% Construct the entire closed transfer function for all three afm axes.
% 
% Arguments
% ----------
%   Gxyz_frf:  An FRD object of the full MIMO frequency response;
%   xdir_cntrl: a struct of the x-direction state-space stuff, obtain this from 
%               get_xdir_standard_control();
%   Dy:         This should be a simple PI compensator (for now, the code works
%               fine with a more general TF, but labview doesnt implement that yet.
%   Dzki:      Z-direction PI controller
%   Dinv:      Z-direction bending mode inversion
%   Dx_ff:     X-direction feedforward
%   Dy_ff:     Y-direction feedforward
%
%  The resulting compensator will have the form (all blocks are transfer
%  matrices):
% 
%    R-->[ M ]--->O---->[ D2 ]-->[ G ]--+---> y
%                 |                     |
%                 +-----[ D1 ]<---------+
%
% y =        G*D1
%      --------------- * M * R
%      I + D2*G*G1


function [HH] = close_three_axis(Gxyz_frf, xdir_cntrl, ydir_cntrl, Dzki, Dinv)
  
  Dssx = xdir_cntrl.D1_loop; % The feedback
  Mx = xdir_cntrl.D2_ss_ff; % The feedforward
  Dx_ff = xdir_cntrl.D_ff;
  Dx = xdir_cntrl.D;       % The transfer-function feedback (if any)
  Dki_x = xdir_cntrl.D_ki;
  
  Dssy = ydir_cntrl.D1_loop; % The state-space feedback (if any)
  My = ydir_cntrl.D2_ss_ff; % The state-space feedforward (if any)
  Dy_ff = ydir_cntrl.D_ff; % the real feedforward. (if any)
  Dy = ydir_cntrl.D;       % The transfer-function feedback (if any)
  Dki_y = ydir_cntrl.D_ki;
  %    R-->[ M ]--->O---->[ DD ]-->[ G ]--+---> y
  %                 |                     |
  %                 +-----[ Dss ]<---------+
  %
  % y =        G*D1
  %      --------------- * M * R
  %      I + D2*G*G1
  % 
  % For state-space, DD = 1, Dss and M are some transfer function.
  % FOr typical TF based feedback, Dss=1, and M=1, unless we designed a
  % feedforward control, like Dx_ff.
  
  
  M = [Mx*Dx_ff, 0, 0;
      0,  My*Dy_ff, 0;
      0, 0, 1];
  
  Dss = [Dssx, 0, 0;
      0,   Dssy, 0;
      0,   0, 1];
  DD = [Dx*Dki_x, 0, 0;
      0, Dy*Dki_y, 0;
      0, 0,  Dzki*Dinv];

  I = eye(3);
  
  HH = (Gxyz_frf*DD / (I + Dss*Gxyz_frf*DD)) * M;
  
end