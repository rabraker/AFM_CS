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


function [HH] = close_three_axis(Gxyz_frf, xdir_cntrl, Dy, Dzki, Dinv, Dx_ff, Dy_ff)
  
  D1x = xdir_cntrl.D1_loop; % The feedback
  Mx = xdir_cntrl.D2_ss_ff; % The feedforward
  
  %    R-->[ M ]--->O---->[ D2 ]-->[ G ]--+---> y
  %                 |                     |
  %                 +-----[ D1 ]<---------+
  %
  % y =        G*D1
  %      --------------- * M * R
  %      I + D2*G*G1
  
  M = [Mx*Dx_ff, 0, 0;
      0,  Dy_ff, 0;
      0, 0, 1];
  
  DD1 = [D1x, 0, 0;
      0,   1, 0;
      0,   0, 1];
  DD2 = [1, 0, 0;
      0, Dy, 0;
      0, 0,  Dzki*Dinv];

  I = eye(3);
  
  HH = (Gxyz_frf*DD2 / (I + DD1*Gxyz_frf*DD2)) * M;
  
end