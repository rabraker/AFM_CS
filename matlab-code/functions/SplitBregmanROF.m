% function img_filt = SplitBregmanROF(img, lambda, tol)
%
% Performs anisotropic denoising on the image in matrix img.
%
% Solves the optimization problem
%   min ||\nabla_x U||_1 + ||\nabla_y U||_1 + (\mu/2) ||U - img||_2
%    U
% by recasting it as the unconstrained problem
% 
%    min     |d_x|_1 + |d_y|_1 + (\mu/2)|U-f|_2 + 
% U,d_x, d_y      (\lambda/2) |dx - \nabla_x U|_2^2 + (\lambda/2)|dy -\nabla_yU|_2^2
%                 
% 
% which is solved via Bregman splitting.
% 
% Inputs
% ------
%   img : (m by n real) noisy image to be de-noised.
%   labda : (real) weight on the TV term.
%   tol :   (real) Iterations stop when ||img^k - img^{k-1}||_2 < tol
% 