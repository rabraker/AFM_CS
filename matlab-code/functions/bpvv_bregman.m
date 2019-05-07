function [Ir] = bpvv_bregman(I, pix_idx, epsilon, mu, lmda, gam, maxiter, sigma)
% [Ir] = BregmanSplitwithVerticalPenalty(I, pix_idx, weight, mu, lmda, gam, maxiter, sigma)
%
% BregmanSplitwithVerticalPenalty is used to solve the following
% problem:
% min |\Psi x|_1 + epsilon/2*|\nabla_y x|_1  s.t. ||\Phi x - y||_2 < \sigma
%
%
%  - The tuning parameters lmbda and gam seem critical to the convergence
% 
% Arguments
% ----------
% I - sampled image with size n*m
% pix_idx - indexes of the sampled data, y = I(pix_idx).
% lmda - tuning parameter
% gam  - tuning parameter
% maxiter - max # of iteration
% sigma - Iterations stop when ||y_k - y0||_2 / ||y0|| < sigma
% 
% output
% ------
% Ir - reconstrcuted image

% Yufan Luo 1/20/2018
% Arnold Braker 5-1-2019
    
    I_vec = CsTools.pixmat2vec(I); 

    [n m] = size(I);
    
    % u = I_vec.*E_vec;
    u = zeros(n^2, 1);
    u(pix_idx) = I_vec(pix_idx);
    
    w = zeros(n^2,1);
    b_w = zeros(n^2,1);
    d_y = zeros(n^2,1);
    b_y = zeros(n^2,1);

    y0 = I_vec(pix_idx);

    Dely = sparse(n^2,n^2);
    Dely = epsilon * spdiags([[ones(n*(n-1), 1); zeros(n,1)], -ones(n^2, 1)], [0, n], Dely);
  
    % A = mu * PHI^T*PHI + \lambda * Dely^T*Dely + \gamma*I
    phi_diag = zeros(n^2, 1);
    phi_diag(pix_idx) = 1;
    PHI = sparse(n^2, n^2);
    PHI = spdiags(phi_diag, 0, PHI);
    
    A = mu*PHI + lmda*(Dely'*Dely) + speye(n^2)*gam;

    diag_A = spdiags(A, 0);
    M = @(x) x./diag_A;
    u_current = u;
    u_next = u_current;
    
    y_k = y0;
	iter = 1;
    while true
 
        u_current = u_next;

        for i = 1:maxiter

            rhs = mu*PhiT_fun(y_k,pix_idx, n^2) + lmda*Dely'*(d_y-b_y)+gam*idct(w-b_w);

            % u_next = A^-1 * rhs;
            [u_next, ~, ~, cgit] = pcg(A, rhs, 1e-8, 200, M); %, L1, L1T);
%             fprintf('cg-iters: %d \n', cgit);
            d_y = shrink_y(Dely*u_current + b_y, 1/lmda);
            w = shrink_y(dct(u_next)+b_w, 1/gam);

            b_y = b_y + Dely*u_next - d_y;

            b_w = b_w + dct(u_next) - w;

        end

        % PHI * u
        Phi_u = u_next(pix_idx);

        y_k = y_k + y0 - Phi_u;
        
        Ir = CsTools.pixvec2mat(u_next, n, m);
        figure(3)
        imagesc(Ir); colormap('gray');
        
       nrm_cur = norm(Phi_u-y0) / norm(y0);
        fprintf('current norm: %f   outer_iteration: %d \n', nrm_cur, iter);
        iter = iter + 1;
        if nrm_cur < sigma 
            break
        end
    end
end

function x = PhiT_fun(y, pix_idx, n)
    x = zeros(n, 1);
    x(pix_idx) = y;
end
function x_shrunk = shrink_y(x, gam)

    x_shrunk = sign(x).* max(abs(x) - gam, 0);
    
end
