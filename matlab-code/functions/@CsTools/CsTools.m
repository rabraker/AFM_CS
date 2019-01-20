classdef CsTools
  
  
  methods (Static)
    xp = l1qc_logbarrier(x0, A, At, b, epsilon, lbtol, mu, ...
                              cgtol, cgmaxiter, verbose);
    [xp, up, niter, cgtot_iter] = l1qc_newton(x0, u0, A, At, b, epsilon, tau, newtontol,...
      newtonmaxiter, cgtol, cgmaxiter, Tii, verbose)
    [x, res, iter] = cgsolve(A, b, tol, maxiter, verbose, x0);
    
    
    vpix = pixmat2vec(Mat);
    
    Mpix = pixvec2mat(vpix, nrows);
    
    [eta, LBRes] = l1qc(x0, b, pix_idx, opts);
    
    pixelifsampled = muPathMaskGen(mupathLength,n,m,samplingRatio,RepeatSamplingFlag);
    
    [Ir] = bregman_split_delxy(I,E,epsy, epsx, mu, lam,gam,maxiter,tol)
    
    [Ir] = bregman_split(I,E, mu, gam,maxiter,tol);
    
    function opts = l1qc_opts(varargin)
    % Will build an options struct for l1qc.
      p = inputParser();
      p.addParameter('epsilon', 0.1);
      p.addParameter('mu', 10);
      p.addParameter('cgtol', 1e-8);
      p.addParameter('cgmaxiter', 200);
      p.addParameter('warm_start_cg', 0);
      p.addParameter('lbtol', 1e-3);
      p.addParameter('newton_tol', 1e-3);
      p.addParameter('newton_max_iter', 50);
      p.addParameter('verbose', 2);
      
      p.parse(varargin{:});
      opts = p.Results();
      
    end
    
    function [ b ] = Afun_dct(eta, pix_idx)
    % Given the CS equation
    % b = E * M * eta
    % where E is the subsampling matrix and M is the idct.
      b = idct(eta);
      b = b(pix_idx);
    end

    function [ eta ] = Atfun_dct(b, pix_idx, len_eta)
    % Computes the adjoint of the CS equation
    % b = E * M * eta
    % where E is the subsampling matrix and M is the idct. That is
    %
    % eta = M'*E'*b
    %
    % N2 is the length eta.
      
      eta = zeros(len_eta, 1);
      eta(pix_idx) = b;
      eta = dct(eta);
      
    end

    function [ eta ] = Atfun_dct2(b, pix_idx, N)
    % Computes the adjoint of the CS equation
    % b = E * M * eta
    % where E is the subsampling matrix and M is the 2D-idct. That is
    %
    % eta = M'*E'*b
    %
    % N2 is the length eta.
    
      z_sparse = zeros(N,1);
      z_sparse(pix_idx) = b;
      Z_sparse_mat = CsTools.pixvec2mat(z_sparse, N);
      
      Eta_mat = dct2(Z_sparse_mat);
      
      eta = CsTools.pixmat2vec(Eta_mat);
    end
    
    function [b] = Afun_dct2(eta, pix_idx, N)
    % Computes the adjoint of the CS equation
    % b = E * M * eta
    % where E is the subsampling matrix and M is the 2D-idct. That is
    %
    
      Eta_mat = CsTools.pixvec2mat(eta, N);
      B_mat = idct2(Eta_mat);
      b = CsTools.pixmat2vec(B_mat);
      b = b(pix_idx);

    end

  end
  
end