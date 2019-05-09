classdef CsSim < handle
% cs_sim = CsSim(Img_original, pix_mask)
% This is a basic class to handle simulating CS reconstructions.
% Basically, I container to do that access the BP algo.
  properties
    Img_original
    pix_mask
    Img_sub_sampled;
    Img_bp
    Img_smp1d
    sample_frac;
    npix;
  end
  
  methods 
    function self = CsSim(Img_original, pix_mask)
      self.Img_original = Img_original;
      self.pix_mask = pix_mask;
      self.npix = size(pix_mask, 1);
      self.sample_frac = sum(pix_mask(:))/self.npix^2;
      self.Img_sub_sampled = pix_mask.*Img_original;
    end
    function bp_vec =  solve_bp_other(self, A, At)
      
      [n m] = size(self.Img_sub_sampled);
      
      tic
      I_vector = PixelMatrixToVector(self.Img_sub_sampled);

      pix_mask_vec = PixelMatrixToVector(self.pix_mask);
      % y, set of measurements. have to remove all the spots we didn't sample.
      I_vector = I_vector(pix_mask_vec>0.5); 

      bp_vec = l1qc_logbarrier(At(I_vector), A, At, I_vector, 0.01);
      bp_vec = real(bp_vec);
      
      time_bp = toc;
      
      fprintf('BP Time: %f\n', time_bp, opts);
    end
    function solve_bp(self, recalc, use_2d, opts)
    % Solve the Basis Pursuit problem in either 1d or 2d. If in 1D, use the mex
    % function. 
    % Options
    % -------
    % recalc : (true|false), default false. Do not optimize if self.Img_bp is
    %          non-empty and non-zero.
    % use_2d : (true|false), default false. If true, compute using 2D-dct.
      if nargin <2
        recalc = false;
      end
      if nargin <3
        use_2d = false;
      end
      if ~recalc && ~isempty(self.Img_bp) && sum(self.Img_bp(:)) ~= 0
        warning(['BP solution already calculated, so skipping optimization.',...
          ' Pass recalc flag to recompute.']);
        return;
      end
      
      [n m] = size(self.Img_sub_sampled);
      
      if ~exist('opts', 'var')
        opts = l1qc_dct_opts('verbose', 2, 'l1_tol', 0); %'epsilon', 0.01);
      end
      pix_idx = find(CsTools.pixmat2vec(self.pix_mask) > 0.5);
      % b, set of measurements. have to remove all the spots we didn't sample.
      b = CsTools.pixmat2vec(self.Img_sub_sampled);
      b = b(pix_idx);
      min_b = min(b);
      max_b = max(b);
      max_diff_b = max_b - min_b;
      b = b/max_diff_b;
      
      tic
      if use_2d
        % A = @(x) CsTools.Afun_dct2(x, pix_idx, n, m);
        % At = @(x) CsTools.Atfun_dct2(x, pix_idx, n, m);
        % x0 = At(b);
        % eta_vec = CsTools.l1qc_logbarrier(x0, A, At, b, opts);
        % self.Img_bp = idct2(CsTools.pixvec2mat(eta_vec, n))*max_diff_b;
        [x_est, LBRes] = l1qc_dct(n, m, b, pix_idx, opts);
        self.Img_bp = CsTools.pixvec2mat(x_est*max_diff_b, n);
      else
        [x_est, LBRes] = l1qc_dct(n*m, 1, b, pix_idx, opts);
        self.Img_bp = CsTools.pixvec2mat(x_est*max_diff_b, n);
      end
      
      time_bp = toc;
      
      fprintf('BP Time: %f\n', time_bp);
      
    end
    
    function solve_nesta(self, recalc, use_2d, opts)
    % Solve the Basis Pursuit problem in either 1d or 2d, using NESTA (my version) 
    % To set opts, use NESTA_opts.m
    % Options
    % -------
    % recalc : (true|false), default false. Do not optimize if self.Img_bp is
    %          non-empty and non-zero.
    % use_2d : (true|false), default false. If true, compute using 2D-dct.
      if nargin <2
        recalc = false;
      end
      if nargin <3
        use_2d = false;
      end
      if ~recalc && ~isempty(self.Img_bp) && sum(self.Img_bp(:)) ~= 0
        warning(['BP solution already calculated, so skipping optimization.',...
          ' Pass recalc flag to recompute.']);
        return;
      end
      
      [n, m] = size(self.Img_sub_sampled);
      
      pix_idx = find(CsTools.pixmat2vec(self.pix_mask) > 0.5);
      
      if use_2d
          M_fun = @(x) CsTools.pixmat2vec(dct2(CsTools.pixvec2mat(x, n)));
          Mt_fun = @(x) CsTools.pixmat2vec(idct2(CsTools.pixvec2mat(x, n)));
      else
          M_fun = @(x) dct(x);
          Mt_fun = @(x) idct(x);
      end
      E_fun = @(x) CsTools.E_fun1(x, pix_idx);
      Et_fun = @(x) CsTools.Et_fun1(x, pix_idx, n, m);
      
      % b, set of measurements. have to remove all the spots we didn't sample.
      b = CsTools.pixmat2vec(self.Img_sub_sampled);
      b = b(pix_idx);
      min_b = min(b);
      max_b = max(b);
      max_diff_b = max_b - min_b;
      b = b/max_diff_b;

      if ~exist('opts', 'var')
        opts = NESTA_opts('Verbose', 10, 'errFcn', @(x)norm(x),...
            'U', M_fun, 'Ut', Mt_fun);
      end
      delta = 1e-2;
      mu = 1e-2;
      
      tic
      [x_est] = NESTA_mine(E_fun, Et_fun, b, mu, delta, opts);

      self.Img_bp = CsTools.pixvec2mat(x_est*max_diff_b, n);

      time_nesta = toc;
      
      fprintf('NESTA Time: %f\n', time_nesta);
      
    end
  end
  
  
  
end

