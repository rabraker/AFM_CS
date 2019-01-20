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
      
      fprintf('BP Time: %f\n', time_bp);
    end
    
    function solve_bp(self, recalc, use_2d)
    % Solve the Basis Pursuit problem in either 1d or 2d. If in 1D, use the mex
    % function. 
      if nargin <2
        recalc = false;
      end
      if nargin <3
        use_2d = false;
      end
      if ~recalc && ~isempty(self.Img_bp) && sum(self.Img_bp(:)) ~= 0
        warning(['BP solution already calculated, so skipping optimization.',...
          'Pass recalc flag to recompute']);
        return;
      end
      
      [n m] = size(self.Img_sub_sampled);
      
      tic
      
      pix_idx = find(CsTools.pixmat2vec(self.pix_mask) > 0.5);
      b = CsTools.pixmat2vec(self.Img_sub_sampled);
      b = b(pix_idx);
      
      % y, set of measurements. have to remove all the spots we didn't sample.
      opts = CsTools.l1qc_opts();
      if use_2d
        A = @(x) CsTools.Afun_dct2(x, pix_idx, n);
        At = @(x) CsTools.Atfun_dct2(x, pix_idx, n);
        x0 = At(b);
        eta_vec = CsTools.l1qc_logbarrier(x0, A, At, b, opts);
        self.Img_bp = idct2(CsTools.pixvec2mat(eta_vec, n));
      else
        % A = @(x) CsTools.Afun_dct(x, pix_idx);
        At = @(b) CsTools.Atfun_dct(b, pix_idx, n*m);
        x0 = At(b);
        eta_vec = CsTools.l1qc(x0, b, pix_idx-1, opts);
        self.Img_bp = CsTools.pixvec2mat(idct(eta_vec), n);
      end
      
      time_bp = toc;
      
      fprintf('BP Time: %f\n', time_bp);
      
    end
    
  end
  
  
  
end
