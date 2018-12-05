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
  
  end
  
  methods 
    function self = CsSim(Img_original, pix_mask)
      self.Img_original = Img_original;
      self.pix_mask = pix_mask;
      self.Img_sub_sampled = pix_mask.*Img_original;
    end

    
    function solve_bp(self, recalc)
      if nargin <2
        recalc = false;
      end
      if ~recalc && ~isempty(self.Img_bp) && sum(self.Img_bp(:)) ~= 0
        warning(['BP solution already calculated. Pass recalc flag ' ...
                 'to recompute'])
      end
      
      [n m] = size(self.Img_sub_sampled);
      
      tic
      I_vector = PixelMatrixToVector(self.Img_sub_sampled);

      pix_mask_vec = PixelMatrixToVector(self.pix_mask);
      % y, set of measurements. have to remove all the spots we didn't sample.
      I_vector = I_vector(find(pix_mask_vec>0.5)); 

      A = @(x) IDCTfun(x,pix_mask_vec);
      At = @(x) DCTfun(x,pix_mask_vec);

      Ir_bp = idct(l1qc_logbarrier(At(I_vector), A, At, I_vector, 0.1));
      Ir_bp = real(Ir_bp);
      
      self.Img_bp = PixelVectorToMatrix(Ir_bp,[n m]);
      time_bp = toc;
      
      fprintf('BP Time: %f\n', time_bp);
      
    end
    
  end
  
  
  
end
