classdef TestMuPathTraj < matlab.unittest.TestCase
  
  properties
%     mpt_inst;
    width_mic;
    mu_len;
    Ts;
    x_rate;
    npix;
    mu_pix;
    pix_mask;
  end
  
  methods
    
    function self = TestMuPathTraj()
      self.width_mic = 10;
      self.Ts = 0.1;
      
      self.x_rate = 1;
      
      npix = 16;
      mu_pix = 4;
      
      self.mu_len = (self.width_mic/npix)*mu_pix;
      pix_mask = zeros(npix,npix);
      idx = 0:mu_pix-1;
      pix_mask(1,1+idx) = 1;
      pix_mask(1, 10+idx) = 1;
      pix_mask(11, 6+idx) = 1;
      self.pix_mask = pix_mask;
   
    end
  end
  
  methods (Test)
  
    function test_adjust_pre_pad(self)
      mpt = MuPathTraj(self.pix_mask, self.width_mic, self.mu_len, self.x_rate, self.Ts);
      mpt.x_rate_mic_per_sec = 5; % adjust so we get 1 volt per sec.
      mpt.pre_pad_samples = 10;
      xr_ = 0;
      yr_ = 0;
      N_ = 0;
      [xr, yr, N] = mpt.adjust_pre_pad(xr_, yr_, 10);
      
      self.assertEqual(yr, yr_);
      self.assertEqual(N, 20);
      self.assertEqual(xr, xr_ - 1)
      
    end
    
    function test_connect_mu_paths(self)
      
      pix_min = 3; % Threshold for connecting two paths.
      
      idx = 0:3;
      r0 = zeros(1,16);  
      xr_idx = [3, 6];
      
      % connect two
      r1 = r0;
      r1(3+idx) = 1;
      r1(6+idx) = 1;
      r1_c = r0; 
      r1_c(3:6+4-1) = 1;
      
      % connect first two, but not last.
      r2 = r0;
      r2(1+idx) = 1;
      r2(4+idx) = 1;
      r2(11+idx) = 1;
      r2_c = r0; 
      r2_c(1:8-1) = 1;
      r2_c(11:11+4-1) = 1;
      
      % Check when only one exists:
      r3 = r0;
      r3(4+idx) = 1;
      r3_c = r3; 

      % connect all three
      r4 = r0;
      r4(1+idx) = 1;
      r4(4+idx) = 1;
      r4(8+idx) = 1;
      r4_c = r0; 
      r4_c(1:12-1) = 1;
      
      % Check when we dont connect the first one, but the next two
      % connect all three
      r5 = r0;
      r5(1+idx) = 1;
      r5(6+idx) = 1;
      r5(9+idx) = 1;
      r5_c = r0; 
      r5_c(1+idx) = 1;
      r5_c(6:9+4-1) = 1;
      
      pix_mask = [r1;r2;r3;r4;r5;zeros(11,16)];
      pix_mask_exp = [r1_c;r2_c;r3_c;r4_c;r5_c;zeros(11,16)];
      
      mpt = MuPathTraj(pix_mask, self.width_mic, self.mu_len, self.x_rate, self.Ts);
      
      mpt.connect_mu_paths(pix_min/mpt.x_rate_pix_per_sec);
      pix_mask_new = mpt.pix_mask;
      
      self.assertEqual(pix_mask_new, pix_mask_exp);
      
    end
    function test_connect_row(self)
      mpt = MuPathTraj(self.pix_mask, self.width_mic, self.mu_len, self.x_rate, self.Ts);
      
      xr_idx = [3, 6];
      mu_pix_s = [4,4];
      pix_min = 3;
      
      [idx_con, mu_pix_con] = mpt.connect_row(xr_idx, mu_pix_s, pix_min);
      self.assertEqual(idx_con, 3)
      self.assertEqual(mu_pix_con, 8)

      % connect first two, but not last.
      xr_idx = [1, 4, 11];
      mu_pix_s = [4, 4, 4];
      
      [idx_con, mu_pix_con] = mpt.connect_row(xr_idx, mu_pix_s, pix_min);
      self.assertEqual(idx_con, [1;11])
      self.assertEqual(mu_pix_con, [8; 4])

      % Check when only one exists:
      xr_idx = [4];
      mu_pix_s = [3];
      
      [idx_con, mu_pix_con] = mpt.connect_row(xr_idx, mu_pix_s, pix_min);
      self.assertEqual(idx_con, [4])
      self.assertEqual(mu_pix_con, [3])

      % connect all three
      xr_idx = [1, 4, 8];
      mu_pix_s = [4, 4, 4];
      
      [idx_con, mu_pix_con] = mpt.connect_row(xr_idx, mu_pix_s, pix_min);
      self.assertEqual(idx_con, [1])
      self.assertEqual(mu_pix_con, [12])
      
      % Check when we dont connect the first one, but the next two
      % connect all three
      xr_idx = [1, 6, 9];
      mu_pix_s = [4, 4, 4];
      
      [idx_con, mu_pix_con] = mpt.connect_row(xr_idx, mu_pix_s, pix_min);
      self.assertEqual(idx_con, [1, 6]')
      self.assertEqual(mu_pix_con, [4, 8]')

    end
    
    
    function test_build_xr_yr_starts(self)
      
      mpt = MuPathTraj(self.pix_mask, self.width_mic, self.mu_len, self.x_rate, self.Ts);
      
      mpt.build_xr_yr_pix_starts();
      
      self.assertEqual(mpt.XR_pix_starts, [1,10, 6]')
      self.assertEqual(mpt.YR_pix_starts, [1,1, 11]')
      
    end
    
  end
  
  
  
end