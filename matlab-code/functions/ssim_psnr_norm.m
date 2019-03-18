function [psn, ssm] = ssim_psnr_norm(Im_master, Im_compare, varargin)
  
  Im_master = Im_master - mean(Im_master(:));
  Im_compare = Im_compare - mean(Im_compare(:));
%   Im_master = Im_master - min(Im_master(:));
%   Im_compare = Im_compare - min(Im_compare(:));
  
  mx1 = max(abs(Im_master(:)));
  mx2 = max(abs(Im_compare(:)));
  mx = max(mx1, mx2);
  
  psn = psnr(Im_master, Im_compare, mx);
  ssm = ssim(Im_master, Im_compare, 'DynamicRange', mx);
  
end