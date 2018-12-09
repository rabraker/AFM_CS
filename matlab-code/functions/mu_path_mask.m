% [pixelifsampled, XR, YR] = mu_path_mask(pix, mu_pix, sub_sample_frac,...
%     width, microns2volts)
%
% Inputs
% ------
%    pix:   (int) number of pixels in width and height of image. We assume a
%    square image.
%    mu_pix:  (int) length of mu-path in pixels
%    sub_sample_frac: (double) sub sampling fraction, e.g., 0.1 for 10%
%    width:           (double) Optional, default=1. Width in microns 
%                     of the image
%    microns2vols:    (double) Optional, default=1. Conversion from 
%                     microns 2 volts, typically 1/5
%    for our stage.
%
% Outputs
% -------
%   pixelifsampled:  pix by pix matrix of ones and zeros, where a 1 indicates
%   the pixel should be sampled.
%   XR, YR:  (double vectors) Vectors in volts of the reference values that the
%   stage should move to when sampling.


function [pixelifsampled, XR, YR] = mu_path_mask(pix, mu_pix, sub_sample_frac,...
    width, microns2volts)
  
  if nargin < 4
    width = 1;
    microns2volts = 1;
  end
  
  pix_per_micron = pix/width;
  pixelifsampled = zeros(pix, pix);
  XR = [];
  YR = [];
  swtch = rand(pix,1);
  for n=1:pix % down rows
    
    if  swtch(n) >0.5
      m = 1;
      while m < pix - mu_pix
        if rand(1,1) < sub_sample_frac/mu_pix
          pixelifsampled(n, m:m+mu_pix) = 1;
          XR = [XR; ( (m - 1) / pix_per_micron) * microns2volts];
          YR = [YR; ( (n - 1) / pix_per_micron) * microns2volts];
          m = m + mu_pix;
        else
          m = m+1;
        end
      end
    else
      m = pix; % reverse direction for odd ones.
      while m > mu_pix
        if rand(1,1) < sub_sample_frac/mu_pix
          pixelifsampled(n, m-mu_pix:m) = 1;
          
          XR = [XR; ( (m - mu_pix) / pix_per_micron) * microns2volts];
          YR = [YR; ( (n - 1) / pix_per_micron) * microns2volts];
          m = m - mu_pix;
        else
          m = m-1;
        end
      end
    end
    
  end

end

% for n=1:pix % down rows
%     m = 1;
%    while m < pix  % accros columns. pix - mu_pix so that the paths
%                            % dont hang off outside the 5-micron square. 
%        if rand(1,1) < sub_sample_frac/mu_pix
%            if m < pix-mu_pix
%               pixelifsampled(n, m:m+mu_pix) = 1;
%               XR = [XR; ( (m - 1) / pix_per_micron) * microns2volts];
%               YR = [YR; ( (n - 1) / pix_per_micron) * microns2volts];
%            else
%                pixelifsampled(n, end-mu_pix:end) = 1;
%                XR = [XR; ( (pix-mu_pix - 1) / pix_per_micron) * microns2volts];
%                YR = [YR; ( (n - 1) / pix_per_micron) * microns2volts];
%            end
%           m = m + mu_pix;
%        else
%            m = m+1;
%        end
%    end
%    m
% end
% 
% for n=1:pix % down rows
%     m = 1;
%    while m < pix - mu_pix  % accros columns. pix - mu_pix so that the paths
%                            % dont hang off outside the 5-micron square. 
%        if rand(1,1) < sub_sample_frac/mu_pix
%           pixelifsampled(n, m:m+mu_pix) = 1;
%           XR = [XR; ( (m - 1) / pix_per_micron) * microns2volts];
%           YR = [YR; ( (n - 1) / pix_per_micron) * microns2volts];
%           m = m + mu_pix;
%        else
%            m = m+1;
%        end
%    end
% end