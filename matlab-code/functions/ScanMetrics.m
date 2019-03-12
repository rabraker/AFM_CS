classdef ScanMetrics
  properties
    damage;
    quality;
    psnr;
    ssim;
    time;
    rate;
    coverage;
    type;
  end
  
  methods
    function self = ScanMetrics(varargin)
      p = inputParser();
      p.addParameter('damage', Inf);
      p.addParameter('quality', Inf);
      p.addParameter('psnr', 0);
      p.addParameter('ssim', 0);
      p.addParameter('rate', 0);
      p.addParameter('time', 0);
      p.addParameter('type', 'na');
      p.addParameter('coverage', 0);
      
      p.parse(varargin{:});
      
      self.damage = p.Results.damage;
      self.quality = p.Results.quality;
      self.psnr = p.Results.psnr;
      self.ssim = p.Results.ssim;
      self.rate = p.Results.rate;
      self.time = p.Results.time;
      self.coverage = p.Results.coverage;
      self.type = p.Results.type;
      
    end
  end

end