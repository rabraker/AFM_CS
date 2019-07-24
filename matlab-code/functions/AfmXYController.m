classdef AfmXYController
   
    
    properties
        xdir_cntrl;
        ydir_cntrl;
    end
    methods
        function self = AfmController()
      p = inputParser();

      p.addParameter('r', 0); % These must not be empty, else simulink fails.
      p.addParameter('w', 0);
      p.addParameter('rp', 0);
      p.addParameter('wp', 0);
      p.addParameter('d', 0);
      p.addParameter('ws', 0);
      p.addParameter('dp', fi_zero);
      p.addParameter('wsp', fi_zero);
      p.addParameter('gdrift',  kd_sys);
      p.addParameter('gdrift_inv', kd_sys);
      p.addParameter('Dx_ff',  kd_sys);
      p.addParameter('Dy_ff',  kd_sys);
      
      p.addParameter('nw', 0);
      p.addParameter('nf', 0);
      p.addParameter('useNbar', false)
      p.parse(varargin{:});
      
      self.thenoise = p.Results.thenoise;
      self.step_amp =  p.Results.step_amp;
      self.r = p.Results.r;
      self.w = p.Results.w;
      self.rp = p.Results.rp;
      self.wp = p.Results.wp;
      self.d = p.Results.d;
      self.ws = p.Results.ws;
      self.dp = p.Results.dp;
      self.wsp = p.Results.wsp;
      self.useNbar = p.Results.useNbar;
            
        end
    end
    
end