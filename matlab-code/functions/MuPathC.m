classdef MuPathC
  
  properties
    xt;
    yt;
    met_idx;
  end
  
  methods
    function self = MuPathC(xt, yt, met_idx)
      self.xt = xt;
      self.yt = yt;
      self.met_idx = met_idx;
    end
  end
    
end