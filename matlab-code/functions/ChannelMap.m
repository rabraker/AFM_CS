classdef ChannelMap
  properties
    names;
    x;
    y;
    ze;
    uz;
    met;
  end
    methods
      function self = ChannelMap(idx_s)
        self.names = {'x', 'y', 'ze', 'uz', 'met'};
        for k=1:length(self.names)
          self.(self.names{k}) = idx_s(k);
        end
      end
    end
      
  
end
  
  