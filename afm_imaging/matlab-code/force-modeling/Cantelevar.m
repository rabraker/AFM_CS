classdef Cantelevar
    %CANTELEVAR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        wo;
        Q
        k;
        m;
        z1;
        force_func;
    end
    
    methods
        function self = Cantelevar(wo, Q, k, m, z1, force_func)
           self.wo = wo;
           self.Q = Q;
           self.k = k;
           self.m = m;
           self.z1 = z1;
           self.force_func = force_func;
        end
        
        
        function xdot = dyn(self, t, x, u)
           dx1 = x(2);
           
           dx2 = -(self.wo/self.Q)*x(2) - (self.k/self.m)*(x(1) - self.z1 + u)...
               + (1/self.m)*self.force_func.force(x(1));
            
           xdot = [dx1; dx2];
            
        end
    end
    
end

