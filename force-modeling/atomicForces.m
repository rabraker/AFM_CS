classdef atomicForces
    %ATOMICFORCES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sigma;
        k1;
        
    end
    
    methods
        function obj = atomicForces(sigma, k1)
            obj.sigma = sigma;
            obj.k1 = k1;
        end
        
        function F_r = force(obj, r)
           F_r = obj.k1*(-(obj.sigma./r).^2 + (1/30)*(obj.sigma./r).^8);
            
        end
    end
    
end

