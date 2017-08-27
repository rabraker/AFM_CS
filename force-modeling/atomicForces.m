classdef atomicForces
% Taken from 
% High Performance feedback for fast scanning atomic force microscopes,
% G. Schitter, P. Menold, H.F. Knapp, F. Allgower, A. Stemmer, Rev. Sci.
% Instruments, vol 72, no 8.
%
% k1 = fo * R
% where
% R = tip radius
% fo = (2/3)pi^2 eps rho1 rho2 sig^2
%
% rho1 = number density of tip  (what does this mean ?)
% rho2 = number density of sample
% eps = minimum energy of lennard-jones potential
% sig = spherical approximaton of the molecule diameter (tip or sample???)
    
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

