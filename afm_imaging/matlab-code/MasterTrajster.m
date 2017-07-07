classdef MasterTrajster
    %MASTERTRAJSTER Summary of this class goes here
    %   Detailed explanation goes here
    
    
    % ME and MVE should be classes of measurement entity and move entity.
    % That means that should be able to be initialized with 
    % scalars (x_ref, y_ref, N, index) and must have a method asVector()
    % which returns the inter-leaved column vector.
    
    
    
    properties
        XR;
        YR;
        num_entities;
    end
    
    methods
        function obj = MasterTrajster(XR, YR, ME, MVE)
           check_xryr(XR, YR);
            obj.XR = XR;
            obj.YR = YR;
            obj.num_entities = length(XR);
            
        end
    end
    
    end


function check_xryr(XR, YR)
    % We must have the same number of x and y CS points.
     if length(XR) ~= length(YR)
        S = sprintf('Expected len(XR) == len(XY) but recieved len(XR) = %d, len(YR) = %d',...
            length(XR), length(YR));
        error(S)
     end
        
end
