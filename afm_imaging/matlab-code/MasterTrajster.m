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
        CSP;
    end
    
    methods
        function obj = MasterTrajster(XR, YR, CSP)
           check_xryr(XR, YR);
            obj.XR = XR;
            obj.YR = YR;
            obj.num_entities = length(XR);
            obj.CSP = CSP;
        end
        
        function masterTraj = asVector(self)
            % INIT: State 2
            N = 2500;
            ME_0 = self.CSP.me_factory(0, 0, N, 1);   
            masterTraj = ME_0.asVector();
            
            for i = 1:length(self.XR)
                % state 4
                MVE_i = self.CSP.mve_factory(self.XR(i), self.YR(i), N);   
                
                masterTraj = [masterTraj; MVE_i.asVector()];
                % State 2: Define this separatly, instead of just holding at the same 
                % XR, YR so that in the future this same stuff can define, eg. a
                % mu-path.
                ME_i = self.CSP.me_factory(self.XR(i), self.YR(i), N, i+1);   
                masterTraj = [masterTraj; ME_i.asVector()];
            end
            
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
