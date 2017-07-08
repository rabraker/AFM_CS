classdef MasterTrajster
    %MASTERTRAJSTER Summary of this class goes here
    %   Detailed explanation goes here
    
    
    % ME and MVE should be classes of measurement entity and move entity.
    % That means that should be able to be initialized with 
    % scalars (x_ref, y_ref, N, index) and must have a method asVector()
    % which returns the inter-leaved column vector.
    
    
    
    properties
        XR;           % row vectory of x-points
        YR;           % row vector of y-points
        num_entities; % number of 'entities'
        MoveEntity;   %
        MeasEntity;   %
    end
    
    methods
        function obj = MasterTrajster(XR, YR, MoveEntity, MeasEntity)
           check_xryr(XR, YR);
            obj.XR = XR;
            obj.YR = YR;
            obj.num_entities = length(XR);
            obj.MoveEntity = MoveEntity;
            obj.MeasEntity = MeasEntity;
            
        end
        
        function masterTraj = asVector(self)
            % INIT: State 2
            ME_0 = self.MeasEntity(0, 0, 1);   
            masterTraj = ME_0.asVector();
            
            for i = 1:length(self.XR)
                % state 4
                MVE_i = self.MoveEntity(self.XR(i), self.YR(i));   
                
                masterTraj = [masterTraj; MVE_i.asVector()];
                % State 2: Define this separatly, instead of just holding at the same 
                % XR, YR so that in the future this same stuff can define, eg. a
                % mu-path.
                ME_i = self.MeasEntity(self.XR(i), self.YR(i), i+1);   
                masterTraj = [masterTraj; ME_i.asVector()];
            end
            
        end
    end
    
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                  LOCAL FUNCTIONS                                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function check_xryr(XR, YR)
    % We must have the same number of x and y CS points.
     if length(XR) ~= length(YR)
        S = sprintf('Expected len(XR) == len(XY) but recieved len(XR) = %d, len(YR) = %d',...
            length(XR), length(YR));
        error(S)
     end
        
end
