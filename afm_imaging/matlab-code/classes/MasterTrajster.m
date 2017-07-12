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
        meta_cell;     % cell-array of meta data for the ith measurement,
                      % e.g., how many samples we should collect 
                      % the current measurement for. 
        num_entities; % number of 'entities', ie, total number of 
                      % measurements.
        MoveEntity;   %
        MeasEntity;   %
    end
    
    methods
        function obj = MasterTrajster(XR, YR, meta_cell,...
                MoveEntity, MeasEntity)
           % Idiot check everything is the same length.
           check_xryr(XR, YR);
           check_xryr(XR, meta_cell);
           
            obj.XR = XR;
            obj.YR = YR;
            obj.meta_cell = meta_cell;
            
            obj.num_entities = length(XR);
            obj.MoveEntity = MoveEntity;
            obj.MeasEntity = MeasEntity;
        end
        
        function masterTraj = asVector(self)
            % INIT: State 2
            masterTraj = [];
            for iter = 1:length(self.XR)
                % state 4
                ith_meta = self.meta_cell{iter};
                MVE_i = self.MoveEntity(self.XR(iter), self.YR(iter));   
                
                masterTraj = [masterTraj; MVE_i.asVector()];
                % State 2: Define this separatly, instead of just holding
                % at the same XR, YR so that in the future this same stuff 
                % can define, eg. a mu-path.
                ME_i = self.MeasEntity(self.XR(iter), self.YR(iter),...
                    ith_meta, iter);   
                masterTraj = [masterTraj; ME_i.asVector()];
            end
        end

            
        function [measCell, moveCell] = asMatrix(self)
            moveCell = {};
            measCell = {};
            for iter = 1:length(self.XR)
                % state 1
                ith_meta = self.meta_cell{iter};
                MVE_i = self.MoveEntity(self.XR(iter), self.YR(iter));   

                moveCell{iter} = MVE_i.asMatrix();
                
                % State 2: Define this separatly, instead of just holding
                % at the same XR, YR so that in the future this same stuff 
                % can define, eg. a mu-path.
                ME_i = self.MeasEntity(self.XR(iter), self.YR(iter),...
                    ith_meta, iter);   
                measCell{iter} = ME_i.asMatrix();
            end

        end
        
        function H = visualize(self)
            [measCell, moveCell] = self.asMatrix();
            N = size(measCell, 2);
            
            % Setup the figure.
            H = figure(1);clf;
            subplot(1,2,1); hold on;
            axx = gca();
            title('x')
            subplot(1,2,2); hold on;
            axy = gca();
            title('y')
            
            cmve = 'r';
            cme = 'k';
            n_ = 0;
            
            % Each iteration is a new entity with two parts: the *move*
            % (MVE_i) and the *measurement* (ME_i).
            for iter=1:N
                % step 1: MOVE
                MVE_i = moveCell{iter};
                x_move = MVE_i(1,:);
                y_move = MVE_i(2,:);
                n_move_i = n_ + length(x_move);
                % The move time vector begins at the end of the previous
                % iterations measure time vector. 
                t_move_i = [n_:1:n_move_i-1];
                
                % step 2: measure
                ME_i = measCell{iter};
                x_meas = ME_i(1,:);
                y_meas = ME_i(2,:);
                % The measure time vector begins at the *end* of the move
                % time vector. 
                n_meas_i = n_move_i + length(x_meas);
                t_meas_i = [n_move_i:1:n_meas_i-1];
                
                
                plot(axx, t_move_i, x_move, cmve, 'linewidth', 2)
                plot(axx, t_meas_i, x_meas, cme, 'linewidth', 2)
                
                plot(axy, t_move_i, y_move, cmve, 'linewidth', 2)
                plot(axy, t_meas_i, y_meas, cme, 'linewidth', 2)
                % Increment the beginning of time.
                n_ = t_meas_i(end);

            end
            
        end % end self.visualize()
            
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
