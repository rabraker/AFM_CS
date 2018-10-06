classdef MasterTrajster

   % Defines the data structure for building up FIFO data for use in
   % labview. See the example `cs-trajectories-example.m` in the examples
   % folder. 
    
    properties
        XR;           % row vectory of x-points
        YR;           % row vector of y-points
        meta_cell;     % cell-array of meta data for the ith measurement,
                      % e.g., how many samples we should collect 
                      % the current measurement for. 
        num_entities; % number of 'entities', ie, total number of 
                      % measurements.
        MoveEntity;   % A handle to a subclass of CsEntity.
        MeasEntity;   % A handle to a subclass of CsEntity.
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
        
        function masterTraj = as_vector(self)
            % Combines the interleaved data from all N ME and MVE into a
            % single, long vector. This is the data we want to use in
            % LabView. 
            
            % INIT: 
            masterTraj = [];
            for iter = 1:length(self.XR)
                ith_meta = self.meta_cell{iter};
                % get the actual data.
                MVE_i = self.MoveEntity(self.XR(iter), self.YR(iter));   
                ME_i = self.MeasEntity(self.XR(iter), self.YR(iter),...
                    ith_meta, iter);   
                % Concatenate in the master vector. 
                masterTraj = [masterTraj; MVE_i.as_vector(); ME_i.as_vector()];
               
            end
        end

            
        function [measCell, moveCell] = as_matrix(self)
            moveCell = {};
            measCell = {};
            for iter = 1:length(self.XR)
                % State 1:
                ith_meta = self.meta_cell{iter};
                MVE_i = self.MoveEntity(self.XR(iter), self.YR(iter));   

                moveCell{iter} = MVE_i.as_matrix();
                
                % State 2: Define this separatly, instead of just holding
                % at the same XR, YR so that in the future this same stuff 
                % can define, eg. a mu-path.
                ME_i = self.MeasEntity(self.XR(iter), self.YR(iter),...
                    ith_meta, iter);   
                measCell{iter} = ME_i.as_matrix();
            end

        end
        
        function write_csv(self, fname)
            % Write the interleaved vector data into a file. This is the
            % link to LabView and is exactly the data we want to stuff into
            % the FIFO buffer.
            csvwrite(fname, self.as_vector());
        end
        function H = visualize_sampling(self)
            % Plots, in the time-domain, the constructed movements and
            % measurements. I think it would be better to move the actually
            % plotting into the classes themselves. This way, for example,
            % the plotting MeasEntitySatic data could result in a dot
            % instead of nothing, like it currently does (since all the
            % points are at the same place).
            [measCell, moveCell] = self.as_matrix();
            N = size(measCell, 2);
            
            % Setup the figure.
            H = figure(2);clf;
            hold on
            axx = gca();
            xlabel('x')
            ylabel('y')

            cme = 'k';
            
            % Each iteration is a new entity with two parts: the *move*
            % (MVE_i) and the *measurement* (ME_i).
            for iter=1:N
                % step 1: MOVE
                MVE_i = moveCell{iter};
                
                % step 2: measure
                ME_i = measCell{iter};
                x_meas = ME_i(1,:);
                y_meas = ME_i(2,:);

                plot(axx, x_meas, y_meas, cme, 'linewidth', 2, 'Marker', 'x')

            end
            
        end % end self.visualize_sampling()
 
        
        function H = visualize(self)
            % Plots, in the time-domain, the constructed movements and
            % measurements. I think it would be better to move the actually
            % plotting into the classes themselves. This way, for example,
            % the plotting MeasEntitySatic data could result in a dot
            % instead of nothing, like it currently does (since all the
            % points are at the same place).
            [measCell, moveCell] = self.as_matrix();
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
