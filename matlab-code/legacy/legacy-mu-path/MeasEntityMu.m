classdef MeasEntityMu < CsEntity
    %MEASENTITYMU Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        xref_vec;
        yref_vec;
        index_vec;
        x_rate; % volts per sample
        y_rate;
    end
    
    methods
        function self = MeasEntityMu(x0, y0, N, index, xy_rate)
            self.x_rate = xy_rate(1);
            self.y_rate = xy_rate(2);
            
            self.xref_vec = MeasEntityMu.gen_ramp(x0, self.x_rate, N);
            self.yref_vec = MeasEntityMu.gen_ramp(y0, self.y_rate, N);
            self.index_vec = ones(N, 1)*index;
            self.index_vec(end) = -1;
        end
        
    end
    
    methods(Static)
        function self = factory(xy_rate)
           self = @(x0, y0, N, ind)MeasEntityMu(x0, y0, N, ind, xy_rate); 
        end

        function x_vec = gen_ramp(x0, x_rate, N)
            % Given a starting point x0, and a ramp rate x_rate and number
            % of samples N, generate a reference trajectory vector of N
            % points. 
            
            % N-1 because x0 is sample 1.
            x_N = x0 + (N-1)*x_rate;
            
            % This will put out a zero length vector if x_rate =0.
            % x_vec = [x0:x_rate:x_N];
            
            x_vec = linspace(x0, x_N, N);
            
        end
    end
end

