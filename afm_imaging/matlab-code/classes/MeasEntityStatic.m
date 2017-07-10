classdef MeasEntityStatic < CsEntity
    %MEASENTITYSTATIC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        xref_vec;
        yref_vec;
        index_vec;
    end
    
    methods
        function obj = MeasEntityStatic(x_ref, y_ref, N, index)
            % All measurements have a DIFFERENT index.
            obj.xref_vec = ones(N, 1)*x_ref;
            obj.yref_vec = ones(N, 1)*y_ref;
            obj.index_vec = ones(N, 1)*index;
        end
    end
    methods(Static)
        function self = factory(N)
            self = @(x, y, ind)MeasEntityStatic(x,y,N,ind);
        end
    end
    
end
