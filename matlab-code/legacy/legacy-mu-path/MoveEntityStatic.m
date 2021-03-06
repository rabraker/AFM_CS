classdef MoveEntityStatic < CsEntity
    %MOVEENTITYSTATIC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        index = 0;
    end
    properties
        xref_vec;
        yref_vec;
        index_vec;
    end
    
    methods
        function self = MoveEntityStatic(x_ref, y_ref, N)
            index = 0;
            % All measurements have a DIFFERENT index.
            self.xref_vec = ones(N, 1)*x_ref;
            self.yref_vec = ones(N, 1)*y_ref;
            self.index_vec = ones(N, 1)*index;
        end
        
    end
    
    methods(Static)
        function self=factory(N)
            self = @(x, y)MoveEntityStatic(x,y,N);
        end
    end

    
end

