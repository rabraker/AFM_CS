classdef MoveEntityStatic < MeasEntityStatic
    %MOVEENTITYSTATIC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        index = 0;
    end
    
    methods
        function obj = MoveEntityStatic(x_ref, y_ref, N)
            obj@MeasEntityStatic(x_ref, y_ref, N, MoveEntityStatic.index);
            
        end
    end
    
end

