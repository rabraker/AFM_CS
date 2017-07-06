classdef MeasEntityStatic
    %MEASENTITYSTATIC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x_ref_s;
        y_ref_s;
        index_s;
    end
    
    methods
        function obj = MeasEntityStatic(x_ref, y_ref, N, index)
            % All measurements have a DIFFERENT index.
            obj.x_ref_s = ones(N, 1)*x_ref;
            obj.y_ref_s = ones(N, 1)*y_ref;
            obj.index_s = ones(N, 1)*index;
        end
        
        function mergedTraj = asVector(obj)
            N = size(obj.x_ref_s, 1);
            mergedTraj = zeros(N*3, 1);
            for ik=[1:3:3*N; 1:1:N]
                i=ik(1);
                k = ik(2);
                mergedTraj(i) = obj.x_ref_s(k);
                mergedTraj(i+1) = obj.y_ref_s(k);
                mergedTraj(i+2) = obj.index_s(k);
           end
            
        end
    end
    
end

