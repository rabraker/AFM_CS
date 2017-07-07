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
        
        function mergedTraj = asVector(obj)
            N = size(obj.xref_vec, 1);
            mergedTraj = zeros(N*3, 1);
            for ik=[1:3:3*N; 1:1:N]
                i=ik(1);
                k = ik(2);
                mergedTraj(i) = obj.xref_vec(k);
                mergedTraj(i+1) = obj.yref_vec(k);
                mergedTraj(i+2) = obj.index_vec(k);
           end
            
        end
    end
    
end

