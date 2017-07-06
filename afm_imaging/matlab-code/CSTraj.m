classdef CSTraj
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        moveIndex = 0;
    end
    
    methods(Static)
        function [x_ref_s, y_ref_s, index_s] = move_entity_static(x_ref, y_ref, N)
            % All moves have a the SAME index. 
            [x_ref_s, y_ref_s, index_s] = CSTraj.meas_entity_static(x_ref, y_ref, N, CSTraj.moveIndex);
        end
        function [x_ref_s, y_ref_s, index_s] = meas_entity_static(x_ref, y_ref, N, index)
            % All measurements have a DIFFERENT index.
            
            x_ref_s = ones(N, 1)*x_ref;
            y_ref_s = ones(N, 1)*y_ref;
            index_s = ones(N, 1)*index;
        end
        
        function mergedTraj = entity2vector(x_ref_s, y_ref_s, index_s)
           N = size(x_ref_s, 1);
            mergedTraj = zeros(N*3, 1);
            
            for ik=[1:3:3*N; 1:1:N]
                i=ik(1);
                k = ik(2);
                mergedTraj(i) = x_ref_s(k);
                mergedTraj(i+1) = y_ref_s(k);
                mergedTraj(i+2) = index_s(k);
           end
        end
    end
    
end

