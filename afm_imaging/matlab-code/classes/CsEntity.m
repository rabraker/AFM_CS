classdef (Abstract) CsEntity
    %CSENTITY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Abstract)
        xref_vec;
        yref_vec;
        index_vec;
    end
    
    methods 
%         asVector(obj)
        function mergedTraj = asVector(self)
            N = size(self.xref_vec, 1);
            mergedTraj = zeros(N*3, 1);
            
            % 3-row by N-col matrix
            mergedTraj = self.asMatrix();
            % This does the interleaving.
            mergedTraj = mergedTraj(:);
        end
        
        function matTraj = asMatrix(self)
            % a(:)' ensures we get a row vector. So this will get us a 
            % 3-row by N-col matrix
            matTraj = [self.xref_vec(:)';
                          self.yref_vec(:)';
                          self.index_vec(:)'];
        end
    end
    
end




% [1] This is a naive loop implementation
%             for ik=[1:3:3*N; 1:1:N]
%                 i=ik(1);
%                 k = ik(2);
%                 mergedTraj(i) = obj.xref_vec(k);
%                 mergedTraj(i+1) = obj.yref_vec(k);
%                 mergedTraj(i+2) = obj.index_vec(k);
%             end
