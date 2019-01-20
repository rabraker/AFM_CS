
function [x_min, row_idx, col_idx] = matrix_min(X)
% [x_max, row_idx, col_idx] = matrix_min(X)  
% Finds the manimum value AND index in a matrix X. The builtin matlab command
% will find min of X via min(X, [], [1,2]), but will not return the index.

if length(size(X)) ~= 2
  error('Expected 2D array')
end

x_min = min(X, [], [1,2]);

[row_idx, col_idx] = find(X == x_min);
  

end