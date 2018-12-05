
function [x_max, row_idx, col_idx] = matrix_max(X)
% [x_max, row_idx, col_idx] = matrix_max(X)
% Finds the maximum value AND index in a matrix X. The builtin matlab command
% will find max of X via max(X, [], [1,2]), but will not return the index.

if length(size(X)) ~= 2
  error('Expected 2D array')
end

x_max = max(X, [], [1,2]);

[row_idx, col_idx] = find(X == x_max);
  

end