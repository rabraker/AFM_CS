function [ v ] = PixelMatrixToVector( X )
% Convert the matrix X to a single vector of pixels, going along the rows of X.
% 
% In other words,
% v = [X(1,:), X(2,:), ...]'

[n, m] = size(X);

for row_idx = 1:n
    for col_idx = 1:m
        v(m*(row_idx-1) + col_idx, 1) = X(row_idx, col_idx);
    end
end



end

