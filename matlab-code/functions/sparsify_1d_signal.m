
function x_sparse = sparsify_1d_signal(x, k, fwd_transform, inv_transform)
% Given a 1d signal x, this function will transform to some other basis and
% "sparsify" the signal by sorted the coefficients according to the magnitude,
% and setting to zero all but the k largest coefficients. 
%
% Arguments
% ------
%   x : (double vector) 
%     Time domain signal to sparsify
%   k : (int) 
%       How many coefficients to keep.
%   fwd_transform: (function handle)
%                   Function to transform the sparsity basis, e.g., @dct
%   inv_transform: (function handle)
%                  Function to transform back to the time domain, e.g., @idct
%
% returns
% ---------
% x_sparse : (double vector)
%            x in the time domain basis, with only k non-zero coefficients in
%            the sparsity basis.
%
% See Also
% ---------
%   sparsify_2d_signal

  k = floor(k); % ensure integer.
  
  x_trans = fwd_transform(x);
  [~,idx_descend] = sort(abs(x_trans), 'descend');
  z = abs(x_trans(idx_descend));
  z_thresh = z(k);
  
  x_trans(abs(x_trans)<z_thresh) = 0;
%   length(
  x_sparse = inv_transform(x_trans);
end
