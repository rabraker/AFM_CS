function [ output ] = DCTfun_parametrizedIDX(v, pix_mask_idx)
% t_vec is a work vector so we avoid creating a new one every
% iteration. It should be all zeros and have length(v).
  % Not sure this will work in a function handle
  
t = zeros(length(v),1);
    t(pix_mask_idx)=v;
  
    output = dct2(t);
     

end
