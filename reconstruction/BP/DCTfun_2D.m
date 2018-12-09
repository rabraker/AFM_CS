function [ output ] = DCTfun_2D( z,E, npix)

  Z_sparse_mat = PixelVectorToMatrix(addzeros(z,E), [npix, npix]);
  
  out_mat = dct2(Z_sparse_mat);
  
  output = PixelMatrixToVector(out_mat);
  
  
%   output = dct2(addzeros(z,E));

end
