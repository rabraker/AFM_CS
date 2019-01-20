function [ output ] = IDCTfun_2D( z,E, npix)
  
  Z_mat =PixelVectorToMatrix(z, [npix, npix]);
  output_mat = idct2(Z_mat);
  
  output = PixelMatrixToVector(output_mat);
  output = output(find(E>0.5));

    
%  output = idct2(z);
%  output = output(find(E>0.5));
    
    
end

