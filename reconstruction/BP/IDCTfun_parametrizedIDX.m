function [ output ] = IDCTfun_parametrizedIDX( z, pix_mask_IDX )
  if any( imag(z) ~= 0)
    keyboard
  end
    output = idct2(z);
    output = output(pix_mask_IDX);

end

