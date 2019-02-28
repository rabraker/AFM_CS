% pixmat = pin_along_column(pixmat, x1, x2)
function pixmat = pin_along_column(pixmat, x1, x2)

  xs = 1:size(pixmat,2);
  for k=1:size(pixmat,1)
    y1 = pixmat(k, x1);
    y2 = pixmat(k, x2);
    
    m = (y2-y1)/(x2-x1);
    b = y1 - m*x1;
    line = m*xs + b;
    pixmat(k,:) = pixmat(k,:) - line;
  end
  

end
