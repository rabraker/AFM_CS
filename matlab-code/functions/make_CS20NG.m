function img_mat = make_CS20NG(x_start, y_start)
  
  % Make these options later
  npix = 512;
  hole_width = npix/20;
  pitch = npix/10;
  
  img_mat = ones(npix, npix);
  m_pitch = ceil(pitch);
  n_hole = ceil(hole_width);
  
  % First lets, create a prototype hole. Then we will insert it into the master
  % image periodically.
  sq = ones(n_hole,n_hole);
  % circle radius
  r = n_hole/2;
  % center
  xc = n_hole/2;
  yc = n_hole/2;
  
  
  for x_pt=1:n_hole
    for y_pt=1:n_hole
      
      rad_pt = sqrt( (x_pt - xc)^2 + (y_pt - yc)^2 );
      
      if rad_pt <= r
        sq(y_pt, x_pt) = 0; % make it black
      end
    end
  end
  
  x_lft = x_start;
  while x_lft+n_hole < npix +n_hole
    y_tp = y_start;
    while y_tp+n_hole <= npix + n_hole
      img_mat(y_tp:y_tp+n_hole-1, x_lft:x_lft+n_hole-1) = sq;
      y_tp = y_tp + m_pitch;
    end
    x_lft = x_lft + m_pitch;
    
  end
  img_mat = img_mat(1:npix, 1:npix);

end