function t = addzeros_idx(v, pix_mask_idx)

    t = zeros(length(v),1);

    t(pix_mask_idx)=v;


end

