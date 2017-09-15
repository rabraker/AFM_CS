function bp_im = bp_cs_recon(I_sampled, pixelifsampled)

    [n m] = size(I_sampled);

    I_vector = PixelMatrixToVector(I_sampled);

    pixelifsampled_vec = PixelMatrixToVector(pixelifsampled);
    I_vector = I_vector(find(pixelifsampled_vec>0.5));

    A = @(x) IDCTfun(x,pixelifsampled_vec);
    At = @(x) DCTfun(x,pixelifsampled_vec);

    Ir_bp = idct(l1qc_logbarrier(At(I_vector), A, At, I_vector, 0.1));
    Ir_bp = real(Ir_bp);
    bp_im = PixelVectorToMatrix(Ir_bp,[n m]);


end