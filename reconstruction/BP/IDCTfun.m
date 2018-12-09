function [ output ] = IDCTfun( z,E )

    output = idct(z);
    output = output(find(E>0.5));

end

