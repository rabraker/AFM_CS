function [ output ] = DCTfun( z,E )

    output = dct(addzeros(z,E));

end
