%-----------------------------------------------------------------------
%  SMP with 1D DCT basis -----------------------------------------------\


%  inputs
%  I_sample        - sampled image (the size should be n^2)
%  pixelifsampled  - knowledgement matrix (the size should be n^2, if pixel i,j is 
%                               known pixelifsampled(i,j) = 1, else pixelifsampled(i,j)=0)
%  maxiter         - max value of # iter

%  output
%  Ir  -  reconstruction result (with size n^2)



function [ Ir ] = SMP_1D( I_sample,pixelifsampled,maxiter)


pixelifsampled = PixelMatrixToVector(pixelifsampled);
I_sample = PixelMatrixToVector(I_sample).*pixelifsampled;


N = length(pixelifsampled);
n = sqrt(N);
m = length(find(pixelifsampled>0.5));


J_hat = dct2(I_sample)*N/m;

J_recover = zeros(N,1);

setCover = [];

k = 1;


while(k<=maxiter)
    
    
    [a b] = max(abs(J_hat));
    
    
    vb = zeros(N,1);
    
    for i = 1:N
        vb(i) = 1;
        if pixelifsampled(i) == 0
            vb(i) = 0;
        elseif (b == 1)
            vb(i) = vb(i)*sqrt(1/N)*cos(pi*(2*i-1)*(b-1)/2/N);
        else
            vb(i) = vb(i)*sqrt(2/N)*cos(pi*(2*i-1)*(b-1)/2/N);
        end
    end

    T = dct2(vb)*(N/m);
    
    J_recover(b) = J_recover(b) + J_hat(b)/T(b);
    if isempty(find(b == setCover))
        setCover(end+1)=b;
    end
    
    J_hat = J_hat - T*J_recover(b);
    J_hat(setCover)=0;
    
    if mod(k, 200) == 0
      fprintf('SMP-1D iter:  %d\n', k);
    end
    k = k+1;
    

end

    Ir = idct2(J_recover);
    Ir = PixelVectorToMatrix(Ir, [n n]);

end


