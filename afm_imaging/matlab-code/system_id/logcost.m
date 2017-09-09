
function J = logcost(theta, P, w)

    w_weight = [w(1);w;w(end)];
    
   logW = log(w_weight(3:end)) - log(w_weight(1:end-2));
%     logW = log(W);

   
    H_frf = H(w, theta);
    J_vec = abs(log(H_frf) - log(P)).^2;
    J_vec = J_vec.*logW;
%     J_vec = J_vec
    J = sum(J_vec);



end


function H_frf = H(w, theta)

g = tf(theta(1:3), [1 theta(4:5)], 40e-6);
H_frf = squeeze(freqresp(g, w));

end