% [ABCD_vec, Ns_p1] = lti2ABCD_vec(sys)
% 
% Convert an lti system into an ABCD_vec suitible for use in the state-space
% based labview compensator.
% 
% 
function [ABCD_vec, Ns_p1] = lti2ABCD_vec(sys)
    sys = ss(sys);
    Ns = length(pole(sys));
    nu = size(sys.b, 2);
    no = size(sys.c, 1);
    if nu >1 || no >1
        error('lti2ABCD_vec requires a SISO system')
    end
    
    if Ns >0
        Ns_p1 = Ns + 1;
        [a, b, c, d] = ssdata(balreal(sys));
        ABCD = [a, b; c, d];
        % reshape into a column vector along rows.
        ABCD_vec = reshape(ABCD', [], 1)';
    else
        Ns_p1 = 0;
        ABCD_vec = [];
    end

end