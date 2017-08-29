function F_r = force2(p,ell, params)

r = ell + p;

% if r > params.zo
    F_r = params.fo*(-(params.sigma./r).^2 + (1/30)*(params.sigma./r).^8);
% else
%    F_r = params.go*(params.zo - r)^(2/3); 
% end



end