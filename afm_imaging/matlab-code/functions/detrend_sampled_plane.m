% [ I_fit ] = detrend_sampled_plane(I, pixelifsampled)
%
% Fits a plane to a sub-sampled data set and removes it.
% inputs
% ------
% I : pix x pix image data
% pixelifsampled : pix x pix binary matrix. entry i,j == 1 if the data
% there was sampled, else (i,j) == 0.
%
%
function [ I_fit ] = detrend_sampled_plane(I, pixelifsampled)


% Fit a plane to the data
% 
% Z = [x1, y1, 1]
%     |x2, y2, 1|[ mx ]
%     |x3, y3, 1|| my |
%     |x4, y4, 1|[ b  ]
%     |:    :  :|

% Z = zeros(length(find(pixelifsampled == 1)), 1);
PHI = [];
Z = [];
for n_row = 1:size(I, 1)
    for m_col = 1:size(I, 2)
        if pixelifsampled(n_row, m_col) == 1
            Z = [Z; I(n_row, m_col)];
            PHI  = [PHI; n_row, m_col, 1];
%             PHI  = [PHI; m_col, n_row, 1];
        end
    end
end

mxmyb = PHI\Z;
mx = mxmyb(1);
my = mxmyb(2);
b = mxmyb(3);

z_fit = PHI(:,1)*mx + PHI(:,2)*my + PHI(:,3)*b;
Z = Z -z_fit;

% Now, rebuild the image

I_fit = I*0;

for k=1:length(Z)
    x = PHI(k, 1);
    y = PHI(k, 2);
    z = Z(k);
    I_fit(x, y) = z;
end

end

