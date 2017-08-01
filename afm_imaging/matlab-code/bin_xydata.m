% Because we are sampling faster than we are moving, we typically end up
% with many data points per pixel. We need a way to average all of these
% together into a single point per pixel. This is complicated by the fact
% that for many pixels we have no data, and the data we read in (in CS) is
% not necessarily in order. 
%
% To construct an image, we need to average the HEIGHT data based on
% binning from the x-y data. For example 
%
%
% [xy     |      |     xy  ]
% [    xy |      |  xy     ]
% [  xy   |      |    xy   ]
% [-------------------------
% [       |      |         ]
% [       |      |  xy     ]
% [       |      |xy       ]
% [------------------------]
% [       | xy   |         ]
% [       |   xy |         ]
% [       |      |         ]

format short
e = 0:1:3;
lenMu = 10;
npix = 32;
edge_x = linspace(0, lenMu, npix);
edge_y = linspace(0, lenMu, npix);

x = [.5, .76, 2, 2.5, 3]
y = [.1, .2, 2.1, 2.2, .1]
z = [10, 5, 3,4, 20];



[~, ~,~,idx, idy] = histcounts2(x_mu, y_mu, edge_x, edge_y);

for i = 1:max(idx)
    for j=1:max(idy)
        a = find(idx == i);
        b = find(idy==j);
        z_index = intersect(a, b);
        if ~isempty(z_index)
            topo(j, i) = mean(z(z_index));
        end
    end
end


