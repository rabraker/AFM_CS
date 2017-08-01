clc
close all
data_root = 'C:\Users\arnold\Documents\labview\afm_imaging\data'
datmat = csvread(fullfile(data_root, 'raster_8-1-2017_v3.csv'));
% datmat = csvread('C:\Users\arnold\Documents\labview\afm_imaging\data\data-out.csv');
datmat = datmat(:,2:end)';
lo = min(min(datmat));
hi = max(max(datmat));
% cmp = contrast(datmat);

figure(1)
 imshow(datmat, [lo, hi]);


[ypix, xpix] = size(datmat);

% Fit a plane to the data and remove. 
xs = [1:1:xpix];
X = reshape(repmat(xs, ypix,1),[],1);
ys = [1:1:ypix]';
Y = reshape(repmat(ys, xpix,1),[],1);

PHI = [X, Y, X*0 + 1;];


Z = reshape(datmat, [],1);

coeffs = PHI\Z;


mx = coeffs(1)
my = coeffs(2)
b = coeffs(3)

z_fit = X*mx + Y*my + b;

datmat = reshape(Z-z_fit, ypix, xpix);

figure(2);
lo = min(min(datmat));
hi = max(max(datmat));
imshow(datmat, [lo, hi]);