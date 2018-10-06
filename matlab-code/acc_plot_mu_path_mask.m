% The middle of our grating has a 500 nm pitch size. The holes themselves 
% appear to be about 1/4 of this. Thus, i think that if we travel for
% 1/2 * 500 nm = 250 nm, we should be able to detect wheather or not we 
% started in a hole. The idea being that if no jump occures in the u_z 
% signal, we assume we are on the flat, and if we see a jump, we can 
% backout which part is the hole. 
% 
% This will essentially come down to edge detection. If an edge is present,
% then we know the lower part of the signal was in a hole. If no edge is
% present, the entire signal was on the flat part of the grating. 
clear, clc
close all
addpath('classes')
addpath('functions')

%-------------- Location of System Model -------------------------
dataroot      = fullfile(getMatPath(), 'AFM_SS',...
                'System_Identification', 'data', 'data_xaxis'); 
expName           = ['22-Jun-2016_exp01'];
modFitName    = [expName, '.mat'];
modFitPath    = fullfile(dataroot, modFitName);
load(modFitPath, 'modelFit')


% FitNum    = 'FitNum001';
PLANT_init_x = ltiFit(modFitPath, 'SS02').sys;
%-------------- Define scan Params -------------------------

Ts = 40e-6;
TsTicks = 1600;
Ki_x = 0.01;

if 1
    width = 20;  % microns
    pix =64;  % image resolution.
    mu_length = 1;  % 1000 nm. length of the horizontal mu-path. 
    sub_sample_frac = 0.15;  % Percent of pixels to subsample. 
end
% Unit conversions.
pix_per_micron = pix/width;
mu_pix = ceil(mu_length*pix_per_micron);

% ************************************************
% ************************************************

Ts = 40e-6;  % AFM sample rate. 
Fs = 1/Ts;
microns2volts = 1/5;  %10/50;
raster_freq = 1;  % hz
raster_period = 1/raster_freq;

microns_per_second = width/(raster_period/2);
pixels_per_second = pix_per_micron * microns_per_second;
volts_per_second = microns_per_second*microns2volts;
volts_per_sample = volts_per_second * Ts;  % Ts = seconds per sample

% Convert mu-path length in pixels to volts
mu_micron = (1/pix_per_micron) * mu_pix;
mu_volts = mu_micron * microns2volts;
mu_Nsamples = ceil(mu_volts / volts_per_sample);

ME = MeasEntityMu.factory([volts_per_sample, 0]);

N_mve = 1;
XR = [];
YR = [];
% ************************************************


pixelifsampled = zeros(pix, pix);
% 
% for n=1:pix % down rows
%     m = 1;
%    while m < pix - mu_pix  % accros columns. pix - mu_pix so that the paths
%                            % dont hang off outside the 5-micron square. 
%        if rand(1,1) < sub_sample_frac/mu_pix
%           pixelifsampled(n, m:m+mu_pix) = 1;
%           XR = [XR; ( (m - 1) / pix_per_micron) * microns2volts];
%           YR = [YR; ( (n - 1) / pix_per_micron) * microns2volts];
%           m = m + mu_pix;
%        else
%            m = m+1;
%        end
%    end
% end
I = ones(pix,pix)-pixelifsampled;
I = pixelifsampled;
I = flipud(I); % why the fuck does it start from the top???
% imshow(I)
F3=figure(210); clf
figwidth = 2.4;
figheight = 2.4;
set(F3, 'Units', 'Inches', 'Position', [-10,3, figwidth, figheight],...
    'PaperUnits', 'Inches', 'PaperSize', [figwidth, figheight])
set(F3, 'Color', 'w');

xpad = .05;
wd = .9
lft1 = .1; % start of first col;
ht_im = .89; % image heights
bt_im = 0.11;
ax1 = axes('Position', [lft1, bt_im, wd, ht_im]);
hold on
xlim([0,pix])
ylim([0,pix])
swtch = rand(pix,1);
for n=1:pix % down rows
    
%     if  swtch(n) >0.5;
        m = 1;
       while m < pix - mu_pix
           if rand(1,1) < sub_sample_frac/mu_pix
              pixelifsampled(n, m:m+mu_pix) = 1;
              
              xs = m:m+mu_pix;
              ys = ones(1, length(xs))*n;
              plot(xs, ys, '-k', 'LineWidth', 2)
              XR = [XR; ( (m - 1) / pix_per_micron) * microns2volts];
              YR = [YR; ( (n - 1) / pix_per_micron) * microns2volts];
              m = m + mu_pix;
           else
               m = m+1;
           end
       end
%     else
%        m = pix; % reverse direction for odd ones. 
%        while m > mu_pix
%            if rand(1,1) < sub_sample_frac/mu_pix
%               pixelifsampled(n, m-mu_pix:m) = 1;
%               xs = m:m+mu_pix;
%               ys = ones(1, length(xs))*n;
%               plot(xs, ys, '-k', 'LineWidth', 2)
%               
%               XR = [XR; ( (m - mu_pix) / pix_per_micron) * microns2volts];
%               YR = [YR; ( (n - 1) / pix_per_micron) * microns2volts];
%               m = m - mu_pix;
%               
%            else
%                m = m-1;
%            end
%        end
%     end

end

ax1.Box = 'on'

actual_sub_sample_frac = length(find(pixelifsampled == 1))/pix^2;

fprintf('Desired sub sample fraction: %f\n', sub_sample_frac)
fprintf('Actual  sub sample fraction: %f\n', actual_sub_sample_frac)

fig_path = fullfile(getfigroot, 'mu_path_mask.eps');
saveEps(F3, fig_path)

%%
close all
I = ones(pix,pix)-pixelifsampled;
I = pixelifsampled;
I = flipud(I); % why the fuck does it start from the top???
% imshow(I)
F3=figure(210); clf
figwidth = 2.4;
figheight = 2.4;
set(F3, 'Units', 'Inches', 'Position', [-10,3, figwidth, figheight],...
    'PaperUnits', 'Inches', 'PaperSize', [figwidth, figheight])
set(F3, 'Color', 'w');

xpad = .05;
wd = .9
lft1 = .1; % start of first col;
ht_im = .89; % image heights
bt_im = 0.11;
ax1 = axes('Position', [lft1, bt_im, wd, ht_im]);


%%
% hold on
% for i_row = 1:size(I,1)
%     the_row = I(i_row, :)
%     xs_inds = find(the_row ==1);
%     plot(xs_inds, the_row(xs_inds))
% 
% end

I_up = ones(2*64, 2*64)
I_up(1:2:end, 1:2:end) = I;
%%
imshow(I_up, [0, 1], 'Parent', ax1)
axis('on')
%%

fig_path = fullfile(getfigroot, 'mu_path_mask.eps');
saveEps(F3, fig_path)
% export_fig(F3, fig_path, '-q101')


    

